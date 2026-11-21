"""Bounded transient Maxwell/CR run with an explicitly perpendicular optical seed.

The narrow laser source is number-conservingly averaged onto bulk nodes, not
represented as a wire parallel to Ey. The optical field is vacuum Gaussian;
its bulk-grid absorption is diagnosed but physical subgrid attenuation and
refraction are not validated by that diagnostic.
Gas motion is handled in a separate local, RF-off afterglow handoff.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numba
import numpy as np

from . import plasma3d as k3
from . import plasma3d_cr as kc
from .breakdown_probability import Setup, seed_per_length
from .coupled3d import RESULTS, _focus_edge, _load_cal, _poynting_z, make_layout
from .fdtd3d import FDTD
from .gas_model import load_gas_model
from .gas_dynamics import chemical_energy
from .geometry3d import TE10Port
from .radiation_diagnostics import cr_emission, free_free_estimate


def perpendicular_seed(active, focus, dx, beam):
    """Return node density and total count for an x-directed beam, perpendicular to Ey."""
    fi, fj, fk = focus
    x = (np.arange(active.shape[0]) - fi) * dx
    count = seed_per_length(beam, beam.energy_j, x) * dx
    density = np.zeros(active.shape)
    for j, weight in ((fj, 0.5), (fj + 1, 0.5)):
        density[:, j, fk] = count * weight / dx ** 3
    density[~active] = 0
    return density, float(density.sum() * dx ** 3)


def coarsen_plane(values):
    indices = [np.arange(0, n, 2) for n in values.shape]
    counts = [np.minimum(2, n - i) for n, i in zip(values.shape, indices)]
    summed = np.add.reduceat(np.add.reduceat(values, indices[0], axis=0), indices[1], axis=1)
    return summed / (counts[0][:, None] * counts[1][None, :])


@numba.njit(cache=True)
def weighted_edge_work(px, py, pz, wp2, work):
    nx, ny, nz = wp2.shape
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                weight = wp2[i, j, k]
                value = 0.0
                if weight > 0:
                    if i > 0:
                        value += px[i-1, j, k] * weight / (weight + wp2[i-1, j, k])
                    if i < nx-1:
                        value += px[i, j, k] * weight / (weight + wp2[i+1, j, k])
                    if j > 0:
                        value += py[i, j-1, k] * weight / (weight + wp2[i, j-1, k])
                    if j < ny-1:
                        value += py[i, j, k] * weight / (weight + wp2[i, j+1, k])
                    if k > 0:
                        value += pz[i, j, k-1] * weight / (weight + wp2[i, j, k-1])
                    if k < nz-1:
                        value += pz[i, j, k] * weight / (weight + wp2[i, j, k+1])
                else:
                    if i < nx-1 and wp2[i+1, j, k] == 0 and px[i, j, k] != 0:
                        raise ArithmeticError("Drude-memory work has no electron population")
                    if j < ny-1 and wp2[i, j+1, k] == 0 and py[i, j, k] != 0:
                        raise ArithmeticError("Drude-memory work has no electron population")
                    if k < nz-1 and wp2[i, j, k+1] == 0 and pz[i, j, k] != 0:
                        raise ArithmeticError("Drude-memory work has no electron population")
                work[i, j, k] = value


@numba.njit(cache=True)
def advance_cells(coords, state, tg, gas_energy, work, dt, n0, prm, tabs, ln_e0, h,
                  ln_i0, hi, ppt_tab, wp2, nu, ledger):
    for index in range(len(coords)):
        i, j, k = coords[index]
        cell_state = state[:, i, j, k].copy()
        neutral = n0 - cell_state[0] - 2 * cell_state[1] - cell_state[3] - 2 * cell_state[4] - cell_state[5]
        if neutral < 0:
            raise ArithmeticError("cavity species exceed helium nuclei inventory")
        out = np.zeros(11)
        kc.kinetic_cell_cr(cell_state, dt, neutral, tg[i, j, k], work[i, j, k] / (k3.QE * dt),
                          0.0, prm, tabs, ln_e0, h, ln_i0, hi, ppt_tab, out)
        gas_energy[i, j, k] += out[6]
        particles = n0 - cell_state[1] - cell_state[4]
        if particles <= 0 or gas_energy[i, j, k] <= 0:
            raise ArithmeticError("invalid cavity heavy-gas energy or particle inventory")
        tg[i, j, k] = gas_energy[i, j, k] / (1.5 * k3.KB * particles)
        state[:, i, j, k] = cell_state
        neutral = n0 - cell_state[0] - 2 * cell_state[1] - cell_state[3] - 2 * cell_state[4] - cell_state[5]
        if neutral < 0:
            raise ArithmeticError("CR update exceeds cavity nuclei inventory")
        ne = cell_state[0] + cell_state[1] + 2 * cell_state[5]
        wp2[i, j, k], nu[i, j, k] = k3.drude_params(ne, cell_state[2], neutral, tabs, ln_e0, h, prm)
        ledger[i] += out


def run(out: Path, field_kv_cm=100.0, duration_ns=18.0, frames_count=180, threads=2,
        block=4, laser=True):
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    if not 0 < field_kv_cm <= 200 or not 2 < duration_ns <= 40 or frames_count < 20:
        raise ValueError("bounded run requires 0<E<=200 kV/cm, 2<t<=40 ns and at least 20 frames")
    if not isinstance(block, int) or not 1 <= block <= 16:
        raise ValueError("coupling block must be an integer from 1 to 16")
    numba.set_num_threads(threads)
    cal = _load_cal()["cavity_dx1mm"]
    dx = 1e-3
    lay = make_layout("cavity", dx, cal["chosen"]["aperture_m"])
    fd = FDTD(lay.grid, lay.node_inside, lay.pml_lo, lay.pml_hi, cal["prefill"]["f_drive_hz"])
    scale = field_kv_cm * 1e5 / cal["prefill"]["focus_peak_per_unit"]
    prefill = dict(np.load(RESULTS / cal["prefill"]["file"]))
    prefill.pop("ey_focus")
    for name in prefill:
        if name != "step_count":
            prefill[name] = prefill[name] * scale
    fd.load_state(prefill)
    loaded_step_count = fd.step_count
    port = TE10Port(lay, cal["prefill"]["f_drive_hz"])
    port.amplitude = scale
    port.envelope = lambda t: 1.0
    focus = _focus_edge(lay)
    fi, fj, fk = focus
    active = lay.node_inside.copy()
    active[:, :, :lay.gas_k_lo] = False
    active[:, :, lay.gas_k_hi + 1:] = False
    beam = Setup(axis="x", energy_j=0.75e-3, n_ref_m3=1e17, n_sat_m3=1e16)
    n0 = beam.n_gas
    gm = load_gas_model("He", beam.wavelength_m, allow_unverified=True, extended=True, cr=True)
    source, count = perpendicular_seed(active, focus, dx, beam)
    if not laser:
        source[:] = 0
        count = 0.0
    kinetic_active = active & (source > 0)
    coords = np.argwhere(kinetic_active)
    seed_ev = 0.5
    state = np.zeros((6,) + active.shape)
    tg = np.full(active.shape, beam.t_k)
    gas_energy = np.full(active.shape, 1.5 * k3.KB * n0 * beam.t_k)
    work = np.zeros(active.shape)
    wp2, nu = np.zeros(active.shape), np.zeros(active.shape)
    ledger = np.zeros((active.shape[0], 11))
    dt = lay.grid.dt
    dtb = dt * block
    steps = math.ceil(duration_ns * 1e-9 / dtb)
    stride = max(1, steps // frames_count)
    initial_em = fd.energy(lay.energy_k_lo, lay.energy_k_hi)
    flux_in = 0.0
    previous_fraction = 0.0
    delivered, seed_input = 0.0, 0.0
    optical_depth_max = 0.0
    nu_optical_max = 0.0
    plasma_optical_max = 0.0
    edge_work = 0.0
    mapped_work = 0.0
    probe_t, probe_ey = [], []
    data = {k: [] for k in ("time_s", "ey_xz", "ne_xz", "tg_xz", "uv_xz", "visible_xz",
                            "field_focus_v_m", "electrons", "em_energy_j", "gas_heat_j",
                            "sink_j", "seed_electrons", "laser_power_w", "optical_depth")}
    handoff = None
    t_laser = 2e-9
    laser_sigma = beam.fwhm_s / math.sqrt(8 * math.log(2))
    peak_power = beam.energy_j / (laser_sigma * math.sqrt(2 * math.pi))
    for b in range(steps):
        for substep in range(block):
            fd.step(port.inject)
            probe_t.append((b * block + substep + 1) * dt)
            probe_ey.append(float(fd.Ey[fi, fj, fk]))
            flux_in += _poynting_z(fd.Ex, fd.Ey, fd.Hx, fd.Hy, lay.meas_k) * dx ** 2 * dt
        t = (b + 1) * dtb
        fraction = 0.5 * (1 + math.erf((t - t_laser) / (math.sqrt(2) * laser_sigma)))
        delta = fraction - previous_fraction
        state[0] += source * delta
        state[2] += source * delta * seed_ev
        previous_fraction = fraction
        delivered = beam.energy_j * fraction if laser else 0.0
        seed_input = count * fraction * (gm.prm[k3.P_IP] + seed_ev) * k3.QE
        if seed_input > delivered:
            raise ArithmeticError("prescribed seed costs more energy than the incident laser")
        total_edge = float(fd.Px.sum() + fd.Py.sum() + fd.Pz.sum())
        scale_edge = float(np.abs(fd.Px).sum() + np.abs(fd.Py).sum() + np.abs(fd.Pz).sum())
        weighted_edge_work(fd.Px, fd.Py, fd.Pz, wp2, work)
        total_mapped = float(work.sum())
        if abs(total_edge - total_mapped) > 1e-12 * scale_edge + 1e-30:
            raise ArithmeticError("edge RF work was not conservatively mapped into plasma populations")
        edge_work += total_edge * dx ** 3
        mapped_work += total_mapped * dx ** 3
        fd.Px[:] = 0
        fd.Py[:] = 0
        fd.Pz[:] = 0
        advance_cells(coords, state, tg, gas_energy, work, dtb, n0, gm.prm, gm.tabs,
                      gm.ln_e0, gm.h, gm.ln_i0, gm.hi, gm.ppt_tab, wp2, nu, ledger)
        fd.set_bulk_plasma(wp2, nu)
        ne = state[0] + state[1] + 2 * state[5]
        opacity = ne[1:-1, fj, fk] * k3.QE ** 2 * nu[1:-1, fj, fk] / (
            k3.ME * k3.EPS0 * k3.C0 * (2 * math.pi * k3.C0 / beam.wavelength_m) ** 2)
        optical_depth = float(opacity.sum() * dx)
        power = peak_power * math.exp(-0.5 * ((t - t_laser) / laser_sigma) ** 2) if laser else 0.0
        if power > 0.01 * peak_power:
            optical_depth_max = max(optical_depth_max, optical_depth)
            omega_optical = 2 * math.pi * k3.C0 / beam.wavelength_m
            nu_optical_max = max(nu_optical_max, float(nu.max()) / omega_optical)
            plasma_optical_max = max(plasma_optical_max, float(wp2.max()) / omega_optical ** 2)
        if b % stride == 0 or b == steps - 1:
            resonance, excimer = cr_emission(state, gm.prm)
            _, visible = free_free_estimate(state)
            bl = ledger.sum(axis=0) * dx ** 3
            data["time_s"].append(t)
            for key, values in (("ey_xz", fd.Ey[:, fj, :]), ("ne_xz", ne[:, fj, :]),
                                ("tg_xz", tg[:, fj, :]), ("uv_xz", (resonance + excimer)[:, fj, :]),
                                ("visible_xz", visible[:, fj, :])):
                data[key].append(coarsen_plane(values).astype(np.float32))
            for key, value in (("field_focus_v_m", fd.Ey[fi, fj, fk]),
                               ("electrons", ne[active].sum() * dx ** 3),
                               ("em_energy_j", fd.energy(lay.energy_k_lo, lay.energy_k_hi)),
                               ("gas_heat_j", bl[6]), ("sink_j", bl[7]),
                               ("seed_electrons", count * fraction), ("laser_power_w", power),
                               ("optical_depth", optical_depth)):
                data[key].append(float(value))
            fraction_excited = (state[0] + 2 * state[1] + state[3] + 2 * state[4] + state[5]) / n0
            if t > t_laser + 3 * beam.fwhm_s and tg[active].max() < 3000 and fraction_excited.max() < 0.01:
                handoff = dict(time_s=t, state_yz=state[:, fi].copy(), tg_yz=tg[fi].copy(),
                               gas_energy_yz=gas_energy[fi].copy())
            print(f"E={field_kv_cm:g} kV/cm t={t*1e9:.2f} ns N={data['electrons'][-1]:.3g}", flush=True)
    if handoff is None:
        raise ArithmeticError("no snapshot meets the declared low-ionization afterglow EOS limits")
    out.mkdir(parents=True)
    map_x = np.arange(0, active.shape[0], 2)
    map_z = np.arange(0, active.shape[2], 2)
    np.savez_compressed(out / "frames.npz", **{k: np.asarray(v) for k, v in data.items()},
                        probe_time_s=np.asarray(probe_t), probe_ey_v_m=np.asarray(probe_ey),
                        x_m=(map_x+np.minimum(map_x+1,active.shape[0]-1))*0.5*dx,
                        z_m=(map_z+np.minimum(map_z+1,active.shape[2]-1))*0.5*dx,
                        active_xz=coarsen_plane(active[:, fj, :].astype(float)))
    np.savez_compressed(out / "handoff.npz", **handoff, n0=n0, dx=dx, focus=np.asarray(focus))
    bl = ledger.sum(axis=0) * dx ** 3
    final_electron = state[2].sum() * k3.QE * dx ** 3
    final_chem = chemical_energy(state, gm.prm).sum() * dx ** 3
    inputs = bl[0] + seed_input
    closure = inputs - final_electron - final_chem - bl[6] - bl[7]
    final_em = fd.energy(lay.energy_k_lo, lay.energy_k_hi)
    rf_residual = flux_in - final_em + initial_em - edge_work
    summary = dict(field_peak_kv_cm=field_kv_cm, field_target_convention="empty-cavity focus peak",
                   pressure_pa=beam.p_pa, temperature_k=beam.t_k,
                   frequency_hz=port.f, dx_m=dx, time_step_s=dt, coupling_step_s=dtb,
                   loaded_prefill_step_count=loaded_step_count, rf_clock_offset_s=loaded_step_count*dt,
                   laser_rf_phase_rad=float((2*math.pi*port.f*(loaded_step_count*dt+t_laser)) % (2*math.pi)),
                   rf_source_convention="sin(2*pi*f*absolute_RF_clock)", laser_enabled=laser,
                   laser_axis="x perpendicular to Ey", laser_energy_j=beam.energy_j,
                   laser_fwhm_s=beam.fwhm_s, laser_waist_m=beam.w0_m, laser_time_s=t_laser,
                   seed_model="empirical saturated I^6 active-site kernel; not measured",
                   seed_time_profile="empirical Gaussian-CDF injection with laser FWHM, not a multiphoton rate",
                   seed_unsaturated_reference_amplitude_m3=beam.n_ref_m3,
                   seed_saturated_reference_density_m3=float(beam.n_sat_m3 * -np.expm1(-beam.n_ref_m3/beam.n_sat_m3)),
                   seed_expected_electrons=count, seed_input_j=seed_input, laser_delivered_j=delivered,
                   optical_depth_during_pulse=optical_depth_max,
                   optical_depth_scope="bulk-grid IB diagnostic; physical beam attenuation/refraction unverified",
                   bulk_grid_nu_over_optical_omega_max=nu_optical_max,
                   bulk_grid_wp2_over_optical_omega2_max=plasma_optical_max,
                   plasma_energy_residual_rel=float(closure / max(inputs, final_chem, 1e-30)),
                   edge_to_node_work_residual_j=float(edge_work-mapped_work),
                   edge_to_kinetic_work_residual_j=float(edge_work-bl[0]),
                   rf_energy_residual_j=float(rf_residual),
                   rf_energy_residual_rel=float(rf_residual/max(initial_em,abs(flux_in),abs(edge_work),1e-30)),
                   rf_residual_relative_plasma_work=float(rf_residual/abs(edge_work)) if edge_work else None,
                   rf_boundary_flux_j=flux_in, initial_em_energy_j=initial_em, final_em_energy_j=final_em,
                   rf_work_interpretation="net J dot E work, not automatically irreversible absorption",
                   rf_net_work_j=float(bl[0]),
                   rf_absorbed_j=float(bl[0]),
                   rf_absorbed_j_note="compatibility alias for rf_net_work_j, not irreversible absorption",
                   gas_heat_j=float(bl[6]), escaped_sink_j=float(bl[7]),
                   modeled_radiative_energy_j=float(bl[9]+bl[10]),
                   unresolved_sink_energy_j=float(bl[7]-bl[9]-bl[10]),
                   handoff_time_s=handoff["time_s"], provenance=gm.provenance,
                   limitations=["bulk-node-averaged optical seed: sub-mm plasma structure unresolved",
                                "gas is isochoric during this short 3D RF stage",
                                "radiation uses existing lumped CR escape rates; visible free-free is conditional",
                                "RF energy/flux diagnostic uses staggered fields and is not an exact discrete identity",
                                "source increment is inserted at block end before a full kinetic block; coupling is first order",
                                "no diffusion in this restricted 3D demonstration",
                                "RF remains on; separate handoff assumes RF quench, not the continuing experiment"])
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--field-kv-cm", type=float, default=100)
    p.add_argument("--duration-ns", type=float, default=18)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--block", type=int, default=4)
    p.add_argument("--no-laser", action="store_true")
    a = p.parse_args()
    print(json.dumps(run(a.output, a.field_kv_cm, a.duration_ns, threads=a.threads,
                         block=a.block, laser=not a.no_laser), indent=2))
