"""3D cavity switching by a laser-timed, collisional ionized-gas sheet.

A short laser pulse enters through the narrow wall at x = 0 and crosses the
guide along +x in the sheet plane. Each point of the sheet starts ionizing when
the pulse reaches it; the local growth to the final density is prescribed, so
this is not a laser-breakdown or sheet-lifetime calculation. The optical field
itself is not part of the RF FDTD. No PEC mask is added at the sheet.
Atomic-physics demonstrations are separate in helium_cr.py.
"""

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import time

import numba
import numpy as np
from scipy.constants import elementary_charge as QE, electron_mass as ME, epsilon_0 as EPS0, k as KB

from .cold_rf import available_memory_bytes, estimated_memory_bytes, memory_budget_bytes
from .coupled3d import RESULTS, REF
from .drude import relative_permittivity, te10_slab_scattering
from .fdtd3d import FDTD, _plane_projection
from .geometry3d import TE10Port, iris_cavity
from .swarm import SWARM_DIR, read_table


@dataclass
class SheetCase:
    pressure_pa: float = 1e6
    gas_temperature_k: float = 300.0
    electron_temperature_ev: float = 1.0
    density_m3: float = 1e22
    thickness_m: float = 6e-3
    field_peak_v_m: float = 1e5
    trigger_s: float = 1.5e-9  # laser pulse peak at the entry wall x = 0
    rise_s: float = 0.6e-9
    laser_fwhm_s: float = 30e-12
    end_s: float = 10e-9
    frame_steps: int = 12


C0 = 299792458.0


def sheet_fraction(t, trigger, rise):
    """Finite smooth turn-on; held density is prescribed, not an afterglow."""
    x = np.clip((np.asarray(t) - trigger) / rise, 0.0, 1.0)
    return 0.5 * (1.0 - np.cos(math.pi * x))


def laser_front_fraction(t, x, trigger, rise):
    """Local turn-on that starts when the laser pulse peak reaches x (pulse moves along +x at c)."""
    return sheet_fraction(t, trigger + np.asarray(x) / C0, rise)


def laser_intensity(t, x, trigger, fwhm):
    """Normalized Gaussian pulse envelope travelling along +x at c (visual/timing only)."""
    tau = (np.asarray(t) - trigger - np.asarray(x) / C0) / fwhm
    return np.exp(-4 * math.log(2) * tau**2)


def collision_frequency(case):
    tab = read_table(SWARM_DIR / "He_maxwell.csv")
    eps = 1.5 * case.electron_temperature_ev
    if not tab["mean_energy_eV"][0] <= eps <= tab["mean_energy_eV"][-1]:
        raise ValueError("prescribed electron temperature is outside the helium swarm table")
    km = np.interp(eps, tab["mean_energy_eV"], tab["km_m3_s"])
    return float(case.pressure_pa / (KB * case.gas_temperature_k) * km)


def analytic_scan(frequency, width, thickness, nu):
    rows = []
    for ne in np.geomspace(1e17, 1e24, 150):
        r, t = te10_slab_scattering(frequency, width, thickness, float(ne), nu)
        rows.append((ne, abs(r)**2, abs(t)**2, 1 - abs(r)**2 - abs(t)**2))
    return np.asarray(rows)


def run(output: Path, case=SheetCase(), threads=2):
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if case.pressure_pa <= 8e5 or case.density_m3 <= 0 or case.rise_s <= 0 or case.laser_fwhm_s <= 0:
        raise ValueError("this demonstration requires p > 800 kPa and positive density/rise/pulse time")
    if case.end_s <= case.trigger_s + 0.1 / C0 + case.rise_s or case.frame_steps < 1:
        raise ValueError("duration must extend beyond the trigger and frame_steps must be positive")
    numba.set_num_threads(threads)
    cal = json.loads((RESULTS / "calibration.json").read_text(encoding="utf-8"))["cavity_dx1mm"]
    dx = 1e-3
    rect = REF["rectangular"]
    lay = iris_cavity(dx, rect["width_m"], rect["height_m"], rect["closed_cavity_length_m"],
                      cal["chosen"]["aperture_m"])
    if estimated_memory_bytes(lay.node_inside.size) > memory_budget_bytes(available_memory_bytes()):
        raise MemoryError("estimated 3D cavity memory exceeds the shared-host budget")
    g = lay.grid
    f0 = cal["prefill"]["f_drive_hz"]
    nu = collision_frequency(case)
    nk = int(round(case.thickness_m / dx))
    if nk < 4:
        raise ValueError("the sheet must occupy at least four 1 mm cells")
    # nk full-density nodes have integrated thickness nk*dx after edge averaging.
    k0 = lay.focus_k - nk // 2
    k1 = k0 + nk
    if not lay.gas_k_lo < k0 < k1 < lay.gas_k_hi:
        raise ValueError("sheet must lie inside the cavity")
    actual_thickness = nk * dx
    z_sheet = (k0 + k1 - 1) * dx / 2
    scan = analytic_scan(f0, lay.feed_width_m, actual_thickness, nu)
    r, tr = te10_slab_scattering(f0, lay.feed_width_m, actual_thickness, case.density_m3, nu)
    eps = relative_permittivity(case.density_m3, nu, f0)
    beta = np.sqrt((2 * math.pi * f0 / 299792458.0)**2 * eps - (math.pi / lay.feed_width_m)**2)
    skin = 1 / beta.imag
    cold = np.load(RESULTS / cal["prefill"]["file"])
    output.mkdir(parents=True)
    np.save(output / "analytic_density_scan.npy", scan)
    frame_step = case.frame_steps
    steps = int(math.ceil(case.end_s / g.dt / frame_step)) * frame_step
    frame_count = steps // frame_step + 1
    histories = {}
    volume_shape = None
    x_nodes = np.arange(g.nx + 1) * dx
    start = time.perf_counter()
    for mode in ("control", "sheet"):
        fd = FDTD(g, lay.node_inside, lay.pml_lo, lay.pml_hi, f0)
        scale = case.field_peak_v_m / cal["prefill"]["focus_peak_per_unit"]
        fd.load_state({key: int(cold[key]) if key == "step_count" else cold[key] * scale
                       for key in fd.state()})
        port = TE10Port(lay, f0)
        port.amplitude = scale
        wp2 = np.zeros(lay.node_inside.shape)
        nun = np.zeros_like(wp2)
        nun[:, :, k0:k1] = nu
        work = np.zeros_like(wp2)
        exz = np.empty((frame_count, g.nx + 1, g.nz + 1), np.float32)
        eyz = np.empty((frame_count, g.ny, g.nz + 1), np.float32)
        evo = np.empty((frame_count,) + fd.Ey[::3, ::3, ::3].shape, np.float32)
        times = np.arange(frame_count) * frame_step * g.dt
        us = np.empty(frame_count)
        ud = np.empty(frame_count)
        heat = np.zeros(frame_count)
        v = np.empty(steps)
        cur = np.empty(steps)
        for frame in range(frame_count):
            if frame:
                # Resolve the source on every RF step, not just at animation frames.
                for sub in range(frame_step):
                    step = (frame - 1) * frame_step + sub
                    if mode == "sheet":
                        frac = laser_front_fraction(step * g.dt, x_nodes, case.trigger_s, case.rise_s)
                        wp2[:, :, k0:k1] = (case.density_m3 * QE**2 / (EPS0 * ME)) * frac[:, None, None]
                        fd.set_bulk_plasma(wp2, nun)
                    fd.step(port.inject)
                    v[step], cur[step] = _plane_projection(
                        fd.Ey, fd.Hx, lay.meas_k, port.profile, port.profile)
                fd.collect_bulk_work(work)
                heat[frame] = heat[frame - 1] + work.sum() * dx**3
            exz[frame] = fd.Ey[:, g.ny // 2, :]
            eyz[frame] = fd.Ey[lay.focus_i, :, :]
            evo[frame] = fd.Ey[::3, ::3, ::3]
            us[frame] = fd.energy(lay.energy_k_lo, lay.energy_k_hi)
            ud[frame] = fd.energy(k1 + 2, lay.energy_k_hi)
        forward, backward = port.modal_amplitudes(v, cur, g.dt)
        t_modal = (np.arange(len(forward)) + 0.5 * (steps - len(forward))) * g.dt
        np.savez_compressed(output / f"{mode}_fields.npz", t=times, ey_xz=exz, ey_yz=eyz,
                            ey_volume=evo, final_Ey=fd.Ey.astype(np.float32), dx=dx)
        histories[mode] = dict(t=times, energy=us, downstream=ud, plasma_work=heat,
                               t_modal=t_modal, forward=forward, backward=backward)
        volume_shape = list(evo.shape[1:])
        print(f"{mode}: {frame_count} frames from {steps} 3D RF steps", flush=True)
    t = histories["sheet"]["t"]
    period = max(1, round(1 / f0 / (frame_step * g.dt)))
    tail = slice(-period, None)
    ratio = float(histories["sheet"]["downstream"][tail].mean()
                  / histories["control"]["downstream"][tail].mean())
    records = {f"{mode}_{key}": val for mode, series in histories.items() for key, val in series.items()}
    np.savez_compressed(output / "histories.npz", **records)
    min_ion_energy = case.density_m3 * lay.feed_width_m * lay.feed_height_m * actual_thickness * 24.587389011 * QE
    summary = dict(
        status="executed_laser_timed_prescribed_density_3d_em_comparison", case=asdict(case),
        laser=dict(entry="narrow wall x = 0 (optical port not modeled; wall remains PEC)",
                   direction="+x across the guide, in the sheet plane",
                   crossing_time_s=lay.feed_width_m / C0, pulse_fwhm_s=case.laser_fwhm_s,
                   role="sets local ionization onset time only; optical field not in RF FDTD"),
        grid=dict(nx=g.nx, ny=g.ny, nz=g.nz, dx_m=dx, dt_s=g.dt),
        frequency_hz=f0, collision_frequency_s=nu, sheet_node_k=[k0, k1],
        sheet_center_m=z_sheet, sheet_thickness_m=actual_thickness,
        cavity_entrance_m=(lay.gas_k_lo - 1) * dx, cavity_end_m=g.nz * dx,
        width_m=lay.feed_width_m, height_m=lay.feed_height_m,
        measurement_z_m=lay.meas_k * dx, volume_shape=volume_shape,
        analytic_isolated_slab=dict(R=abs(r)**2, T=abs(tr)**2,
                                    A=1-abs(r)**2-abs(tr)**2, eps_real=eps.real,
                                    eps_imag=eps.imag, attenuation_length_m=skin),
        downstream_energy_ratio_last_period=ratio, minimum_sheet_ionization_energy_J=min_ion_energy,
        initial_stored_energy_J=float(histories["control"]["energy"][0]),
        elapsed_s=time.perf_counter() - start,
        limitations=[
            "Laser-timed prescribed density: ionization onset follows the pulse arrival at each x; growth to the final density (rise time) and 1 eV Maxwellian collision closure are prescribed, not calculated breakdown or lifetime.",
            "No new metal/PEC is inserted: the sheet consists only of Drude electrons in helium.",
            "Primary EM animation does not couple the separately demonstrated collisional-radiative network.",
            "1 mm mesh: finite skin depth is only marginally resolved; cavity result is qualitative.",
            "Isolated analytic slab R,T,A are not the reflection/transmission of this shorted cavity.",
            "Transient port backward/forward power can exceed one because the cavity releases stored energy.",
            "RF walls are PEC; neutral-gas permittivity, material dispersion and wall losses are omitted.",
        ],
        collision_source="cases/swarm/He_maxwell.csv, conductivity-equivalent momentum rate from IST-Lisbon/LXCat via BOLOS",
    )
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--threads", type=int, default=2)
    args = p.parse_args()
    print(json.dumps(run(args.output.resolve(), threads=args.threads), indent=2))


if __name__ == "__main__":
    main()
