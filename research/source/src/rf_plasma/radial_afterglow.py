"""Local cylindrical RF-off afterglow: Euler expansion coupled to existing CR chemistry.

The handoff conserves integrated electron, excitation and excess heavy-gas
energy per unit beam length. It is an axisymmetric projection of a coarse
3D slice, not a resolved optical filament or continuation of the driven cavity.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from numba import njit

from . import plasma3d as k3
from . import plasma3d_cr as kc
from .gas_dynamics import Euler1D, M_HE, chemical_energy
from .gas_model import load_gas_model
from .radiation_diagnostics import cr_emission, free_free_estimate


def project_slice(source, nr, rmax, prm):
    dx = float(source["dx"])
    n0 = float(source["n0"])
    _, fj, fk = source["focus"]
    state = source["state_yz"]
    yy = (np.arange(state.shape[1]) - fj - 0.5) * dx
    zz = (np.arange(state.shape[2]) - fk) * dx
    radius = np.hypot(yy[:, None], zz[None, :])
    coarse = np.arange(0, rmax + 0.5 * dx, dx)
    if coarse[-1] < rmax:
        coarse = np.append(coarse, rmax)
    area = np.pi * np.diff(coarse ** 2)
    profiles = []
    for field in state:
        hist = np.histogram(radius, bins=coarse, weights=field * dx ** 2)[0]
        profiles.append(hist / area)
    excess = source["gas_energy_yz"] - 1.5 * k3.KB * n0 * 300.0
    heat = np.histogram(radius, bins=coarse, weights=excess * dx ** 2)[0] / area
    solver = Euler1D(np.linspace(0, rmax, nr + 1), cylindrical=True, plasma_prm=prm)
    # Exact annular overlap, rather than interpolated point values.
    inner = np.maximum(solver.edges[:-1, None], coarse[None, :-1])
    outer = np.minimum(solver.edges[1:, None], coarse[None, 1:])
    overlap = np.pi * np.maximum(outer ** 2 - inner ** 2, 0)
    u = np.zeros((9, nr))
    u[0] = n0 * M_HE
    u[3:] = np.asarray(profiles) @ overlap.T / solver.volumes
    gas = 1.5 * k3.KB * n0 * 300.0 + overlap @ heat / solver.volumes
    u[2] = gas + u[5] * k3.QE + chemical_energy(u[3:], prm)
    original = np.asarray(profiles) @ area
    mapped = (u[3:] * solver.volumes).sum(axis=1)
    np.testing.assert_allclose(mapped, original, rtol=2e-12, atol=1e-20)
    np.testing.assert_allclose(((gas - 1.5 * k3.KB * n0 * 300) * solver.volumes).sum(),
                               float(heat @ area), rtol=2e-12, atol=1e-12)
    solver.primitive(u)
    return solver, u


@njit(cache=True)
def chemistry_step(u, dt, prm, tabs, ln_e0, h, ln_i0, hi, ppt_tab):
    escaped = np.zeros(u.shape[1])
    gas_heat = np.zeros(u.shape[1])
    photons = np.zeros(u.shape[1])
    for i in range(u.shape[1]):
        n0 = u[0, i] / M_HE
        state = u[3:, i].copy()
        chem = k3.QE * (prm[k3.P_IP] * state[0] + (prm[k3.P_IP] - prm[kc.P_DION]) * state[1]
                       + prm[kc.P_ESTAR] * state[3] + prm[kc.P_EX] * state[4]
                       + (prm[k3.P_IP] + prm[kc.P_I2]) * state[5])
        gas = u[2, i] - 0.5 * u[1, i] ** 2 / u[0, i] - state[2] * k3.QE - chem
        particles = n0 - state[1] - state[4]
        neutral = n0 - state[0] - 2 * state[1] - state[3] - 2 * state[4] - state[5]
        if gas <= 0 or particles <= 0 or neutral < 0:
            raise ArithmeticError("invalid afterglow gas energy or heavy-particle density")
        tg = gas / (1.5 * k3.KB * particles)
        out = np.zeros(11)
        kc.kinetic_cell_cr(state, dt, neutral, tg, 0.0, 0.0, prm, tabs, ln_e0, h,
                          ln_i0, hi, ppt_tab, out, 0.0)
        u[3:, i] = state
        u[2, i] -= out[7]
        escaped[i] = out[7]
        gas_heat[i] = out[6]
        photons[i] = out[9] + out[10]
    return escaped, gas_heat, photons


def run(source_path: Path, out: Path, nr=160, duration_s=1e-6, samples=81, trap_radius_m=5e-3,
        source_step_s=5e-9):
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    if nr < 40 or nr > 640 or not 0 < duration_s <= 2e-6 or samples < 2:
        raise ValueError("invalid bounded afterglow settings")
    if not 0 < source_step_s <= 5e-9:
        raise ValueError("chemistry splitting step must lie between zero and 5 ns")
    with np.load(source_path) as data:
        source = {k: data[k] for k in data.files}
    gm = load_gas_model("He", 532e-9, allow_unverified=True, extended=True, cr=True)
    if not 0 < trap_radius_m <= 0.1:
        raise ValueError("invalid fixed trapping radius")
    old_escape = gm.optics["resonance_escape_factor"]
    new_escape = min(1.0, old_escape * np.sqrt(17e-3 / trap_radius_m))
    gm.prm[kc.P_ARES] *= new_escape / old_escape
    solver, u = project_slice(source, nr, 5e-3, gm.prm)
    initial = (u[:3] * solver.volumes).sum(axis=1)
    baseline_gas = 1.5 * k3.KB * float(source["n0"]) * 300 * solver.volumes.sum()
    initial_excess = initial[2] - baseline_gas
    times = np.linspace(0, duration_s, samples)
    t, escaped, heat, photon_energy = 0.0, 0.0, 0.0, 0.0
    frames = []
    step_count = 0
    for target in times:
        while t < target * (1 - 1e-12):
            dt = min(target - t, solver.stable_dt(u), source_step_s)
            loss, gain, photons = chemistry_step(u, 0.5 * dt, gm.prm, gm.tabs, gm.ln_e0, gm.h,
                                        gm.ln_i0, gm.hi, gm.ppt_tab)
            escaped += float(loss @ solver.volumes)
            heat += float(gain @ solver.volumes)
            photon_energy += float(photons @ solver.volumes)
            u = solver.step(u, dt)
            loss, gain, photons = chemistry_step(u, 0.5 * dt, gm.prm, gm.tabs, gm.ln_e0, gm.h,
                                        gm.ln_i0, gm.hi, gm.ppt_tab)
            escaped += float(loss @ solver.volumes)
            heat += float(gain @ solver.volumes)
            photon_energy += float(photons @ solver.volumes)
            solver.primitive(u)
            t += dt
            step_count += 1
        rho, vel, pressure, gas = solver.primitive(u)
        particles = rho / M_HE - u[4] - u[7]
        if np.any(particles <= 0):
            raise ArithmeticError("nonpositive heavy-particle count in diagnostic")
        tg = gas / (1.5 * k3.KB * particles)
        ne = u[3] + u[4] + 2 * u[8]
        resonance, excimer = cr_emission(u[3:], gm.prm)
        _, visible = free_free_estimate(u[3:])
        frames.append(dict(rho=rho.copy(), velocity=vel.copy(), pressure=pressure.copy(),
                           gas_temperature=tg, electron_density=ne, uv_power_density=resonance + excimer,
                           visible_continuum_estimate=visible, sink_energy_j_per_m=escaped,
                           radiative_energy_j_per_m=photon_energy))
    final = (u[:3] * solver.volumes).sum(axis=1)
    mass_residual = (final[0] + solver.boundary_out[0] - initial[0]) / initial[0]
    energy_error = final[2] + solver.boundary_out[2] + escaped - initial[2]
    signal_energy = max(abs(initial_excess), escaped, 1e-12 * initial[2])
    energy_residual = energy_error / signal_energy
    roundoff_bound = 32 * np.finfo(float).eps * max(step_count, 1) * abs(initial[2])
    if abs(mass_residual) > 1e-10 or abs(energy_error) > 1e-8 * signal_energy + roundoff_bound:
        raise ArithmeticError(f"afterglow conservation failed: mass={mass_residual}, energy={energy_residual}")
    packed = {key: np.asarray([frame[key] for frame in frames]) for key in frames[0]}
    if any(not np.isfinite(values).all() for values in packed.values()):
        raise ArithmeticError("non-finite afterglow diagnostic")
    out.mkdir(parents=True)
    np.savez_compressed(out / "afterglow.npz", time_s=times, r_m=solver.centres, **packed)
    summary = dict(source=str(source_path), handoff_time_s=float(source["time_s"]),
                   nr=nr, duration_s=duration_s, maximum_source_step_s=source_step_s,
                   geometry="cylindrical, per unit laser-axis length",
                   boundaries="reflecting axis, outflow outer boundary at 5 mm",
                   mass_residual_rel=float(mass_residual), energy_residual_rel=float(energy_residual),
                   energy_residual_j_per_m=float(energy_error), roundoff_bound_j_per_m=float(roundoff_bound),
                   plasma_energy_signal_resolved=bool(signal_energy > 100 * roundoff_bound),
                   step_count=step_count,
                   initial_excess_energy_j_per_m=float(initial_excess),
                   sink_energy_j_per_m=float(escaped), radiative_energy_j_per_m=float(photon_energy),
                   unresolved_sink_energy_j_per_m=float(escaped - photon_energy),
                   trapping_radius_m=trap_radius_m, resonance_escape_factor=float(new_escape),
                   chemistry_to_gas_j_per_m=float(heat),
                   maximum_gas_temperature_k=float(packed["gas_temperature"].max()),
                   maximum_velocity_m_s=float(np.abs(packed["velocity"]).max()),
                   limitations=["axisymmetric local handoff; not a 3D cavity-flow continuation",
                                "RF is assumed quenched at handoff; actual 5-us drive is not continued",
                                "1-mm source resolution remains even when Euler mesh is refined",
                                "inviscid monatomic gas closure restricted to low ionization",
                                "no drift diffusion, heat conduction, line-resolved emission or wall reflections",
                                "visible continuum is an optically-thin Maxwellian estimate, not a camera signal"])
    summary["limitations"].append("fixed-radius Holstein escape reduction; excimer radiative rate is unverified")
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nr", type=int, default=160)
    parser.add_argument("--source-step-ns", type=float, default=5)
    parser.add_argument("--trap-radius-mm", type=float, default=5)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output, nr=args.nr, source_step_s=args.source_step_ns*1e-9,
                         trap_radius_m=args.trap_radius_mm*1e-3), indent=2))
