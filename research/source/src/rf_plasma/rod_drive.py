"""Bridged plasma rod under continued RF drive in the WR-284 TE101 cavity (sealed 200 psig helium).

After rod3d the channel spans the 34 mm gap along E at the field maximum. From then on the rod is an
inductive post with a plasma resistance. This module couples

* Cavity: transmission-line view of the TE101 cavity at the post plane (two shorted stubs of length d/2)
  -> near-resonance parallel resonator with gap voltage V = E b, C = d / (2 Z0 v_g) (= eps0 a d / 4b),
  Z0 = 2 b Z_TE / a (power-consistent with V), external conductance G_ext = omega C / Q_L (iris; walls
  PEC). Slowly varying envelope: 2C dV/dt = I_d - (G_ext + Y_rod) V, drive I_d = G_ext V0 (Norton
  equivalent of the generator behind the iris; P_forward = G_ext V0^2 / 8).
* Rod: shunt post admittance Y_rod = 1 / (R + j X); X from Marcuvitz (Waveguide Handbook, sec. 5.11,
  centred inductive post) for the instantaneous conducting radius; R = b / int sigma dA with
  Re sigma = e^2 n_e nu / (m_e (nu^2 + omega^2)), nu = N k_m(eps) + nu_ei (skin depth >> radius checked).
* Plasma: radial cylindrical Euler (gas_dynamics.Euler1D, passive plasma tracers in the energy) plus the
  collisional-radiative helium cell (plasma3d_cr.kinetic_cell_cr) with the Joule power density
  p = Re(sigma) |E_ax|^2 / 2, E_ax = I_rod R / b uniform in r (tangential continuity). Strang splitting.
* Electrode sheaths, wall erosion, radial diffusion and heat conduction (<10 um per 100 ns at 1.5 MPa) and
  3D current redistribution are neglected; the CR closure is restricted to moderate ionization.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from numba import njit

from . import plasma3d as k3
from . import plasma3d_cr as kc
from .gas_dynamics import Euler1D, M_HE, chemical_energy
from .gas_model import load_gas_model

C0 = 299792458.0
EPS0 = 8.8541878128e-12
MU0 = 4e-7 * math.pi
ETA0 = math.sqrt(MU0 / EPS0)
QE = 1.602176634e-19
ME = 9.1093837015e-31
KB = 1.380649e-23
PSI = 6894.757
ATM = 101325.0
A_WG, B_WG, D_CAV = 0.072136, 0.034036, 0.06927
OUT = Path(__file__).resolve().parents[2] / "results" / "rod_drive"


@dataclass
class DriveConfig:
    name: str = "drive_E200"
    p_pa: float = 200 * PSI + ATM
    t_k: float = 300.0
    f_hz: float = 2.95e9
    q_l: float = 331.0
    e0_v_m: float = 2.0e7            # empty-cavity steady peak field (sets the forward power)
    e_start_v_m: float = 0.0         # field at bridging (0: = e0)
    drive_off_s: float = -1.0        # time after bridging when the generator stops (<0: never)
    n0_m3: float = 3e23              # rod3d / streamer2d channel density
    r0_m: float = 100e-6             # Gaussian 1/e radius of the bridged channel
    eps0_ev: float = 1.5
    r_max_m: float = 3e-3
    nr: int = 150
    t_end_s: float = 300e-9
    dt_max_s: float = 0.05e-9
    samples: int = 301

    @property
    def n_gas(self):
        return self.p_pa / (KB * self.t_k)


def cavity_constants(c: DriveConfig):
    w = 2 * math.pi * c.f_hz
    lam = C0 / c.f_hz
    ratio = math.sqrt(1 - (lam / (2 * A_WG)) ** 2)
    z_te = ETA0 / ratio
    z0 = 2 * B_WG * z_te / A_WG
    vg = C0 * ratio
    cap = D_CAV / (2 * z0 * vg)
    g_ext = w * cap / c.q_l
    lam_g = lam / ratio
    return dict(omega=w, z0=z0, cap=cap, g_ext=g_ext, lam_g=lam_g, c_eps=EPS0 * A_WG * D_CAV / (4 * B_WG))


def post_reactance(r_post, f_hz):
    """Marcuvitz centred inductive post: X/Z0 = (a/2 lam_g)[ln(4a/(pi d)) - 2 + 2 sum_{n=3,5,..}
    (1/sqrt(n^2 - (2a/lam)^2) - 1/n)], d = 2 r_post (thin post, d << a)."""
    lam = C0 / f_hz
    lam_g = lam / math.sqrt(1 - (lam / (2 * A_WG)) ** 2)
    k = 2 * A_WG / lam
    s = sum(1 / math.sqrt(n * n - k * k) - 1 / n for n in range(3, 400, 2))
    return (A_WG / (2 * lam_g)) * (math.log(4 * A_WG / (math.pi * 2 * r_post)) - 2 + 2 * s)


@njit(cache=True)
def conductivity(u, omega, prm, tabs, ln_e0, h):
    """Re sigma per cell from the CR state (electron energy from the state, consistent closure)."""
    n = u.shape[1]
    sig = np.zeros(n)
    for i in range(n):
        nxp, nx2p, w, nm, nx, na = u[3, i], u[4, i], u[5, i], u[6, i], u[7, i], u[8, i]
        ne = nxp + nx2p + 2 * na
        if ne <= 1e6:
            continue
        n0 = u[0, i] / M_HE
        neu = max(n0 - nxp - 2 * nx2p - na - nm - 2 * nx, 0.0)
        eps = max(w / ne, k3.EPS_LO)
        wmx = kc._wsel(ne, neu, eps, tabs, ln_e0, h, prm)
        km = kc._rx(tabs, 0, wmx, ln_e0, h, eps, prm)
        nu = neu * km + k3.nu_ei(ne, eps / 1.5)
        sig[i] = QE * QE * ne * nu / (ME * (nu * nu + omega * omega))
    return sig


@njit(cache=True)
def chemistry_step(u, dt, p_ext, prm, tabs, ln_e0, h, ln_i0, hi, ppt_tab):
    """Kinetics with Joule input p_ext (eV m^-3 s^-1) per cell; returns deposited, sink, radiated (J/m^3)."""
    n = u.shape[1]
    dep = np.zeros(n); sink = np.zeros(n); rad = np.zeros(n)
    for i in range(n):
        n0 = u[0, i] / M_HE
        state = u[3:, i].copy()
        chem = k3.QE * (prm[k3.P_IP] * state[0] + (prm[k3.P_IP] - prm[kc.P_DION]) * state[1]
                        + prm[kc.P_ESTAR] * state[3] + prm[kc.P_EX] * state[4]
                        + (prm[k3.P_IP] + prm[kc.P_I2]) * state[5])
        gas = u[2, i] - 0.5 * u[1, i] ** 2 / u[0, i] - state[2] * k3.QE - chem
        particles = n0 - state[1] - state[4]
        neutral = n0 - state[0] - 2 * state[1] - state[3] - 2 * state[4] - state[5]
        if gas <= 0 or particles <= 0 or neutral < 0:
            raise ArithmeticError("invalid rod gas energy or heavy-particle density")
        tg = gas / (1.5 * k3.KB * particles)
        out = np.zeros(11)
        kc.kinetic_cell_cr(state, dt, neutral, tg, p_ext[i], 0.0, prm, tabs, ln_e0, h,
                           ln_i0, hi, ppt_tab, out, 0.0)
        u[3:, i] = state
        u[2, i] += out[0] - out[7]
        dep[i] = out[0]; sink[i] = out[7]; rad[i] = out[9] + out[10]
    return dep, sink, rad


def initial_state(c: DriveConfig, prm):
    edges = c.r_max_m * np.linspace(0, 1, c.nr + 1) ** 1.6
    solver = Euler1D(edges, cylindrical=True, plasma_prm=prm)
    r = solver.centres
    n = c.n_gas
    ne = c.n0_m3 * np.exp(-(r / c.r0_m) ** 2)
    u = np.zeros((9, c.nr))
    u[0] = n * M_HE
    u[4] = ne                     # He2+ (conversion He+ + 2He -> He2+ + He is ~1e10 s^-1 at 1.5 MPa)
    u[5] = ne * c.eps0_ev
    u[2] = 1.5 * KB * n * c.t_k + u[5] * QE + chemical_energy(u[3:], prm)
    solver.primitive(u)
    return solver, u


def run(c: DriveConfig, out_root: Path = OUT, save=True, verbose=True):
    gm = load_gas_model("He", 532e-9, allow_unverified=True, extended=True, cr=True)
    prm, tabs = gm.prm, gm.tabs
    cc = cavity_constants(c)
    w, cap, g_ext = cc["omega"], cc["cap"], cc["g_ext"]
    z0 = cc["z0"]
    v0 = c.e0_v_m * B_WG
    i_d = g_ext * v0
    p_fwd = g_ext * v0 ** 2 / 8
    solver, u = initial_state(c, prm)
    area = solver.volumes                      # annulus areas (per unit rod length)
    v = complex((c.e_start_v_m or c.e0_v_m) * B_WG, 0.0)
    e_init = (u[:3] * area).sum(axis=1)
    w_cav0 = 0.5 * cap * abs(v) ** 2
    t = 0.0
    times = np.linspace(0, c.t_end_s, c.samples)
    rec = {k: [] for k in ("t_s", "e_gap_v_m", "e_axial_v_m", "r_ohm", "x_ohm", "p_rod_w", "p_net_in_w",
                           "p_refl_w", "w_cav_j", "r_cond_m", "ne_axis_m3", "tg_axis_k", "te_axis_ev",
                           "rho_axis_rel", "i_rod_a")}
    prof = {k: [] for k in ("ne", "tg", "te", "rho", "p", "vel", "sigma")}
    e_dep = e_sink = e_rad = 0.0
    w_in = w_rod = 0.0
    steps = 0
    drive_on = True

    def rod_impedance(uu):
        sig = conductivity(uu, w, prm, tabs, gm.ln_e0, gm.h)
        gl = float(sig @ area)                       # conductance per unit length x length -> S m
        r_ohm = B_WG / max(gl, 1e-30)
        # conducting radius: radius enclosing 90 % of the conductance
        cum = np.cumsum(sig * area)
        r_c = float(np.interp(0.9 * cum[-1], cum, solver.edges[1:])) if cum[-1] > 0 else c.r0_m
        x_ohm = post_reactance(max(r_c, 5e-6), c.f_hz) * z0
        return sig, r_ohm, x_ohm, r_c

    sig, r_ohm, x_ohm, r_c = rod_impedance(u)
    for target in times:
        while t < target * (1 - 1e-12):
            dt = min(target - t, c.dt_max_s, solver.stable_dt(u))
            if c.drive_off_s >= 0 and t >= c.drive_off_s:
                drive_on = False
            idr = i_d if drive_on else 0.0
            z_rod = complex(r_ohm, x_ohm)
            y_rod = 1 / z_rod
            v1 = (2 * cap / dt * v + idr) / (2 * cap / dt + g_ext + y_rod)
            vm = 0.5 * (v + v1)
            i_rod = vm * y_rod
            p_rod = 0.5 * abs(i_rod) ** 2 * r_ohm
            p_in = 0.5 * (idr * vm.conjugate()).real - 0.5 * g_ext * abs(vm) ** 2
            e_ax = abs(i_rod) * r_ohm / B_WG
            p_ext = 0.5 * sig * e_ax ** 2 / QE        # eV m^-3 s^-1 (time-averaged Joule)
            d1, s1, r1 = chemistry_step(u, 0.5 * dt, p_ext, prm, tabs, gm.ln_e0, gm.h, gm.ln_i0, gm.hi,
                                        gm.ppt_tab)
            u = solver.step(u, dt)
            d2, s2, r2 = chemistry_step(u, 0.5 * dt, p_ext, prm, tabs, gm.ln_e0, gm.h, gm.ln_i0, gm.hi,
                                        gm.ppt_tab)
            e_dep += float((d1 + d2) @ area); e_sink += float((s1 + s2) @ area); e_rad += float((r1 + r2) @ area)
            w_in += p_in * dt
            w_rod += p_rod * dt
            v = v1
            t += dt
            steps += 1
            sig, r_ohm, x_ohm, r_c = rod_impedance(u)
        rho, vel, pres, gas = solver.primitive(u)
        particles = rho / M_HE - u[4] - u[7]
        tg = gas / (1.5 * KB * particles)
        ne = u[3] + u[4] + 2 * u[8]
        te = np.where(ne > 1e10, u[5] / np.maximum(ne, 1e-30) / 1.5, 0.0)
        z_rod = complex(r_ohm, x_ohm)
        i_rod = v / z_rod
        p_rod = 0.5 * abs(i_rod) ** 2 * r_ohm
        idr = i_d if drive_on else 0.0
        p_in = 0.5 * (idr * v.conjugate()).real - 0.5 * g_ext * abs(v) ** 2
        p_f = p_fwd if drive_on else 0.0
        rec["t_s"].append(t); rec["e_gap_v_m"].append(abs(v) / B_WG)
        rec["e_axial_v_m"].append(abs(i_rod) * r_ohm / B_WG)
        rec["r_ohm"].append(r_ohm); rec["x_ohm"].append(x_ohm); rec["p_rod_w"].append(p_rod)
        rec["p_net_in_w"].append(p_in); rec["p_refl_w"].append(p_f - p_in); rec["w_cav_j"].append(0.5 * cap * abs(v) ** 2)
        rec["r_cond_m"].append(r_c); rec["ne_axis_m3"].append(float(ne[0])); rec["tg_axis_k"].append(float(tg[0]))
        rec["te_axis_ev"].append(float(te[0])); rec["rho_axis_rel"].append(float(rho[0] / (c.n_gas * M_HE)))
        rec["i_rod_a"].append(abs(i_rod))
        for k, val in (("ne", ne), ("tg", tg), ("te", te), ("rho", rho / (c.n_gas * M_HE)), ("p", pres),
                       ("vel", vel), ("sigma", sig)):
            prof[k].append(np.asarray(val, float).copy())
        if verbose and len(rec["t_s"]) % 20 == 1:
            print(f"  t={t*1e9:7.2f} ns E_gap={abs(v)/B_WG/1e5:7.2f} kV/cm R={r_ohm:9.3g} X={x_ohm:6.0f} "
                  f"P_rod={p_rod/1e3:9.3g} kW ne0={ne[0]:.2e} Te0={te[0]:.2f} Tg0={tg[0]:.0f} rc={r_c*1e6:.0f}um "
                  f"steps={steps}", flush=True)
    e_fin = (u[:3] * area).sum(axis=1)
    w_cav = 0.5 * cap * abs(v) ** 2
    plasma_gain = e_fin[2] + solver.boundary_out[2] - e_init[2]
    summary = dict(config=asdict(c), cavity=cc, p_forward_w=p_fwd, w_cav0_j=w_cav0,
                   x_post_initial_ohm=rec["x_ohm"][0], r_initial_ohm=rec["r_ohm"][0],
                   w_rod_circuit_j=w_rod, w_rod_plasma_j=e_dep * B_WG, w_in_j=w_in, w_cav_end_j=w_cav,
                   circuit_energy_residual_rel=(w_in - w_rod - (w_cav - w_cav0)) / max(abs(w_in) + w_cav0, 1e-30),
                   plasma_energy_residual_rel=(plasma_gain - (e_dep - e_sink)) / max(e_dep, 1e-30),
                   coupling_residual_rel=(e_dep * B_WG - w_rod) / max(w_rod, 1e-30),
                   radiated_j=e_rad * B_WG, sink_j=e_sink * B_WG, steps=steps,
                   e_gap_min_kv_cm=min(rec["e_gap_v_m"]) / 1e5, tg_axis_max_k=max(rec["tg_axis_k"]),
                   ne_axis_max_m3=max(rec["ne_axis_m3"]), r_min_ohm=min(rec["r_ohm"]),
                   r_end_ohm=rec["r_ohm"][-1], r_cond_end_m=rec["r_cond_m"][-1],
                   skin_depth_end_m=float(math.sqrt(2 / (w * MU0 * max(sig.max(), 1e-30)))),
                   limitations=["single bridged rod at the field maximum; post reactance from the thin-post "
                                "Marcuvitz formula at the 90 % conductance radius",
                                "slowly varying envelope at the drive frequency (cavity detuning by the post "
                                "enters through Im Y_rod)",
                                "no sheaths, electrode erosion, diffusion, heat conduction or 3D current paths",
                                "CR closure and Euler gas restricted to moderate ionization degree"])
    if save:
        out = out_root / c.name
        out.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out / "rod_drive.npz", r_m=solver.centres, **{k: np.array(v_) for k, v_ in rec.items()},
                            **{"prof_" + k: np.array(v_) for k, v_ in prof.items()})
        (out / "summary.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    if verbose:
        print(json.dumps({k: summary[k] for k in ("w_cav0_j", "w_rod_circuit_j", "w_rod_plasma_j", "w_in_j",
                                                  "circuit_energy_residual_rel", "plasma_energy_residual_rel",
                                                  "coupling_residual_rel", "e_gap_min_kv_cm", "tg_axis_max_k",
                                                  "r_end_ohm")}, default=float), flush=True)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(prog="rod_drive")
    ap.add_argument("--name", required=True)
    ap.add_argument("--e0-kv-cm", type=float, default=200.0)
    ap.add_argument("--drive-off-ns", type=float, default=-1.0)
    ap.add_argument("--t-end-ns", type=float, default=300.0)
    ap.add_argument("--n0", type=float, default=3e23)
    ap.add_argument("--r0-um", type=float, default=100.0)
    ap.add_argument("--nr", type=int, default=150)
    a = ap.parse_args(argv)
    c = DriveConfig(name=a.name, e0_v_m=a.e0_kv_cm * 1e5,
                    drive_off_s=a.drive_off_ns * 1e-9 if a.drive_off_ns >= 0 else -1.0,
                    t_end_s=a.t_end_ns * 1e-9, n0_m3=a.n0, r0_m=a.r0_um * 1e-6, nr=a.nr,
                    samples=int(a.t_end_ns / 1.0) + 1)
    run(c)


if __name__ == "__main__":
    main()
