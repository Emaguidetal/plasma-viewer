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
  Re sigma = e^2 n_e nu / (m_e (nu^2 + omega^2)), nu = N k_m(eps) + nu_ei. Uniform E_ax needs skin depth >>
  radius: the skin depth is recorded and the first time it falls below 2 r_c is reported, separately for the
  hand-off transient (t <= 20 ns) and for the sustained run (t > 10 ns: t_skin_invalid_s).
* Plasma: radial cylindrical Euler (gas_dynamics.Euler1D, passive plasma tracers in the energy) plus the
  collisional-radiative helium cell (plasma3d_cr.kinetic_cell_cr) with the Joule power density
  p = Re(sigma) |E_ax|^2 / 2, E_ax = I_rod R / b uniform in r (tangential continuity). Strang splitting.
* Electrode sheaths, wall erosion, radial diffusion and heat conduction (<10 um per 100 ns in cold gas, but
  about 30 um per 100 ns at 13 kK, comparable to the contracted arc) and 3D current redistribution are
  neglected; the CR closure is restricted to moderate ionization. Values after contraction are qualitative.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, replace
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
    # 0.1 ns sampling over the first 20 ns resolves the rod-loaded (low-Q) field collapse after bridging
    early = np.linspace(0, min(20e-9, c.t_end_s), int(min(20e-9, c.t_end_s) / 0.1e-9) + 1)
    times = np.unique(np.round(np.concatenate([early, np.linspace(0, c.t_end_s, c.samples)]), 15))
    rec = {k: [] for k in ("t_s", "e_gap_v_m", "e_axial_v_m", "r_ohm", "x_ohm", "p_rod_w", "p_net_in_w",
                           "p_refl_w", "w_cav_j", "r_cond_m", "ne_axis_m3", "tg_axis_k", "te_axis_ev",
                           "rho_axis_rel", "i_rod_a", "skin_depth_m", "q_loaded")}
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
            # Crank-Nicolson: discrete energy closes exactly with the midpoint power bookkeeping below
            ya = 0.5 * (g_ext + y_rod)
            v1 = ((2 * cap / dt - ya) * v + idr) / (2 * cap / dt + ya)
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
        rec["skin_depth_m"].append(math.sqrt(2 / (w * MU0 * max(sig.max(), 1e-30))))
        rec["q_loaded"].append(w * cap / (g_ext + (1 / z_rod).real))
        for k, val in (("ne", ne), ("tg", tg), ("te", te), ("rho", rho / (c.n_gas * M_HE)), ("p", pres),
                       ("vel", vel), ("sigma", sig)):
            prof[k].append(np.asarray(val, float).copy())
        if verbose and len(rec["t_s"]) % 20 == 1:
            print(f"  t={t*1e9:7.2f} ns E_gap={abs(v)/B_WG/1e5:7.2f} kV/cm R={r_ohm:9.3g} X={x_ohm:6.0f} "
                  f"P_rod={p_rod/1e3:9.3g} kW ne0={ne[0]:.2e} Te0={te[0]:.2f} Tg0={tg[0]:.0f} rc={r_c*1e6:.0f}um "
                  f"steps={steps}", flush=True)
    e_fin = (u[:3] * area).sum(axis=1)
    w_cav = 0.5 * cap * abs(v) ** 2
    ts = np.array(rec["t_s"]); rcs = np.array(rec["r_cond_m"]); eg = np.array(rec["e_gap_v_m"])
    # arc onset = thermal contraction: conducting radius below half its running maximum (after the start-up)
    late = ts > 10e-9
    contr = late & (rcs < 0.5 * np.maximum.accumulate(np.where(late, rcs, 0.0)))
    if c.drive_off_s >= 0:
        contr &= ts < c.drive_off_s
    skin_bad = np.array(rec["skin_depth_m"]) < 2 * rcs
    rx = np.array(rec["r_ohm"]) < np.array(rec["x_ohm"])
    first = lambda m: float(ts[np.argmax(m)]) if m.any() else None
    early_m = ts <= 20e-9
    # the hand-off at the full pre-bridge field gives a ~ns over-ionisation transient; report it separately
    # from the sustained (t > 10 ns) R < X and skin-depth onsets
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
                   t_contraction_s=first(contr), t_r_below_x_s=first(rx & late),
                   t_skin_invalid_s=first(skin_bad & late),
                   transient_r_below_x_s=first(rx & early_m), transient_skin_invalid_s=first(skin_bad & early_m),
                   q_loaded_min_first_20ns=float(np.array(rec["q_loaded"])[early_m].min()),
                   e_gap_min_first_20ns_kv_cm=float(eg[early_m].min() / 1e5),
                   t_e_gap_min_first_20ns_s=float(ts[early_m][np.argmin(eg[early_m])]),
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


def from_rod3d(name, t_off_after_laser_s=100e-9, root=None):
    """Hand-off from a bridged rod3d run: cavity field, channel density and radius at the bridge time.
    Returns DriveConfig kwargs; drive_off_s is converted to time after bridging."""
    root = root or Path(__file__).resolve().parents[2] / "results" / "rod3d"
    s = json.loads((root / name / "summary.json").read_text(encoding="utf-8"))
    if s["t_bridge_ns"] is None:
        raise ValueError(f"rod3d case {name} never bridges")
    with np.load(root / name / "rod3d.npz") as z:
        h = {k: z[k] for k in z.files}
    tb = s["t_bridge_ns"] * 1e-9
    k = int(np.searchsorted(h["t_s"], tb, side="right")) - 1   # last frame at or before the bridge (as the viewer)
    k = min(max(k, 0), len(h["t_s"]) - 1)
    m = int(h["m_frame"][k])
    kind, rod, lam, rad = h["kind"][:m], h["rod"][:m], h["lam"][k, :m], h["rad"][:m]
    ch = kind == 1
    rods = np.unique(rod[ch & (rod >= 0)])
    ext = {r: np.ptp(h["pos"][:m][rod == r, 1]) for r in rods}
    rb = max(ext, key=ext.get)                       # bridged rod = largest extent along E
    sel = ch & (rod == rb)
    r0 = float(np.average(rad[sel], weights=lam[sel]))
    # series resistance ~ sum ds / lambda: the harmonic mean line density preserves it (uniform mobility)
    lam_h = float(1.0 / np.mean(1.0 / np.maximum(lam[sel], 1e-30)))
    e_start = float(np.interp(tb, h["hist_t_s"], h["hist_e0_v_m"]))
    return dict(e0_v_m=s["config"]["e0_v_m"], e_start_v_m=e_start,
                n0_m3=lam_h / (math.pi * r0 ** 2), r0_m=r0,
                drive_off_s=max(t_off_after_laser_s - tb, 0.0) if t_off_after_laser_s >= 0 else -1.0,
                t_bridge_s=tb)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="rod_drive")
    ap.add_argument("--name", required=True)
    ap.add_argument("--e0-kv-cm", type=float, default=200.0)
    ap.add_argument("--drive-off-ns", type=float, default=-1.0)
    ap.add_argument("--t-end-ns", type=float, default=300.0)
    ap.add_argument("--n0", type=float, default=3e23)
    ap.add_argument("--r0-um", type=float, default=None, help="initial radius (default 100 um, or rod3d hand-off)")
    ap.add_argument("--nr", type=int, default=150)
    ap.add_argument("--dt-max-ns", type=float, default=0.05)
    ap.add_argument("--power-scale", type=float, default=1.0,
                    help="scale forward power and stored cavity energy (fields by its square root), e.g. 1/3.91^2 "
                         "for an illustrative reduced-power sensitivity; line density is kept")
    ap.add_argument("--from-rod3d", default=None,
                    help="start from the bridged rod of this rod3d case; --drive-off-ns is then after the laser")
    a = ap.parse_args(argv)
    c = DriveConfig(name=a.name, e0_v_m=a.e0_kv_cm * 1e5,
                    drive_off_s=a.drive_off_ns * 1e-9 if a.drive_off_ns >= 0 else -1.0,
                    t_end_s=a.t_end_ns * 1e-9, n0_m3=a.n0, r0_m=(a.r0_um or 100.0) * 1e-6, nr=a.nr,
                    dt_max_s=a.dt_max_ns * 1e-9, samples=int(a.t_end_ns / 1.0) + 1)
    if a.from_rod3d:
        hand = from_rod3d(a.from_rod3d, a.drive_off_ns * 1e-9 if a.drive_off_ns >= 0 else -1.0)
        print("hand-off from rod3d:", json.dumps(hand), flush=True)
        hand.pop("t_bridge_s")
        if a.r0_um:                                   # keep the handed-off line density at a new radius
            hand["n0_m3"] *= (hand["r0_m"] / (a.r0_um * 1e-6)) ** 2
            hand["r0_m"] = a.r0_um * 1e-6
        c = replace(c, **hand)
    if a.power_scale != 1.0:
        s = math.sqrt(a.power_scale)
        c = replace(c, e0_v_m=c.e0_v_m * s, e_start_v_m=(c.e_start_v_m or c.e0_v_m) * s)
    run(c)


if __name__ == "__main__":
    main()
