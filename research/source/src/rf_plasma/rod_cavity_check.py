"""Test of the lumped cavity coupling in rod_drive.py against the project's 3D FDTD (fdtd3d.py).

Question: a fixed straight rod of radius r0 bridges the full 34 mm gap along E at the TE101 field maximum of
the iris-coupled WR-284 cavity. Does rod_drive's lumped resonator (C, G_ext = omega C / Q_L, shunt post
Z_rod = R + j X_Marcuvitz) predict the loaded resonance shift, loaded Q and the fraction of forward power
absorbed in the rod?

What is solved
* FDTD (fdtd3d.FDTD, PEC walls, z-CPML behind the feed) with a TE10 soft source. The rod is a thin wire on
  the Ey edges (i0 = na/2, all j, k0 = cavity centre node):
  - subcell model ("modified permeability"): the four H components that circle the rod edge get their
    Ey-difference coefficient multiplied by f(r0/dx). "lattice" (default) uses the exact 5-point lattice
    Green's function, G(0)-G(r) = (ln r + gamma + 1.5 ln 2)/(2 pi), so that the near-wire bonds carry the
    continuum inductance (mu0/2pi) ln(R/r0) instead of the lattice value: f = (pi/2)/(ln(dx/r0) - C_L),
    C_L = gamma + 1.5 ln 2 - pi/2 = 0.0461 (f = 1 at r0 = 0.1985 dx, the equivalent radius of a plain PEC
    Yee edge). "tu" is the Umashankar-Taflove 1/rho contour model f = 2/ln(dx/r0); "none" is a plain edge.
  - rod material: PEC (edge masked) or the solver's own Crank-Nicolson Drude current on the rod edge with
    (n_e, nu) and weight pi r0^2/dx^2 (uniform current over the rod cross section, no skin effect, as in
    rod_drive).
* The full-height post problem is y-invariant (TE_m0 fields, inductive iris, full-gap rod), so the same 3D
  solver is run on a 2-cell-high box ("thin"); "full" runs the 34-cell height for an identity check.
* Excitation: modulated Gaussian pulse; outputs from one run:
  - complex poles of the cavity field (matrix pencil on a probe after the pulse) -> f, Q_L = pi f tau;
  - DFT spectra of forward/backward TE10 waves at the feed plane and of E, J on the rod edge ->
    P_rod/P_fwd(f) for steady drive at any f (linear system), and the closure P_fwd = P_back + P_rod;
  - an optional CW steady-drive run checks the DFT fraction directly.

Predictions compared
* rod_drive (its own functions, f_hz and q_l set to the FDTD empty cavity at the same dx): pole of the
  envelope equation 2C dV/dt = -(G_ext + Y_rod)V, i.e. df = -Im Y/(4 pi C), Q = omega C/(G_ext + Re Y)
  (= recorded q_loaded), absorbed fraction 4 G_ext Re Y/|G_ext + Y|^2 at the drive frequency; X evaluated
  at rod_drive's 90 % conductance radius (sqrt(0.9) r0 for a uniform rod).
* "tl" (proposed replacement, not applied to rod_drive): single-mode transmission line feed | inductive
  iris | line to the post | post R + jX(omega) | line to the short, fitted to the same empty f0 and Q_L;
  complex natural frequency and 1 - |Gamma|^2 at the drive frequency.

CLI: validate | empty | case | sweep | steady | report   (python -m rf_plasma.rod_cavity_check -h)
"""
from __future__ import annotations

import argparse
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import numba
import numpy as np
from numba import njit
from scipy.constants import epsilon_0 as EPS0, mu_0 as MU0, c as C0
from scipy.optimize import fsolve
from scipy.signal import hilbert

from . import rod_drive as rd
from .fdtd3d import FDTD, _plane_projection, _update_e, _update_h
from .geometry3d import TE10Port, iris_cavity, straight_guide

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "rod_cavity_check"
REF = json.loads((ROOT / "cases" / "reference.json").read_text(encoding="utf-8"))
QE, ME = 1.602176634e-19, 9.1093837015e-31
ETA0 = math.sqrt(MU0 / EPS0)
APERTURE_M = 0.020                       # calibrated iris (results/coupled3d/calibration.json)
LATTICE_C = 0.5772156649015329 + 1.5 * math.log(2.0) - 0.5 * math.pi
NE_CASES = (1e21, 1e22, 1e23, 1e24, 1e25)
R0_CASES = (0.10e-3, 0.15e-3, 0.20e-3, 0.30e-3)
R0DX_MIN, R0DX_MAX = 0.04, 0.65         # validated thin-wire range (validation_post*.json)
EPS_EV = 1.5                             # rod_drive initial mean electron energy (DriveConfig.eps0_ev)


# ---------------------------------------------------------------- thin-wire subcell model
def wire_factor(r0: float, dx: float, model: str = "lattice") -> float:
    """Multiplier of the Ey-difference term in the four H updates circling the wire edge."""
    if model == "none":
        return 1.0
    if model == "tu":
        return 2.0 / math.log(dx / r0)
    if model == "lattice":
        den = math.log(dx / r0) - LATTICE_C
        if den <= 0.05:
            raise ValueError(f"r0/dx = {r0 / dx:.3f} is outside the thin-wire range")
        return 0.5 * math.pi / den
    raise ValueError(model)


@njit(cache=True)
def _wire_h(Ey, Hx, Hz, i0, k0, coef):
    """Add (f-1) x the standard Ey-difference term to the H components circling Ey[i0, :, k0]."""
    for j in range(Ey.shape[1]):
        Hz[i0, j, k0] -= coef * (Ey[i0 + 1, j, k0] - Ey[i0, j, k0])
        Hz[i0 - 1, j, k0] -= coef * (Ey[i0, j, k0] - Ey[i0 - 1, j, k0])
        Hx[i0, j, k0] += coef * (Ey[i0, j, k0 + 1] - Ey[i0, j, k0])
        Hx[i0, j, k0 - 1] += coef * (Ey[i0, j, k0] - Ey[i0, j, k0 - 1])


@dataclass
class Rod:
    r0_m: float
    kind: str = "pec"            # "pec" | "drude" | "none" (no wire: empty reference)
    ne_m3: float = 0.0           # Drude rod density
    nu_s: float = 0.0            # momentum-transfer collision frequency
    model: str = "lattice"       # subcell coefficient model

    def sigma(self, f_hz: float) -> complex:
        w = 2 * math.pi * f_hz
        return QE * QE * self.ne_m3 / (ME * (self.nu_s + 1j * w)) if self.kind == "drude" else complex("inf")

    def resistance(self, f_hz: float, length_m: float = rd.B_WG) -> float:
        """rod_drive convention: R = b / (Re sigma * area) (0 for PEC, inf for no rod)."""
        if self.kind == "pec":
            return 0.0
        if self.kind == "none":
            return math.inf
        return length_m / (self.sigma(f_hz).real * math.pi * self.r0_m ** 2)


class WireFDTD(FDTD):
    """FDTD.step with the thin-wire H correction between the H and E updates."""

    wire_coef = 0.0
    wire_ik = None

    def install_rod(self, i0: int, k0: int, rod: Rod):
        g = self.grid
        self.wire_ik = (int(i0), int(k0))
        self.rod_f = 1.0 if rod.kind == "none" else wire_factor(rod.r0_m, g.dx, rod.model)
        self.wire_coef = self.ch * (self.rod_f - 1.0)
        if rod.kind == "pec":
            self.mEy[i0, :, k0] = 0
            self.Ey[i0, :, k0] = 0.0
        elif rod.kind == "drude":
            wp2 = QE * QE * rod.ne_m3 / (EPS0 * ME) * math.pi * rod.r0_m ** 2 / g.dx ** 2
            d = 1.0 + 0.5 * rod.nu_s * g.dt
            self.aEy[i0, :, k0] = (1.0 - 0.5 * rod.nu_s * g.dt) / d
            self.gEy[i0, :, k0] = g.dt * g.dt * wp2 / (4.0 * d)

    def step(self, source=None):
        _update_h(self.Ex, self.Ey, self.Ez, self.Hx, self.Hy, self.Hz, self.ch,
                  self.psi_hx, self.psi_hy, self.bh, self.chz)
        if self.wire_coef != 0.0:
            _wire_h(self.Ey, self.Hx, self.Hz, self.wire_ik[0], self.wire_ik[1], self.wire_coef)
        _update_e(self.Ex, self.Ey, self.Ez, self.Hx, self.Hy, self.Hz, self.ce, self.grid.dt,
                  self.mEx, self.mEy, self.mEz, self.aEx, self.aEy, self.aEz,
                  self.gEx, self.gEy, self.gEz, self.Jx, self.Jy, self.Jz,
                  self.Px, self.Py, self.Pz, self.psi_ex, self.psi_ey, self.be, self.cez)
        if source is not None:
            source(self)
        self.step_count += 1


def courant_for(f: float) -> float:
    """Default 0.98 of the 3D limit; reduced when the subcell factor stiffens the wire bonds (f > 1.4)."""
    return 0.98 * min(1.0, math.sqrt(1.4 / f)) if f > 1.0 else 0.98


# ---------------------------------------------------------------- geometry
def _width(dx):
    return 2 * round(REF["rectangular"]["width_m"] / (2 * dx)) * dx   # even cell count: rod edge on x = a/2


def cavity_layout(dx: float, height: str = "thin", courant: float = 0.98):
    r = REF["rectangular"]
    h = 2 * dx if height == "thin" else r["height_m"]
    lay = iris_cavity(dx, _width(dx), h, r["closed_cavity_length_m"], APERTURE_M)
    _set_courant(lay, courant)
    k_iris1 = lay.gas_k_lo - 1
    return lay, k_iris1


def guide_layout(dx: float, courant: float = 0.98, length_m: float = 0.15):
    lay = straight_guide(dx, _width(dx), 2 * dx, length_m)
    _set_courant(lay, courant)
    return lay


def _set_courant(lay, courant):
    g = lay.grid
    g.courant = courant
    g.dt = courant * g.dx / (C0 * math.sqrt(3.0))


# ---------------------------------------------------------------- signal processing
def dft(sig: np.ndarray, dt: float, freqs: np.ndarray, shift: float = 0.0, chunk: int = 20000) -> np.ndarray:
    """sum_n sig[n] exp(-i w (n + shift) dt) at each frequency (exact, chunked)."""
    out = np.zeros(len(freqs), complex)
    w = 2 * math.pi * np.asarray(freqs)
    for s in range(0, len(sig), chunk):
        n = np.arange(s, min(s + chunk, len(sig))) + shift
        out += np.exp(-1j * np.outer(w, n * dt)) @ sig[s:s + len(n)]
    return out


def matrix_pencil(y: np.ndarray, dt: float, max_order: int = 24, rel_tol: float = 1e-7):
    """Damped complex exponentials y[n] = sum a_m z_m^n (Hua & Sarkar 1990). Returns f, alpha (1/s), amp."""
    n = len(y)
    L = n // 3
    Y = np.lib.stride_tricks.sliding_window_view(y, L + 1)
    _, s, vh = np.linalg.svd(Y, full_matrices=False)
    m = int(min(max_order, np.sum(s > rel_tol * s[0])))
    v = vh[:m].conj().T
    z = np.linalg.eigvals(np.linalg.pinv(v[:-1]) @ v[1:])
    zz = z[None, :] ** np.arange(n)[:, None]
    a = np.linalg.lstsq(zz, y.astype(complex), rcond=None)[0]
    return np.angle(z) / (2 * math.pi * dt), -np.log(np.abs(z)) / dt, a


def dominant_mode(sig, dt, t0_idx, f_lo, f_hi, f_s_target=36e9, max_samples=1500):
    """Pencil fit on a decimated window after index t0_idx; dominant pole with f in [f_lo, f_hi]."""
    m = max(1, int(round(1 / (f_s_target * dt))))
    y = sig[t0_idx::m]
    env = np.abs(hilbert(y))
    keep = np.nonzero(env > 1e-9 * env.max())[0]
    y = y[: (keep[-1] + 1 if len(keep) else len(y))][:max_samples]
    f, alpha, a = matrix_pencil(y, m * dt)
    sel = (f > f_lo) & (f < f_hi) & (alpha > 0)
    if not sel.any():
        return dict(f_hz=float("nan"), tau_s=float("nan"), q=float("nan"), window_samples=len(y))
    # weight by the energy of each term over the window
    energy = np.abs(a) ** 2 / np.maximum(1 - np.exp(-2 * alpha * m * dt), 1e-12)
    i = np.nonzero(sel)[0][np.argmax(energy[sel])]
    tau = 1 / alpha[i]
    return dict(f_hz=float(f[i]), tau_s=float(tau), q=float(math.pi * f[i] * tau), window_samples=len(y),
                window_s=float(len(y) * m * dt), decimation=m,
                others=[(float(f[k]), float(math.pi * f[k] / alpha[k])) for k in np.nonzero(sel)[0] if k != i])


def gaussian_pulse(port, f_c, sigma_t, t_c):
    port.f = f_c
    port.amplitude = 1.0
    port.envelope = lambda t: math.exp(-0.5 * ((t - t_c) / sigma_t) ** 2)


# ---------------------------------------------------------------- one cavity run
@dataclass
class CaseSpec:
    dx_m: float = 1e-3
    r0_m: float = 0.2e-3
    kind: str = "pec"          # pec | drude | none
    ne_m3: float = 0.0
    nu_s: float = 0.0
    label: str = ""
    model: str = "lattice"
    height: str = "thin"
    t_max_s: float = 300e-9
    stop_rel: float = 3e-3     # stop when the rod-plane field envelope < stop_rel x its maximum
    f_c_hz: float = 3.3e9
    sigma_t_s: float = 0.3e-9
    mode_shape: bool = False

    @property
    def name(self):
        h = "" if self.height == "thin" else "_full"
        m = "" if self.model == "lattice" else f"_{self.model}"
        return f"dx{self.dx_m * 1e3:g}_r{self.r0_m * 1e6:.0f}_{self.label or self.kind}{m}{h}"


def run_cavity(spec: CaseSpec, threads: int = 4, f_eval=None, verbose=True) -> dict:
    numba.set_num_threads(threads)
    t_wall = time.time()
    rod = Rod(spec.r0_m, spec.kind, spec.ne_m3, spec.nu_s, spec.model)
    fac = 1.0 if spec.kind == "none" else wire_factor(spec.r0_m, spec.dx_m, spec.model)
    lay, k_iris1 = cavity_layout(spec.dx_m, spec.height, courant_for(fac))
    g = lay.grid
    fd = WireFDTD(g, lay.node_inside, lay.pml_lo, lay.pml_hi, 3e9)
    i0, k0 = g.nx // 2, lay.focus_k
    fd.install_rod(i0, k0, rod)
    port = TE10Port(lay, spec.f_c_hz)
    t_c = 5 * spec.sigma_t_s
    gaussian_pulse(port, spec.f_c_hz, spec.sigma_t_s, t_c)
    kq = k_iris1 + (k0 - k_iris1) // 2                      # probe a quarter cavity from the iris
    n_max = int(spec.t_max_s / g.dt)
    rec = {k: np.zeros(n_max) for k in ("probe", "probe_x", "e_rod", "j_rod", "v", "i", "e_nb")}
    per = int(round(1 / (3e9 * g.dt)))
    chk = 50 * per
    i_x = i0 + g.nx // 4
    j0 = 0
    plane, plane_n, plane_every, f_plane = None, 0, 8, None
    if spec.mode_shape:
        f_plane = f_eval[0]
        plane = np.zeros((g.nx + 1, g.nz + 1), complex)
    t_pulse_end = t_c + 5 * spec.sigma_t_s
    n = 0
    env_max = 0.0
    while n < n_max:
        fd.step(port.inject)
        rec["probe"][n] = fd.Ey[i0, j0, kq]
        rec["probe_x"][n] = fd.Ey[i_x, j0, k0]
        rec["e_rod"][n] = fd.Ey[i0, j0, k0]
        rec["e_nb"][n] = fd.Ey[i0 + 1, j0, k0]
        rec["j_rod"][n] = fd.Jy[i0, j0, k0]
        rec["v"][n], rec["i"][n] = _plane_projection(fd.Ey, fd.Hx, lay.meas_k, port.profile, port.profile)
        if plane is not None and n % plane_every == 0 and (n + 1) * g.dt > t_pulse_end:
            plane += fd.Ey[:, j0, :] * np.exp(-2j * math.pi * f_plane * (n + 1) * g.dt)
            plane_n += 1
        n += 1
        if n % chk == 0 and n * g.dt > t_pulse_end + 2e-9:
            # envelope of the field near the rod and at the quarter probe
            a = max(np.abs(rec["probe"][n - per:n]).max(), np.abs(rec["e_nb"][n - per:n]).max())
            env_max = max(env_max, np.abs(rec["probe"][:n]).max(), np.abs(rec["e_nb"][:n]).max())
            if a < spec.stop_rel * env_max:
                break
    for k in rec:
        rec[k] = rec[k][:n]
    dt = g.dt
    res = dict(spec=asdict(spec), name=spec.name, nx=g.nx, ny=g.ny, nz=g.nz, dt_s=dt, courant=g.courant,
               wire_factor=fac, steps=n, t_end_s=n * dt, i0=i0, k0=k0, k_iris1=int(k_iris1),
               rod_z_offset_m=float((k0 - k_iris1) * g.dx - (g.nz - k_iris1) * g.dx / 2),
               a_grid_m=g.nx * g.dx, cavity_len_grid_m=(g.nz - k_iris1) * g.dx)
    env_all = max(np.abs(rec["probe"]).max(), 1e-300)
    res["final_envelope_rel"] = float(np.abs(rec["probe"][-per:]).max() / env_all)
    # poles: probe and (if field exists there) the rod-adjacent edge
    t0 = int((t_pulse_end + 1e-9) / dt)
    res["mode"] = dominant_mode(rec["probe"], dt, t0, 2.3e9, 4.6e9)
    res["mode_x"] = dominant_mode(rec["probe_x"], dt, t0, 2.3e9, 4.6e9)
    # spectra
    freqs = np.unique(np.concatenate([np.asarray(f_eval or [], float), np.linspace(2.4e9, 4.4e9, 201)]))
    res["spectrum"] = spectra(rec, lay, port, freqs, dt, g, spec.kind)
    if plane is not None:
        res["mode_shape"] = mode_shape(plane, lay, k_iris1, i0, k0)
    res["wall_s"] = time.time() - t_wall
    if verbose:
        m = res["mode"]
        print(f"{spec.name}: f {m['f_hz'] / 1e9:.5f} GHz Q {m['q']:.1f} steps {n} "
              f"env_end {res['final_envelope_rel']:.1e} wall {res['wall_s']:.0f} s", flush=True)
    return res


def spectra(rec, lay, port0, freqs, dt, g, kind):
    """Forward/backward TE10 power and rod absorption per frequency (same DFT scaling)."""
    V = dft(rec["v"], dt, freqs)
    I = dft(rec["i"], dt, freqs, shift=-0.5)
    out = dict(f_hz=freqs.tolist(), p_fwd=[], p_back=[], p_rod=[])
    if kind == "drude":
        E = dft(rec["e_rod"], dt, freqs)
        J = dft(rec["j_rod"], dt, freqs)
        # CN: J and E are co-located in time; per edge volume dx^3, ny identical edges
        prod = 0.5 * (E * J.conj()).real * g.dx ** 3 * g.ny
    else:
        prod = np.zeros(len(freqs))
    for q, f in enumerate(freqs):
        p = TE10Port(lay, f)
        ve = V[q] / p.norm
        vh = I[q] / p.norm / math.cos(p.beta * g.dx / 2)
        fw, bw = 0.5 * (ve - p.z_te * vh), 0.5 * (ve + p.z_te * vh)
        out["p_fwd"].append(float(p.power(np.array([fw]))[0]))
        out["p_back"].append(float(p.power(np.array([bw]))[0]))
        out["p_rod"].append(float(prod[q]))
    pf, pb, pr = (np.array(out[k]) for k in ("p_fwd", "p_back", "p_rod"))
    out["rod_fraction"] = (pr / pf).tolist()
    out["closure_residual"] = ((pf - pb - pr) / pf).tolist()
    return out


def mode_shape(plane, lay, k_iris1, i0, k0):
    """Phase-aligned real mode shape; deviations of E_y(x) and E_y(z) from the TE101 sine forms."""
    g = lay.grid
    ph = plane[i0, k0] / abs(plane[i0, k0])
    e = (plane / ph)
    cav = e[:, k_iris1:]
    imag_rel = float(np.abs(cav.imag).max() / np.abs(cav.real).max())
    en = np.abs(plane / plane[i0, k0]) ** 2
    a_g, d_g = g.nx * g.dx, (g.nz - k_iris1) * g.dx
    # effective capacitance relative to rod_drive's eps0 a d / 4b: C_eff/C = int |e|^2 dA / (a d / 4), e = 1 at
    # the rod (gap voltage V = E_rod b); cavity only, and cavity plus iris and evanescent feed region (16 mm)
    k_feed = max(lay.pml_lo + 1, k_iris1 - int(round(0.018 / g.dx)))
    c_ratio_cav = float(en[:, k_iris1:].sum() * g.dx ** 2 / (a_g * d_g / 4))
    c_ratio_ext = float(en[:, k_feed:].sum() * g.dx ** 2 / (a_g * d_g / 4))
    e = e.real / e.real[i0, k0]
    x = np.arange(g.nx + 1) * g.dx
    ex = e[:, k0]
    sx = np.sin(math.pi * x / (g.nx * g.dx))
    kk = np.arange(k_iris1, g.nz + 1)
    z = (kk - k_iris1) * g.dx
    ez = e[i0, kk]
    d = (g.nz - k_iris1) * g.dx
    sz = np.sin(math.pi * z / d)
    # best-fit sine in z: E = A sin(pi (z - z1)/L), free z1 and L (iris end is not a perfect short)
    from scipy.optimize import curve_fit
    popt, _ = curve_fit(lambda zz, A, z1, L: A * np.sin(math.pi * (zz - z1) / L), z, ez, p0=(1, 0, d))
    fit = popt[0] * np.sin(math.pi * (z - popt[1]) / popt[2])
    kmax = int(np.argmax(ez))
    return dict(x_m=x.tolist(), ey_x=ex.tolist(), z_m=z.tolist(), ey_z=ez.tolist(),
                dev_x_max=float(np.abs(ex - sx).max()), dev_x_rms=float(np.sqrt(np.mean((ex - sx) ** 2))),
                dev_z_max=float(np.abs(ez - sz).max()), dev_z_rms=float(np.sqrt(np.mean((ez - sz) ** 2))),
                dev_z_fit_max=float(np.abs(ez - fit).max()), z_fit_offset_m=float(popt[1]),
                z_fit_length_m=float(popt[2]), z_peak_m=float(z[kmax]), z_rod_m=float((k0 - k_iris1) * g.dx),
                ey_at_rod_rel_peak=float(ez[k0 - k_iris1] / ez.max()), quadrature_rel_cavity=imag_rel,
                c_eff_over_c_cavity=c_ratio_cav, c_eff_over_c_with_iris_region=c_ratio_ext)


# ---------------------------------------------------------------- empty-cavity checks
def hilbert_ringdown(dx: float, height: str = "thin", t_ring_s: float = 30e-9, threads: int = 8) -> dict:
    """Exactly the coupled3d.ringdown procedure (F_NOMINAL pulse, Hilbert slope) on this layout."""
    numba.set_num_threads(threads)
    lay, _ = cavity_layout(dx, height)
    g = lay.grid
    fd = FDTD(g, lay.node_inside, lay.pml_lo, lay.pml_hi, 3e9)
    port = TE10Port(lay, 3e9)
    port.amplitude = 1.0
    t_c, sig = 2.5e-9, 0.6e-9
    port.envelope = lambda t: math.exp(-0.5 * ((t - t_c) / sig) ** 2)
    i, j, k = g.nx // 2, g.ny // 2 - (1 if g.ny > 2 else 0), lay.focus_k
    n = int((2 * t_c + t_ring_s) / g.dt)
    ey = np.empty(n)
    for s in range(n):
        fd.step(port.inject)
        ey[s] = fd.Ey[i, j, k]
    t = np.arange(n) * g.dt
    sel = t > 2 * t_c + 3e-9
    an = hilbert(ey[sel])
    env, ph = np.log(np.abs(an)), np.unwrap(np.angle(an))
    tt = t[sel]
    trim = slice(len(tt) // 10, -len(tt) // 10)
    slope = np.polyfit(tt[trim], env[trim], 1)[0]
    f0 = np.polyfit(tt[trim], ph[trim], 1)[0] / (2 * math.pi)
    return dict(dx_m=dx, height=height, ny=g.ny, f0_hz=float(f0), tau_amp_s=float(-1 / slope),
                q_loaded=float(-math.pi * f0 / slope))


def prefill_profiles():
    """Full-height 3D steady-state snapshot from coupled3d prefill: y uniformity and x/z sine deviation."""
    p = ROOT / "results" / "coupled3d" / "prefill_cavity_dx1mm.npz"
    if not p.exists():
        return None
    ey = np.load(p)["Ey"]
    lay, k_iris1 = cavity_layout(1e-3, "full")
    g = lay.grid
    i0, k0 = g.nx // 2, lay.focus_k
    # standing wave snapshot: normalise by the cavity maximum of |Ey| on the mid-height plane
    jm = g.ny // 2
    pl = ey[:, jm, :]
    s = pl[i0, k0]
    x = np.arange(g.nx + 1) * g.dx
    ex = pl[:, k0] / s
    z = (np.arange(k_iris1, g.nz + 1) - k_iris1) * g.dx
    ez = pl[i0, k_iris1:] / s
    ey_y = ey[i0, :, k0] / s
    d = (g.nz - k_iris1) * g.dx
    return dict(snapshot_ey_focus=float(s), dev_x_max=float(np.abs(ex - np.sin(math.pi * x / (g.nx * g.dx))).max()),
                dev_z_max=float(np.abs(ez - np.sin(math.pi * z / d)).max()),
                y_nonuniformity=float(np.abs(ey_y - 1).max()))


# ---------------------------------------------------------------- waveguide post validation
def guide_post(dx: float, rod: Rod | None, threads: int = 4, freqs=(2.6e9, 2.95e9, 3.3e9, 3.6e9, 3.9e9),
               t_end_s: float = 24e-9, courant: float | None = None) -> dict:
    numba.set_num_threads(threads)
    fac = 1.0 if rod is None or rod.kind == "none" else wire_factor(rod.r0_m, dx, rod.model)
    lay = guide_layout(dx, courant if courant is not None else courant_for(fac))
    g = lay.grid
    fd = WireFDTD(g, lay.node_inside, lay.pml_lo, lay.pml_hi, 3e9)
    k_post = (lay.meas_k + lay.trans_k) // 2
    if rod is not None:
        fd.install_rod(g.nx // 2, k_post, rod)
    port = TE10Port(lay, 3.3e9)
    gaussian_pulse(port, 3.3e9, 0.3e-9, 1.5e-9)
    n = int(t_end_s / g.dt)
    v = np.zeros(n); cur = np.zeros(n)
    for s in range(n):
        fd.step(port.inject)
        v[s], cur[s] = _plane_projection(fd.Ey, fd.Hx, lay.trans_k, port.profile, port.profile)
    f = np.asarray(freqs)
    V, I = dft(v, g.dt, f), dft(cur, g.dt, f, shift=-0.5)
    fw = []
    for q, fq in enumerate(f):
        p = TE10Port(lay, fq)
        fw.append(0.5 * (V[q] / p.norm - p.z_te * I[q] / p.norm / math.cos(p.beta * g.dx / 2)))
    return dict(dx_m=dx, a_m=g.nx * g.dx, f_hz=f.tolist(), fwd=np.array(fw), wire_factor=fac, courant=g.courant,
                last_rel=float(np.abs(v[-200:]).max() / np.abs(v).max()))


def marcuvitz_x(r0: float, f, a: float) -> np.ndarray:
    """Normalised shunt reactance of a centred thin inductive post (rod_drive.post_reactance, any a, f)."""
    f = np.asarray(f, complex)
    lam = C0 / f
    lam_g = lam / np.sqrt(1 - (lam / (2 * a)) ** 2)
    k = 2 * a / lam
    s = sum(1 / np.sqrt(n * n - k * k) - 1 / n for n in range(3, 400, 2))
    return (a / (2 * lam_g)) * (np.log(4 * a / (math.pi * 2 * r0)) - 2 + 2 * s)


def validate(dxs=(1e-3, 0.5e-3, 0.25e-3), r0s=(0.05e-3, 0.1e-3, 0.15e-3, 0.2e-3, 0.3e-3),
             models=("lattice", "tu"), threads=8, out=OUT, tag="", r0dx_max=R0DX_MAX):
    """Post in a matched guide: z_n = 2/S21 - ... -> compare with Marcuvitz and the Drude series resistance."""
    rows = []
    for dx in dxs:
        refs = {}
        cases = [("none", Rod(0.1985 * dx, "pec", model="none"))]
        for r0 in r0s:
            if not R0DX_MIN <= r0 / dx <= r0dx_max:
                continue
            for m in models:
                cases.append((m, Rod(r0, "pec", model=m)))
            # resistive checks (rod_drive-type plasma): n_e 1e23 and 1e24 at 1.5 eV
            if any(math.isclose(r0, q) for q in (0.1e-3, 0.2e-3, 0.3e-3)):
                for ne in (1e23, 1e24):
                    s_re, nu = rod_sigma(ne)
                    cases.append(("lattice", Rod(r0, "drude", ne, nu, "lattice")))
        for m, rod in cases:
            try:
                r = guide_post(dx, rod, threads)
            except ValueError as e:
                print("skip", e)
                continue
            if r["courant"] not in refs:                # reference at the same dt (same grid dispersion)
                refs[r["courant"]] = guide_post(dx, None, threads, courant=r["courant"])
            s21 = r["fwd"] / refs[r["courant"]]["fwd"]
            y = 2 / s21 - 2                      # shunt element normalised admittance
            z = 1 / y
            a = r["a_m"]
            xm = marcuvitz_x(rod.r0_m, r["f_hz"], a).real
            row = dict(dx_m=dx, r0_m=rod.r0_m, r0_over_dx=rod.r0_m / dx, model=m, kind=rod.kind,
                       ne_m3=rod.ne_m3, f_hz=r["f_hz"], wire_factor=r["wire_factor"], courant=r["courant"],
                       x_fdtd=z.imag.tolist(), r_fdtd=z.real.tolist(), x_marcuvitz=xm.tolist(),
                       x_err_rel=((z.imag - xm) / xm).tolist())
            if rod.kind == "drude":
                fz = np.asarray(r["f_hz"])
                lay = guide_layout(dx)
                zte = np.array([TE10Port(lay, fq).z_te for fq in fz])
                sig = np.array([rod.sigma(fq) for fq in fz])
                zp = 1 / (sig * math.pi * rod.r0_m ** 2)               # ohm per metre (complex Drude)
                zn = zp * a / (2 * zte)                                  # R/Z0, Z0 = 2 b Z_TE / a
                row["r_expected"] = zn.real.tolist()
                row["x_expected"] = (xm + zn.imag).tolist()
                row["r_err_rel"] = ((z.real - zn.real) / zn.real).tolist()
                row["x_err_rel"] = ((z.imag - row["x_expected"]) / np.array(row["x_expected"])).tolist()
            rows.append(row)
            k = r["f_hz"].index(2.95e9)
            msg = f"dx {dx*1e3:g} r0 {rod.r0_m*1e3:.4f} ({rod.r0_m/dx:.3f} dx) {m:7s} {rod.kind:5s} " \
                  f"x {z.imag[k]:.4f} vs {row.get('x_expected', xm)[k]:.4f} ({row['x_err_rel'][k]*100:+.2f} %)"
            if rod.kind == "drude":
                msg += f" r {z.real[k]:.4f} vs {row['r_expected'][k]:.4f} ({row['r_err_rel'][k]*100:+.2f} %)"
            print(msg, flush=True)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"validation_post{tag}.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
    return rows


# ---------------------------------------------------------------- predictions
_SIG_CACHE = {}


def rod_sigma(ne: float, eps_ev: float = EPS_EV):
    """rod_drive.conductivity for a uniform He2+ plasma at its configured gas state -> (Re sigma, nu)."""
    key = (ne, eps_ev)
    if key not in _SIG_CACHE:
        from .gas_dynamics import M_HE
        from .gas_model import load_gas_model
        gm = load_gas_model("He", 532e-9, allow_unverified=True, extended=True, cr=True)
        c = rd.DriveConfig()
        w = 2 * math.pi * c.f_hz
        u = np.zeros((9, 1)); u[0] = c.n_gas * M_HE; u[4] = ne; u[5] = ne * eps_ev
        s = float(rd.conductivity(u, w, gm.prm, gm.tabs, gm.ln_e0, gm.h)[0])
        A = QE * QE * ne / (ME * s)
        nu = 0.5 * (A + math.sqrt(A * A - 4 * w * w))       # Re sigma = e^2 n nu / (m (nu^2 + w^2))
        _SIG_CACHE[key] = (s, nu)
    return _SIG_CACHE[key]


def lumped_prediction(r0: float, R: float, f0: float, q_l: float, f_drive: float | None = None) -> dict:
    """rod_drive's envelope resonator: pole and steady absorbed fraction (its own constants/functions)."""
    c = rd.DriveConfig(f_hz=f0, q_l=q_l)
    cc = rd.cavity_constants(c)
    w, cap, g = cc["omega"], cc["cap"], cc["g_ext"]
    r_c = math.sqrt(0.9) * r0                         # rod_drive: radius enclosing 90 % of the conductance
    X = rd.post_reactance(r_c, f0) * cc["z0"]
    if math.isinf(R):
        Y = 0j
    else:
        Y = 1 / complex(R, X)
    df = -Y.imag / (2 * cap) / (2 * math.pi)
    q = w * cap / (g + Y.real)
    frac = 4 * g * Y.real / abs(g + Y) ** 2
    return dict(f_hz=f0 + df, df_hz=df, q=q, rod_fraction=frac, x_ohm=X, r_ohm=R, r_c_m=r_c,
                linearisation_df_over_f=df / f0)


class TLModel:
    """Single-mode TE10 transmission line: matched feed | iris jb(w) | cavity line with shunt post | short.

    The iris is an inductive susceptance b(w) = b0 (Z_TE(w)/Z_TE(w0)) (w0/w) (normalised to the guide), the
    cavity length d is effective; (b0, d) are fitted so that the empty natural frequency equals the FDTD f0,
    Q_L. The post sits at z_p (from the iris plane) with z_n = R/Z0(w) + j x_Marcuvitz(w) (Z0 = 2 b Z_TE/a).
    """

    def __init__(self, f0, q_l, a=rd.A_WG, b=rd.B_WG, d_guess=rd.D_CAV, post_offset_m=0.0):
        self.a, self.b, self.off = a, b, post_offset_m
        self.w0 = 2 * math.pi * f0
        self.q_l = q_l
        self.wc = math.pi * C0 / a
        target = self.w0 * (1 + 0.5j / q_l)
        best = None
        for b0 in (-5.0, -10.0, -20.0, -40.0):
            sol, info, ier, _ = fsolve(lambda p: self._res(p, target), (b0, d_guess), full_output=True)
            if ier == 1 and (best is None or sol[1] > 0):
                best = sol
                if 0.5 * d_guess < sol[1] < 1.5 * d_guess:
                    break
        self.b0, self.d = best
        self.zte0 = self._zte(self.w0)

    def _res(self, p, w):
        self.b0, self.d = p
        self.zte0 = self._zte(self.w0)
        v = self.char(w, None)
        return [v.real, v.imag]

    def _beta(self, w):
        return np.sqrt(w * w - self.wc ** 2 + 0j) / C0

    def _zte(self, w):
        return ETA0 * w / (C0 * self._beta(w))

    def z0(self, w):
        return 2 * self.b * self._zte(w) / self.a

    def y_load(self, w, R, r0):
        """Normalised admittance seen from the feed at the iris plane."""
        beta = self._beta(w)
        biris = self.b0 * (self._zte(w) / self.zte0) * (self.w0 / w)
        zp = 0.5 * self.d + self.off
        l2 = self.d - zp
        y_s = -1j / np.tan(beta * l2)
        if R is None or math.isinf(R):
            y_l = y_s
        else:
            f = w / (2 * math.pi)
            zn = R / self.z0(w) + 1j * marcuvitz_x(r0, f, self.a)
            y_l = y_s + 1 / zn
        t = np.tan(beta * zp)
        y_in = (y_l + 1j * t) / (1 + 1j * y_l * t)
        return 1j * biris + y_in

    def char(self, w, R, r0=1e-4):
        """Pole-free characteristic function (zero <=> source-free solution with a matched feed).

        Normalised (V, I) are carried from the short (V = j sin(beta l2), I = cos(beta l2)) through the
        post (scaled by z_post so that no division occurs) and the line to the iris; the natural
        frequencies satisfy I_in + j b V_in + V_in = 0 (input admittance + 1 = 0)."""
        beta = self._beta(w)
        biris = self.b0 * (self._zte(w) / self.zte0) * (self.w0 / w)
        zp = 0.5 * self.d + self.off
        l2 = self.d - zp
        V, I = 1j * np.sin(beta * l2), np.cos(beta * l2)
        if R is not None and not math.isinf(R):
            f = w / (2 * math.pi)
            zn = R / self.z0(w) + 1j * marcuvitz_x(r0, f, self.a)
            V, I = zn * V, zn * I + V
        c, s = np.cos(beta * zp), np.sin(beta * zp)
        V, I = V * c + 1j * I * s, I * c + 1j * V * s
        return I + 1j * biris * V + V

    def pole(self, r0, R, n_path=60):
        """Complex natural frequency, tracked by continuation in R from the empty cavity (R = 1e9 ohm)."""
        from scipy.optimize import newton
        w = self.w0 * (1 + 0.5j / self.q_l)
        r_end = max(R, 1e-3) if not math.isinf(R) else 1e9
        path = list(np.logspace(9, math.log10(r_end), n_path)) + ([0.0] if R == 0 else [])
        for Rk in path:
            w = newton(lambda ww: complex(self.char(ww, Rk, r0)), w, tol=1e-2, maxiter=100)
        return w

    def predict(self, r0, R, f_drive, guess_hz=None):
        w = self.pole(r0, R)
        f = w.real / (2 * math.pi)
        yl = complex(self.y_load(2 * math.pi * f_drive, R, r0))
        gam = (1 - yl) / (1 + yl)
        return dict(f_hz=f, q=w.real / (2 * w.imag), rod_fraction=1 - abs(gam) ** 2)


# ---------------------------------------------------------------- sweep / report
def case_specs(dx: float, r0s=R0_CASES, ne_list=NE_CASES, extra_r_ohm=(10.0,), height="thin"):
    specs = []
    for r0 in r0s:
        specs.append(CaseSpec(dx, r0, "pec", label="pec", height=height))
        for ne in ne_list:
            s, nu = rod_sigma(ne)
            specs.append(CaseSpec(dx, r0, "drude", ne, nu, label=f"ne{ne:.0e}".replace("+", ""), height=height))
        for R in extra_r_ohm:
            nu = 2e13
            w = 2 * math.pi * 2.95e9
            s_need = rd.B_WG / (R * math.pi * r0 ** 2)
            ne = s_need * ME * (nu * nu + w * w) / (QE * QE * nu)
            specs.append(CaseSpec(dx, r0, "drude", ne, nu, label=f"R{R:g}", height=height))
    return specs


def _worker(spec_dict, threads, f_eval, out):
    spec = CaseSpec(**spec_dict)
    res = run_cavity(spec, threads, f_eval, verbose=True)
    p = out / "runs" / f"{spec.name}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(res, default=float), encoding="utf-8")
    return spec.name, res["wall_s"]


def empty_reference(dx, height="thin", out=OUT, threads=8):
    """Empty cavity at dx: pencil f0/Q from the same pulse procedure and the mode-shape check."""
    p = out / "runs" / f"{CaseSpec(dx, 0.0, 'none', label='empty', height=height).name}.json"
    if p.exists():
        return json.loads(p.read_text(encoding="utf-8"))
    spec = CaseSpec(dx, 0.0, "none", label="empty", height=height, mode_shape=True)
    # mode shape DFT at the empty resonance needs f0 first: run once with a nominal value, then refine
    res = run_cavity(spec, threads, [2.9504e9], verbose=True)
    f0 = res["mode"]["f_hz"]
    if abs(f0 - 2.9504e9) > 2e6:
        res = run_cavity(spec, threads, [f0], verbose=True)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(res, default=float), encoding="utf-8")
    return res


def sweep(dx, jobs=4, threads=4, out=OUT, r0s=R0_CASES, ne_list=NE_CASES, extra=(10.0,), height="thin",
          model="lattice", skip_existing=True):
    emp = empty_reference(dx, out=out, threads=jobs * threads)
    f_eval = [emp["mode"]["f_hz"]]
    specs = []
    for s in case_specs(dx, r0s, ne_list, extra, height):
        s.model = model
        if not R0DX_MIN <= s.r0_m / dx <= R0DX_MAX:
            continue
        if skip_existing and (out / "runs" / f"{s.name}.json").exists():
            continue
        specs.append(s)
    print(f"dx {dx*1e3:g} mm: {len(specs)} runs, empty f0 {f_eval[0]/1e9:.5f} GHz", flush=True)
    if jobs <= 1:
        for s in specs:
            _worker(asdict(s), threads, f_eval, out)
        return
    with ProcessPoolExecutor(jobs) as ex:
        futs = [ex.submit(_worker, asdict(s), threads, f_eval, out) for s in specs]
        for fu in as_completed(futs):
            name, wall = fu.result()
            print(f"  done {name} ({wall:.0f} s)", flush=True)


def steady_check(dx, r0, ne, threads=8, out=OUT, n_tau=10.0):
    """CW drive at the empty f0 to steady state: rod power / forward power (direct, no DFT)."""
    numba.set_num_threads(threads)
    emp = empty_reference(dx, out=out, threads=threads)
    f0 = emp["mode"]["f_hz"]
    s, nu = rod_sigma(ne)
    rod = Rod(r0, "drude", ne, nu)
    fac = wire_factor(r0, dx)
    lay, _ = cavity_layout(dx, "thin", courant_for(fac))
    g = lay.grid
    fd = WireFDTD(g, lay.node_inside, lay.pml_lo, lay.pml_hi, f0)
    i0, k0 = g.nx // 2, lay.focus_k
    fd.install_rod(i0, k0, rod)
    port = TE10Port(lay, f0)
    port.amplitude = 1.0
    port.envelope = lambda t: 0.5 * (1 - math.cos(math.pi * min(t, 2e-9) / 2e-9))
    pulse = json.loads((out / "runs" / f"{CaseSpec(dx, r0, 'drude', ne, nu, label=f'ne{ne:.0e}'.replace('+', '')).name}.json").read_text())
    tau = pulse["mode"]["tau_s"]
    n = int((2e-9 + n_tau * tau) / g.dt)
    per = int(round(1 / (f0 * g.dt)))
    nav = 10 * per
    v = np.zeros(n); cur = np.zeros(n)
    w_mark = None
    for s_ in range(n):
        fd.step(port.inject)
        v[s_], cur[s_] = _plane_projection(fd.Ey, fd.Hx, lay.meas_k, port.profile, port.profile)
        if s_ == n - nav - 1:
            w_mark = float(fd.Py[i0, :, k0].sum()) * g.dx ** 3
    w_end = float(fd.Py[i0, :, k0].sum()) * g.dx ** 3
    p_rod = (w_end - w_mark) / (nav * g.dt)
    fw, bw = port.modal_amplitudes(v[-nav - per:], cur[-nav - per:], g.dt)
    p_f = float(port.power(np.array([fw.mean()]))[0])
    p_b = float(port.power(np.array([bw.mean()]))[0])
    res = dict(dx_m=dx, r0_m=r0, ne_m3=ne, f_drive_hz=f0, t_end_s=n * g.dt, tau_s=tau, n_tau=n_tau,
               p_rod_over_fwd=p_rod / p_f, p_back_over_fwd=p_b / p_f, closure=(p_f - p_b - p_rod) / p_f,
               dft_fraction=pulse["spectrum"]["rod_fraction"][pulse["spectrum"]["f_hz"].index(f0)])
    print(json.dumps(res), flush=True)
    (out / f"steady_dx{dx*1e3:g}_r{r0*1e6:.0f}_ne{ne:.0e}.json".replace("+", "")).write_text(json.dumps(res, indent=1))
    return res


def _load_runs(out=OUT):
    runs = {}
    for p in sorted((out / "runs").glob("*.json")):
        r = json.loads(p.read_text(encoding="utf-8"))
        runs[r["name"]] = r
    return runs


def compare(out=OUT, f_tol=0.10, q_tol=0.20):
    """Build the FDTD vs rod_drive vs TL table from stored runs."""
    runs = _load_runs(out)
    empties = {r["spec"]["dx_m"]: r for r in runs.values() if r["spec"]["kind"] == "none" and r["spec"]["height"] == "thin"}
    rows = []
    tls = {}
    for name, r in runs.items():
        sp = r["spec"]
        if sp["kind"] == "none":
            continue
        e = empties.get(sp["dx_m"])
        if e is None:
            continue
        f0, q0 = e["mode"]["f_hz"], e["mode"]["q"]
        rod = Rod(sp["r0_m"], sp["kind"], sp["ne_m3"], sp["nu_s"], sp["model"])
        R = rod.resistance(f0)
        lp = lumped_prediction(sp["r0_m"], R, f0, q0)
        key = sp["dx_m"]
        if key not in tls:
            tls[key] = TLModel(f0, q0, a=e["a_grid_m"], post_offset_m=e["rod_z_offset_m"])
        tl = tls[key].predict(sp["r0_m"], R, f0, lp["f_hz"] if R < 1e30 else f0)
        spec = r["spectrum"]
        k = int(np.argmin(np.abs(np.array(spec["f_hz"]) - f0)))
        m = r["mode"]
        df_fd = m["f_hz"] - f0
        lw = m["f_hz"] / (2 * m["q"])                                        # half-power half-width
        floor = 0.1 * lw
        df_ok = abs(lp["df_hz"] - df_fd) <= max(f_tol * abs(df_fd), floor)
        q_ok = abs(lp["q"] / m["q"] - 1) <= q_tol
        rows.append(dict(
            name=name, dx_mm=sp["dx_m"] * 1e3, r0_mm=sp["r0_m"] * 1e3, case=sp["label"], model=sp["model"],
            height=sp["height"], ne_m3=sp["ne_m3"], nu_s=sp["nu_s"], R_ohm=R, X_rod_drive_ohm=lp["x_ohm"],
            f0_empty_hz=f0, q0_empty=q0,
            f_fdtd_hz=m["f_hz"], q_fdtd=m["q"], df_fdtd_hz=df_fd, frac_fdtd=spec["rod_fraction"][k],
            closure_fdtd=spec["closure_residual"][k], f_drive_hz=spec["f_hz"][k],
            f_rod_drive_hz=lp["f_hz"], df_rod_drive_hz=lp["df_hz"], q_rod_drive=lp["q"], frac_rod_drive=lp["rod_fraction"],
            f_tl_hz=tl["f_hz"], df_tl_hz=tl["f_hz"] - f0, q_tl=tl["q"], frac_tl=tl["rod_fraction"],
            df_err_rel=(lp["df_hz"] - df_fd) / df_fd if df_fd else float("nan"),
            q_err_rel=lp["q"] / m["q"] - 1, df_tl_err_rel=(tl["f_hz"] - f0 - df_fd) / df_fd if df_fd else float("nan"),
            q_tl_err_rel=tl["q"] / m["q"] - 1, df_floor_hz=floor, df_pass=bool(df_ok), q_pass=bool(q_ok),
            final_envelope_rel=r["final_envelope_rel"], wire_factor=r["wire_factor"], wall_s=r["wall_s"],
            mode_x_f_hz=r["mode_x"]["f_hz"], mode_x_q=r["mode_x"]["q"]))
    rows.sort(key=lambda d: (d["height"], d["dx_mm"], d["r0_mm"], d["R_ohm"]))
    return rows, empties, tls


def report(out=OUT):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows, empties, tls = compare(out)
    val = [row for p in sorted(out.glob("validation_post*.json")) for row in json.loads(p.read_text())]
    steady = [json.loads(p.read_text()) for p in out.glob("steady_*.json")]
    hil = json.loads((out / "empty_checks.json").read_text()) if (out / "empty_checks.json").exists() else {}
    main_rows = [r for r in rows if r["height"] == "thin" and r["model"] == "lattice"]
    # convergence: same (r0, case) across dx
    conv = {}
    for r in main_rows:
        conv.setdefault((r["r0_mm"], r["case"]), {})[r["dx_mm"]] = dict(
            df=r["df_fdtd_hz"], q=r["q_fdtd"], frac=r["frac_fdtd"], q0=r["q0_empty"], df_tl=r["df_tl_hz"],
            q_tl=r["q_tl"], df_err_rd=r["df_err_rel"], q_err_rd=r["q_err_rel"])
    conv_rows = []
    for (r0, case), d in sorted(conv.items()):
        if len(d) < 2:
            continue
        dxs = sorted(d)
        fine, coarse = d[dxs[0]], d[dxs[-1]]
        conv_rows.append(dict(r0_mm=r0, case=case, dx_mm=dxs, df_mhz=[d[x]["df"] / 1e6 for x in dxs],
                              q=[d[x]["q"] for x in dxs], frac=[d[x]["frac"] for x in dxs],
                              q0_empty=[d[x]["q0"] for x in dxs],
                              df_fdtd_over_tl_minus1=[d[x]["df"] / d[x]["df_tl"] - 1 if d[x]["df_tl"] else None for x in dxs],
                              q_fdtd_over_tl_minus1=[d[x]["q"] / d[x]["q_tl"] - 1 for x in dxs],
                              rod_drive_df_err=[d[x]["df_err_rd"] for x in dxs], rod_drive_q_err=[d[x]["q_err_rd"] for x in dxs],
                              df_change_rel=(coarse["df"] - fine["df"]) / fine["df"] if fine["df"] else None,
                              q_change_rel=coarse["q"] / fine["q"] - 1))
    summary = dict(
        description=__doc__.split("\n\n")[0], criteria=dict(df_rel=0.10, q_rel=0.20,
                                                           df_floor="0.1 x FDTD half-power half-width"),
        empty=dict({f"dx{k*1e3:g}mm": dict(f0_hz=v["mode"]["f_hz"], q=v["mode"]["q"],
                                            mode_shape={kk: vv for kk, vv in v.get("mode_shape", {}).items()
                                                        if not isinstance(vv, list)},
                                            rod_z_offset_m=v["rod_z_offset_m"], a_grid_m=v["a_grid_m"],
                                            cavity_len_grid_m=v["cavity_len_grid_m"])
                    for k, v in empties.items()}, **hil),
        tl_fit={f"dx{k*1e3:g}mm": dict(b0=float(v.b0), d_eff_m=float(v.d)) for k, v in tls.items()},
        rows=rows, convergence=conv_rows, steady_checks=steady,
        validation_post=[{k: v for k, v in r.items()} for r in val],
        verdict=dict(
            n=len(main_rows),
            df_pass=sum(r["df_pass"] for r in main_rows), q_pass=sum(r["q_pass"] for r in main_rows),
            both_pass=sum(r["df_pass"] and r["q_pass"] for r in main_rows)),
        software=dict(numpy=np.__version__, numba=numba.__version__, python=__import__("sys").version.split()[0]))
    (out / "rod_cavity_check.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")

    # ---- figure
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    cols = dict(zip(sorted({r["r0_mm"] for r in main_rows}), ["C0", "C1", "C2", "C3", "C4"]))
    dx_mark = {2.0: "D", 1.0: "o", 0.5: "s", 0.25: "^"}
    Rplot = lambda R: max(R, 1.0)   # PEC drawn at 1 ohm
    for r in main_rows:
        c = cols[r["r0_mm"]]
        mk = dx_mark.get(r["dx_mm"], "x")
        ax[0, 0].plot(Rplot(r["R_ohm"]), r["df_fdtd_hz"] / 1e6, mk, color=c, mfc="none" if r["dx_mm"] != 0.5 else c, ms=6)
        if r["dx_mm"] == 1.0:                      # Q and fraction depend on the dx-dependent iris Q0
            ax[0, 1].plot(Rplot(r["R_ohm"]), r["q_fdtd"], mk, color=c, mfc="none", ms=6)
            ax[1, 0].plot(Rplot(r["R_ohm"]), r["frac_fdtd"], mk, color=c, mfc="none", ms=6)
    kref = 1e-3 if 1e-3 in empties else min(empties)
    f0, q0 = empties[kref]["mode"]["f_hz"], empties[kref]["mode"]["q"]
    Rg = np.logspace(0, 6.2, 160)
    for r0 in cols:
        lp = [lumped_prediction(r0 * 1e-3, R, f0, q0) for R in Rg]
        tlm = tls.get(kref)
        ax[0, 0].plot(Rg, [p["df_hz"] / 1e6 for p in lp], "-", color=cols[r0], lw=1.2)
        ax[0, 1].plot(Rg, [p["q"] for p in lp], "-", color=cols[r0], lw=1.2)
        ax[1, 0].plot(Rg, [p["rod_fraction"] for p in lp], "-", color=cols[r0], lw=1.2, label=f"r0 = {r0:g} mm")
        if tlm is not None:
            tp, guess = [], f0
            for R, p in zip(Rg[::-1], lp[::-1]):
                try:
                    q = tlm.predict(r0 * 1e-3, R, f0, guess)
                    guess = q["f_hz"]
                except RuntimeError:
                    q = dict(f_hz=np.nan, q=np.nan, rod_fraction=np.nan)
                tp.append(q)
            tp = tp[::-1]
            ax[0, 0].plot(Rg, [(p["f_hz"] - f0) / 1e6 for p in tp], "--", color=cols[r0], lw=1)
            ax[0, 1].plot(Rg, [p["q"] for p in tp], "--", color=cols[r0], lw=1)
            ax[1, 0].plot(Rg, [p["rod_fraction"] for p in tp], "--", color=cols[r0], lw=1)
    for a_, yl in ((ax[0, 0], "resonance shift f - f0 (MHz)"), (ax[0, 1], "loaded Q"),
                   (ax[1, 0], "P_rod / P_fwd at drive f = f0(empty)")):
        a_.set_xscale("log"); a_.set_xlabel("rod resistance R = b / (Re sigma pi r0^2)  (ohm; PEC plotted at 1)")
        a_.set_ylabel(yl); a_.grid(alpha=0.3)
    ax[0, 1].set_yscale("log")
    ax[0, 0].set_title("Shift: FDTD (markers) vs rod_drive lumped (solid) vs TL proposal (dashed)", fontsize=9)
    ax[0, 1].set_title("Loaded Q (FDTD matrix-pencil pole), dx = 1 mm (calibrated iris)", fontsize=9)
    ax[1, 0].set_title("P_rod / P_fwd for steady drive at the empty-cavity f0, dx = 1 mm", fontsize=9)
    h = [plt.Line2D([], [], color=c, label=f"r0 = {r0:g} mm") for r0, c in cols.items()]
    h += [plt.Line2D([], [], color="k", marker=m, ls="none", mfc="none" if d != 0.5 else "k", label=f"FDTD dx = {d:g} mm")
          for d, m in dx_mark.items() if any(r["dx_mm"] == d for r in main_rows)]
    h += [plt.Line2D([], [], color="k", ls="-", label="rod_drive lumped"), plt.Line2D([], [], color="k", ls="--", label="TL (proposed)")]
    ax[1, 0].legend(handles=h, fontsize=7, loc="upper left")
    # validation panel
    a_ = ax[1, 1]
    vm = dict(none="k", tu="C3", lattice="C0")
    for r in val:
        if r["kind"] != "pec":
            continue
        k = r["f_hz"].index(2.95e9)
        a_.plot(r["r0_over_dx"], 100 * r["x_err_rel"][k], dx_mark.get(r["dx_m"] * 1e3, "x"), color=vm[r["model"]],
                mfc="none" if r["dx_m"] != 0.5e-3 else vm[r["model"]])
    a_.axhline(0, color="0.5", lw=0.8)
    a_.set_xlabel("r0 / dx"); a_.set_ylabel("post reactance error vs Marcuvitz at 2.95 GHz (%)")
    a_.set_title("Subcell validation: PEC post in matched WR-284 (FDTD)\nblue lattice-corrected, red Umashankar-Taflove, black plain edge (r_eq = 0.1985 dx)", fontsize=9)
    a_.grid(alpha=0.3)
    fig.suptitle("Bridged rod in the iris-coupled TE101 cavity: 3D FDTD (y-invariant, PEC walls) vs rod_drive lumped coupling", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "rod_cavity_check.png", dpi=130)
    plt.close(fig)
    print_table(rows)
    return summary


def print_table(rows):
    print(f"{'dx':>4} {'r0':>5} {'case':>8} {'R ohm':>9} | {'df FDTD':>8} {'df rd':>8} {'err':>6} | "
          f"{'Q FDTD':>7} {'Q rd':>7} {'err':>6} | {'Pfr FDTD':>8} {'Pfr rd':>7} | {'df TL':>7} {'Q TL':>7} {'Pfr TL':>7} | pass")
    for r in rows:
        print(f"{r['dx_mm']:4.2g} {r['r0_mm']:5.2f} {r['case']:>8} {r['R_ohm']:9.3g} | {r['df_fdtd_hz']/1e6:8.2f} "
              f"{r['df_rod_drive_hz']/1e6:8.2f} {100*r['df_err_rel']:+5.0f}% | {r['q_fdtd']:7.1f} {r['q_rod_drive']:7.1f} "
              f"{100*r['q_err_rel']:+5.0f}% | {r['frac_fdtd']:8.3f} {r['frac_rod_drive']:7.3f} | {r['df_tl_hz']/1e6:7.1f} "
              f"{r['q_tl']:7.1f} {r['frac_tl']:7.3f} | {'f' if r['df_pass'] else '-'}{'Q' if r['q_pass'] else '-'} "
              f"{r['height'] if r['height'] != 'thin' else ''}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="rod_cavity_check", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="PEC/Drude post in a matched guide vs Marcuvitz (subcell model)")
    v.add_argument("--dx-mm", type=float, nargs="+", default=[1.0, 0.5, 0.25])
    v.add_argument("--r0-mm", type=float, nargs="+", default=[0.05, 0.1, 0.15, 0.2, 0.3])
    v.add_argument("--tag", default="", help="output file suffix: validation_post<tag>.json")
    v.add_argument("--threads", type=int, default=8)
    e = sub.add_parser("empty", help="empty cavity: pencil f0/Q + mode shape; Hilbert ringdown as coupled3d")
    e.add_argument("--dx-mm", type=float, nargs="+", default=[1.0, 0.5])
    e.add_argument("--full", action="store_true", help="also the full-height Hilbert ringdown at 1 mm")
    e.add_argument("--threads", type=int, default=8)
    c = sub.add_parser("case", help="single cavity run")
    c.add_argument("--dx-mm", type=float, default=1.0)
    c.add_argument("--r0-mm", type=float, default=0.2)
    c.add_argument("--ne", type=float, default=0.0, help="Drude rod density (0: PEC)")
    c.add_argument("--model", default="lattice", choices=["lattice", "tu", "none"])
    c.add_argument("--full", action="store_true", help="full 34 mm height instead of the y-invariant 2-cell box")
    c.add_argument("--threads", type=int, default=8)
    s = sub.add_parser("sweep", help="r0 x conductance matrix at one dx (parallel processes)")
    s.add_argument("--dx-mm", type=float, default=1.0)
    s.add_argument("--r0-mm", type=float, nargs="+", default=[r * 1e3 for r in R0_CASES])
    s.add_argument("--ne", type=float, nargs="+", default=list(NE_CASES))
    s.add_argument("--extra-r-ohm", type=float, nargs="*", default=[10.0])
    s.add_argument("--no-pec", action="store_true")
    s.add_argument("--jobs", type=int, default=4)
    s.add_argument("--threads", type=int, default=4)
    st = sub.add_parser("steady", help="CW steady-drive check of the DFT absorbed fraction")
    st.add_argument("--dx-mm", type=float, default=1.0)
    st.add_argument("--r0-mm", type=float, default=0.2)
    st.add_argument("--ne", type=float, default=1e23)
    st.add_argument("--threads", type=int, default=8)
    sub.add_parser("report", help="comparison table, JSON summary and PNG")
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.cmd == "validate":
        validate([d * 1e-3 for d in a.dx_mm], [r * 1e-3 for r in a.r0_mm], threads=a.threads, tag=a.tag)
    elif a.cmd == "empty":
        chk = {}
        for d in a.dx_mm:
            r = empty_reference(d * 1e-3, threads=a.threads)
            print(json.dumps({k: v for k, v in r["mode_shape"].items() if not isinstance(v, list)}, indent=1))
            chk[f"hilbert_thin_dx{d:g}mm"] = hilbert_ringdown(d * 1e-3, "thin", threads=a.threads)
            print(chk[f"hilbert_thin_dx{d:g}mm"], flush=True)
        if a.full:
            chk["hilbert_full_dx1mm"] = hilbert_ringdown(1e-3, "full", threads=a.threads)
            print(chk["hilbert_full_dx1mm"])
        chk["prefill_full_dx1mm_profiles"] = prefill_profiles()
        p = OUT / "empty_checks.json"
        old = json.loads(p.read_text()) if p.exists() else {}
        old.update(chk)
        p.write_text(json.dumps(old, indent=1))
    elif a.cmd == "case":
        if a.ne > 0:
            s_, nu = rod_sigma(a.ne)
            spec = CaseSpec(a.dx_mm * 1e-3, a.r0_mm * 1e-3, "drude", a.ne, nu, label=f"ne{a.ne:.0e}".replace("+", ""),
                            model=a.model, height="full" if a.full else "thin")
        else:
            spec = CaseSpec(a.dx_mm * 1e-3, a.r0_mm * 1e-3, "pec", label="pec", model=a.model,
                            height="full" if a.full else "thin")
        emp = empty_reference(a.dx_mm * 1e-3, threads=a.threads)
        _worker(asdict(spec), a.threads, [emp["mode"]["f_hz"]], OUT)
    elif a.cmd == "sweep":
        r0s = [r * 1e-3 for r in a.r0_mm]
        if a.no_pec:
            raise SystemExit("--no-pec is not supported; PEC is the reference")
        sweep(a.dx_mm * 1e-3, a.jobs, a.threads, r0s=r0s, ne_list=a.ne, extra=a.extra_r_ohm)
    elif a.cmd == "steady":
        steady_check(a.dx_mm * 1e-3, a.r0_mm * 1e-3, a.ne, a.threads)
    elif a.cmd == "report":
        report()


if __name__ == "__main__":
    main()
