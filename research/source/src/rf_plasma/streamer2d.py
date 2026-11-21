"""Axisymmetric (r, z) drift-diffusion-Poisson channel-growth side model.

Standalone model (NOT coupled to the 3D FDTD solver) of a laser-seeded plasma
filament in helium driven by a uniform RF field E_z(t) = E0 sin(2 pi f t) parallel
to the channel axis ("laser parallel to RF field" configuration, so axisymmetry
holds).

Equations (SI):
    dn_e/dt + div(G_e) = (k_iz(|E|/N) N - alpha n_i) n_e,  G_e = -mu_e n_e E - D_e grad n_e
    dn_i/dt + div(G_i) = (k_iz(|E|/N) N - alpha n_e) n_e,  G_i = +mu_i n_i E - D_i grad n_i
    -div(eps0 grad phi) = e (n_i - n_e),   E = -grad phi + E0 sin(2 pi f t) z_hat

Local-field approximation: mu_e N, k_iz and mean energy from the BOLOS/LXCat He
table (cases/swarm/He_boltzmann.csv), log-interpolated in E/N.  D_e = mu_e Te with
Te = (2/3) mean energy (Einstein relation).  Ions: He2+ with K0 = 16.7 cm^2/Vs
(ASSUMPTION, not re-verified) scaled by N0/N, D_i = mu_i k Tg / e.

Numerics: finite volumes on a nonuniform axisymmetric grid, Scharfetter-Gummel
fluxes, forward-Euler transport under the exact positivity (CFL+diffusion) bound,
and a semi-implicit Poisson step: the charge-flux linearisation
J^{n+1} = J^n + sigma_f (E^{n+1} - E^n), with sigma_f the exact derivative of the SG
charge flux with respect to the face field, is inserted into Gauss's law. This
removes the explicit dielectric-relaxation limit (tau_d ~ 1e-14 s at 1e23 m^-3)
without changing the conservative density update.  Production mode ("frozen"):
Lambda_f = h dt max(sigma_f, sigma_env) with h = 1.5 and sigma_env an envelope built from
the table's maximum mu_e and max(n_L, n_R); the SPD operator -div((eps0 + Lambda) grad)
is LU-factorised (SuperLU) and reused for up to 100 steps, refactorised earlier if
dt sigma_f > 1.5 (eps0 + Lambda_f) anywhere.  Over-estimating Lambda only adds damping
of the relaxation mode (0-D amplification 1 - x/(1+y), stable for x <= 2(1+y)); steady
and quasi-steady states are unchanged.  Verified against "exact" mode (refactorised every
step) on a coarse grid: L2 difference 0.1 % in |E| and 0.5 % in n_e at 1 ns.  The vacuum
Laplacian is factorised once per grid (initial solve and tests).

CLI:  python -m rf_plasma.streamer2d list | run CASE [CASE ...] | run --all | figures
"""

from __future__ import annotations

import os

if __name__ == "__main__":
    # CLI only (python -m rf_plasma.streamer2d): limit BLAS/OpenMP pools before numpy is
    # imported so parallel case runs stay at 1 thread each.  Not done on plain import, so
    # importing this module never alters the threading environment of other code/tests.
    # NUMBA_NUM_THREADS is deliberately not touched (kernels here are serial njit).
    for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(_v, "1")

import argparse
import json
import math
import platform
import sys
import time
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import numba
from numba import njit

E_CH = 1.602176634e-19
EPS0 = 8.8541878128e-12
KB = 1.380649e-23
N_LOSCHMIDT = 2.686780111e25  # m^-3 at 273.15 K, 101.325 kPa (reference for K0)
TD = 1e-21  # V m^2

ROOT = Path(__file__).resolve().parents[2]
SWARM_FILE = ROOT / "cases" / "swarm" / "He_boltzmann.csv"
RESULTS = ROOT / "results" / "streamer2d"


# --------------------------------------------------------------------------- config
@dataclass
class Config:
    name: str = "test"
    pressure_Pa: float = 0.2e6
    T_gas_K: float = 300.0
    E0_V_m: float = 2.23e6          # amplitude of E_z(t) = E0 sin(2 pi f t)
    freq_Hz: float = 2.95e9
    n0_m3: float = 1e23
    a_m: float = 10e-6
    L_m: float = 1e-3
    n_bg_m3: float = 1e9
    alpha_rec_m3_s: float = 1e-14   # UNVERIFIED two-body He2+ + e coefficient
    K0_ion_m2_Vs: float = 16.7e-4   # He2+ in He, ASSUMPTION
    ionization: bool = True
    swarm_file: str = str(SWARM_FILE)
    # grid
    r_max_m: float = 3e-3
    dr_min_m: float = 1e-6
    r_uniform_m: float = 20e-6
    dr_growth: float = 1.05
    dr_mid_m: float = 8e-6
    r_mid_m: float = 0.5e-3
    dr_growth_out: float = 1.15
    z_half_m: float = 4e-3
    dz_m: float = 5e-6
    # time
    t_end_s: float = 10e-9
    safety: float = 0.45            # fraction of the forward-Euler positivity bound
    src_dt_frac: float = 0.2        # max (nu_iz + alpha n) dt
    steps_per_rf_period: int = 40
    diag_dt_s: float = 20e-12
    snapshot_times_s: tuple = (0.0, 1e-9, 2e-9, 3e-9, 5e-9, 7.5e-9, 10e-9, 15e-9, 20e-9)
    stop_margin_m: float = 0.25e-3  # stop when the 1e20 contour gets this close to a boundary
    wall_limit_s: float = 4 * 3600.0
    # Poisson (semi-implicit): "frozen" = Lambda_f = h dt sigma_f held fixed between LU
    # refactorisations; "exact" = Lambda_f = dt sigma_f refactorised every step (verification);
    # "explicit" = no semi-implicit term (only valid for dt << tau_d).
    si_mode: str = "frozen"
    si_mu_envelope: str = "max"      # "max": Lambda uses the table's maximum mu_e; "local": current mu_e
    si_headroom: float = 1.5
    si_refresh_steps: int = 100
    si_violation: float = 1.5
    thresholds_m3: tuple = (1e19, 1e20, 1e21)

    @property
    def N(self) -> float:
        return self.pressure_Pa / (KB * self.T_gas_K)

    @property
    def EN_Td(self) -> float:
        return self.E0_V_m / self.N / TD


def _cases() -> dict:
    N02 = 0.2e6 / (KB * 300.0)
    en_hi = 2.23e6 / N02  # E/N of 22.3 kV/cm at 0.2 MPa (~46 Td)
    out = {}

    def add(name, p, E0, **kw):
        out[name] = Config(name=name, pressure_Pa=p, E0_V_m=E0, **kw)

    add("He0p2MPa_E16p6", 0.2e6, 1.66e6)
    add("He0p2MPa_E22p3", 0.2e6, 2.23e6)
    for p, tag in ((0.446e6, "He0p446MPa_46Td"), (1.0e6, "He1MPa_46Td"), (1.825e6, "He1p825MPa_46Td")):
        add(tag, p, en_hi * p / (KB * 300.0))
    add("He1MPa_E80", 1.0e6, 8.0e6)
    # optional variants
    add("He0p2MPa_E22p3_n1e24_a20", 0.2e6, 2.23e6, n0_m3=1e24, a_m=20e-6)
    add("He0p2MPa_E22p3_bg1e12", 0.2e6, 2.23e6, n_bg_m3=1e12)
    return out


CASES = _cases()


# --------------------------------------------------------------------------- swarm
def read_swarm(path) -> dict:
    lines = [ln for ln in Path(path).read_text(encoding="utf-8").splitlines() if ln.strip() and not ln.startswith("#")]
    head = lines[0].split(",")
    arr = np.array([[float(x) for x in ln.split(",")] for ln in lines[1:]])
    return {h: arr[:, i] for i, h in enumerate(head)}


class Swarm:
    """LFA transport on a fine uniform ln(E/N) table (log-log interpolation of the CSV)."""

    def __init__(self, path, N: float, n_fine: int = 4001):
        d = read_swarm(path)
        en = d["E_N_Td"]
        x = np.log(en)
        self.x0, self.x1 = float(x[0]), float(x[-1])
        xf = np.linspace(self.x0, self.x1, n_fine)
        self.xf = xf
        self.mu = np.exp(np.interp(xf, x, np.log(d["muN_1_Vms"]))) / N
        te = (2.0 / 3.0) * d["mean_energy_eV"]
        self.Te = np.exp(np.interp(xf, x, np.log(te)))
        # k_iz: the table holds solver noise (|k| < 1e-34, both signs) at low E/N.
        # Use the contiguous high-E/N branch with k_iz > 1e-32 m^3/s; zero below it.
        k = d["kiz_m3_s"]
        good = k > 1e-32
        i0 = len(k) - 1
        while i0 > 0 and good[i0 - 1]:
            i0 -= 1
        self.kiz_EN_min_Td = float(en[i0])
        kf = np.exp(np.interp(xf, x[i0:], np.log(k[i0:])))
        kf[xf < x[i0]] = 0.0
        self.nu = kf * N
        self.D = self.mu * self.Te
        self.N = N

    def at(self, en_td):
        """(mu_e, D_e, nu_iz) at reduced field en_td (Td); clamped to the table range."""
        xx = np.clip(np.log(np.maximum(en_td, 1e-30)), self.x0, self.x1)
        return np.interp(xx, self.xf, self.mu), np.interp(xx, self.xf, self.D), np.interp(xx, self.xf, self.nu)


# --------------------------------------------------------------------------- grid
def radial_faces(cfg: Config) -> np.ndarray:
    f = [0.0]
    dr = cfg.dr_min_m
    while f[-1] < cfg.r_max_m - 1e-15:
        r = f[-1]
        if r >= cfg.r_mid_m:
            dr = max(dr, cfg.dr_mid_m) * cfg.dr_growth_out
        elif r >= cfg.r_uniform_m:
            dr = min(dr * cfg.dr_growth, cfg.dr_mid_m)
        f.append(r + dr)
    f[-1] = max(f[-1], cfg.r_max_m)
    return np.array(f)


class Grid:
    def __init__(self, rf: np.ndarray, zf: np.ndarray):
        self.rf, self.zf = rf, zf
        self.nr, self.nz = len(rf) - 1, len(zf) - 1
        self.rc = 0.5 * (rf[1:] + rf[:-1])
        self.zc = 0.5 * (zf[1:] + zf[:-1])
        self.dzc = np.diff(zf)
        ring = np.pi * (rf[1:] ** 2 - rf[:-1] ** 2)
        self.V = ring[:, None] * self.dzc[None, :]
        self.Ar = 2 * np.pi * rf[:, None] * self.dzc[None, :]       # (nr+1, nz)
        self.Az = np.repeat(ring[:, None], self.nz + 1, axis=1)     # (nr, nz+1)
        hr = np.empty(self.nr + 1)
        hr[0] = rf[1] - 0.0
        hr[1:-1] = np.diff(self.rc)
        hr[-1] = rf[-1] - self.rc[-1]
        hz = np.empty(self.nz + 1)
        hz[0] = self.zc[0] - zf[0]
        hz[1:-1] = np.diff(self.zc)
        hz[-1] = zf[-1] - self.zc[-1]
        self.hr, self.hz = hr, hz
        self.cr = self.Ar / hr[:, None]       # geometric conductance A/h of r faces
        self.cz = self.Az / hz[None, :]
        self.cr[0, :] = 0.0                    # axis: no flux
        self.indptr, self.indices = _csr_pattern(self.nr, self.nz)
        self.shape = (self.nr, self.nz)

    @classmethod
    def from_config(cls, cfg: Config):
        nz = int(round(2 * cfg.z_half_m / cfg.dz_m))
        return cls(radial_faces(cfg), np.linspace(-cfg.z_half_m, cfg.z_half_m, nz + 1))


def sr_dummy(g: Grid):
    return np.zeros((g.nr + 1, g.nz))


def sz_dummy(g: Grid):
    return np.zeros((g.nr, g.nz + 1))


@njit(cache=True)
def _csr_pattern(nr, nz):
    n = nr * nz
    indptr = np.zeros(n + 1, np.int64)
    indices = np.empty(5 * n, np.int64)
    p = 0
    for i in range(nr):
        for j in range(nz):
            k = i * nz + j
            if i > 0:
                indices[p] = k - nz; p += 1
            if j > 0:
                indices[p] = k - 1; p += 1
            indices[p] = k; p += 1
            if j < nz - 1:
                indices[p] = k + 1; p += 1
            if i < nr - 1:
                indices[p] = k + nz; p += 1
            indptr[k + 1] = p
    return indptr, indices[:p]


@njit(cache=True)
def _csr_data(nr, nz, wr, wz, nnz):
    """SPD operator sum_f w_f (phi_k - phi_nb); Dirichlet phi=0 on outer r and both z faces."""
    data = np.empty(nnz)
    p = 0
    for i in range(nr):
        for j in range(nz):
            d = wr[i, j] + wr[i + 1, j] + wz[i, j] + wz[i, j + 1]
            if i > 0:
                data[p] = -wr[i, j]; p += 1
            if j > 0:
                data[p] = -wz[i, j]; p += 1
            data[p] = d; p += 1
            if j < nz - 1:
                data[p] = -wz[i, j + 1]; p += 1
            if i < nr - 1:
                data[p] = -wr[i + 1, j]; p += 1
    return data


def assemble(g: Grid, wr, wz) -> sp.csr_matrix:
    n = g.nr * g.nz
    data = _csr_data(g.nr, g.nz, wr, wz, len(g.indices))
    return sp.csr_matrix((data, g.indices, g.indptr), shape=(n, n))


def factor(A: sp.spmatrix):
    return spla.splu(A.tocsc(), permc_spec="MMD_AT_PLUS_A")


def apply_op(g: Grid, wr, wz, phi):
    """Same operator as assemble(), applied matrix-free (used for the explicit phi^n term)."""
    out = (wr[1:, :] + wr[:-1, :] + wz[:, 1:] + wz[:, :-1]) * phi
    out[1:, :] -= wr[1:-1, :] * phi[:-1, :]
    out[:-1, :] -= wr[1:-1, :] * phi[1:, :]
    out[:, 1:] -= wz[:, 1:-1] * phi[:, :-1]
    out[:, :-1] -= wz[:, 1:-1] * phi[:, 1:]
    return out


def face_fields(g: Grid, phi, Eapp):
    """Face-normal fields from the space-charge potential (phi=0 on far boundaries) + applied E_z."""
    Er = np.zeros((g.nr + 1, g.nz))
    Er[1:-1, :] = -(phi[1:, :] - phi[:-1, :]) / g.hr[1:-1, None]
    Er[-1, :] = phi[-1, :] / g.hr[-1]
    Ez = np.empty((g.nr, g.nz + 1))
    Ez[:, 1:-1] = -(phi[:, 1:] - phi[:, :-1]) / g.hz[None, 1:-1]
    Ez[:, 0] = -(phi[:, 0] - 0.0) / g.hz[0]
    Ez[:, -1] = phi[:, -1] / g.hz[-1]
    Ez += Eapp
    return Er, Ez


def cell_field(Er, Ez):
    erc = 0.5 * (Er[1:, :] + Er[:-1, :])
    ezc = 0.5 * (Ez[:, 1:] + Ez[:, :-1])
    return erc, ezc, np.hypot(erc, ezc)


# --------------------------------------------------------------------------- SG kernels
@njit(cache=True, inline="always")
def _bern_pair(x):
    """Return B(x), B'(x) for the Bernoulli function B(x) = x/(exp(x)-1) with one expm1.

    Identities used by callers: B(-x) = x + B(x),  B'(-x) = -1 - B'(x).
    """
    if abs(x) < 1e-4:
        return 1.0 - 0.5 * x + x * x / 12.0, -0.5 + x / 6.0
    if x > 700.0:
        return 0.0, 0.0
    if x < -700.0:
        return -x, -1.0
    em = math.expm1(x)
    b = x / em
    return b, (em - x * (em + 1.0)) / (em * em)


@njit(cache=True)
def sg_fluxes(n, mu_c, D_c, sign, Er, Ez, hr, hz, V, Ar, Az, want_deriv):
    """SG fluxes on all faces, d(flux)/d(v) on faces, and per-cell outflow rate (1/s).

    Flux on an interior face (left cell L, right cell R, spacing h, Pe = v h / D):
        F = D/h [B(-Pe) n_L - B(Pe) n_R],    dF/dv = -B'(-Pe) n_L - B'(Pe) n_R >= 0.
    sign = -1 electrons (v = -mu E), +1 ions. Outer faces: outflow-only drift (no inflow,
    zero diffusive flux). Axis: zero flux. 'out' is the forward-Euler positivity rate.
    """
    nr, nz = n.shape
    Fr = np.zeros((nr + 1, nz))
    Fz = np.zeros((nr, nz + 1))
    Gr = np.zeros((nr + 1, nz))
    Gz = np.zeros((nr, nz + 1))
    out = np.zeros((nr, nz))
    for i in range(1, nr + 1):
        for j in range(nz):
            if i < nr:
                mu = 0.5 * (mu_c[i - 1, j] + mu_c[i, j])
                D = 0.5 * (D_c[i - 1, j] + D_c[i, j])
                v = sign * mu * Er[i, j]
                h = hr[i]
                pe = v * h / D
                bp, dbp = _bern_pair(pe)
                bm = pe + bp
                Fr[i, j] = D / h * (bm * n[i - 1, j] - bp * n[i, j])
                if want_deriv:
                    Gr[i, j] = (1.0 + dbp) * n[i - 1, j] - dbp * n[i, j]
                    out[i - 1, j] += Ar[i, j] * D / h * bm / V[i - 1, j]
                    out[i, j] += Ar[i, j] * D / h * bp / V[i, j]
            else:
                v = sign * mu_c[nr - 1, j] * Er[nr, j]
                if v > 0.0:
                    Fr[nr, j] = v * n[nr - 1, j]
                    out[nr - 1, j] += Ar[nr, j] * v / V[nr - 1, j]
                Gr[nr, j] = n[nr - 1, j]
    for i in range(nr):
        for j in range(nz + 1):
            if 0 < j < nz:
                mu = 0.5 * (mu_c[i, j - 1] + mu_c[i, j])
                D = 0.5 * (D_c[i, j - 1] + D_c[i, j])
                v = sign * mu * Ez[i, j]
                h = hz[j]
                pe = v * h / D
                bp, dbp = _bern_pair(pe)
                bm = pe + bp
                Fz[i, j] = D / h * (bm * n[i, j - 1] - bp * n[i, j])
                if want_deriv:
                    Gz[i, j] = (1.0 + dbp) * n[i, j - 1] - dbp * n[i, j]
                    out[i, j - 1] += Az[i, j] * D / h * bm / V[i, j - 1]
                    out[i, j] += Az[i, j] * D / h * bp / V[i, j]
            elif j == 0:
                v = sign * mu_c[i, 0] * Ez[i, 0]
                if v < 0.0:
                    Fz[i, 0] = v * n[i, 0]
                    out[i, 0] += -Az[i, 0] * v / V[i, 0]
                Gz[i, 0] = n[i, 0]
            else:
                v = sign * mu_c[i, nz - 1] * Ez[i, nz]
                if v > 0.0:
                    Fz[i, nz] = v * n[i, nz - 1]
                    out[i, nz - 1] += Az[i, nz] * v / V[i, nz - 1]
                Gz[i, nz] = n[i, nz - 1]
    return Fr, Fz, Gr, Gz, out


def divergence(g: Grid, Fr, Fz):
    """Integrated outward flux per cell (same units as flux * area)."""
    return (g.Ar[1:, :] * Fr[1:, :] - g.Ar[:-1, :] * Fr[:-1, :]
            + g.Az[:, 1:] * Fz[:, 1:] - g.Az[:, :-1] * Fz[:, :-1])


def boundary_outflow(g: Grid, Fr, Fz):
    """Particles per second leaving through outer r and both z boundaries."""
    return float(np.sum(g.Ar[-1, :] * Fr[-1, :]) + np.sum(g.Az[:, -1] * Fz[:, -1])
                 - np.sum(g.Az[:, 0] * Fz[:, 0]))


# fused kernels (same arithmetic as face_fields/cell_field/Swarm.at/divergence, fewer memory passes)
@njit(cache=True)
def k_fields(phi, Eapp, hr, hz, inv_N_td, x0, x1, dxf, mu_t, D_t, nu_t):
    nr, nz = phi.shape
    Er = np.zeros((nr + 1, nz))
    Ez = np.empty((nr, nz + 1))
    for i in range(1, nr):
        for j in range(nz):
            Er[i, j] = -(phi[i, j] - phi[i - 1, j]) / hr[i]
    for j in range(nz):
        Er[nr, j] = phi[nr - 1, j] / hr[nr]
    for i in range(nr):
        Ez[i, 0] = -phi[i, 0] / hz[0] + Eapp
        for j in range(1, nz):
            Ez[i, j] = -(phi[i, j] - phi[i, j - 1]) / hz[j] + Eapp
        Ez[i, nz] = phi[i, nz - 1] / hz[nz] + Eapp
    erc = np.empty((nr, nz)); ezc = np.empty((nr, nz)); Eabs = np.empty((nr, nz))
    mu = np.empty((nr, nz)); D = np.empty((nr, nz)); nu = np.empty((nr, nz))
    nt = mu_t.shape[0]
    for i in range(nr):
        for j in range(nz):
            a = 0.5 * (Er[i, j] + Er[i + 1, j])
            b = 0.5 * (Ez[i, j] + Ez[i, j + 1])
            e = math.sqrt(a * a + b * b)
            erc[i, j] = a; ezc[i, j] = b; Eabs[i, j] = e
            x = math.log(max(e * inv_N_td, 1e-30))
            x = min(max(x, x0), x1)
            s = (x - x0) / dxf
            k = min(int(s), nt - 2)
            f = s - k
            mu[i, j] = mu_t[k] + f * (mu_t[k + 1] - mu_t[k])
            D[i, j] = D_t[k] + f * (D_t[k + 1] - D_t[k])
            nu[i, j] = nu_t[k] + f * (nu_t[k + 1] - nu_t[k])
    return Er, Ez, erc, ezc, Eabs, mu, D, nu


@njit(cache=True)
def k_sigma_env(ne, ni, mu_e, mu_i):
    """Upwind-independent face conductivity e(mu_e,f max(ne_L,ne_R) + mu_i max(ni_L,ni_R)).

    Used only to set the frozen semi-implicit weight Lambda_f, so that a reversal of the
    face field (which switches the SG upwind cell) does not invalidate the factorisation.
    """
    nr, nz = ne.shape
    sr = np.zeros((nr + 1, nz))
    sz = np.zeros((nr, nz + 1))
    for i in range(1, nr + 1):
        for j in range(nz):
            if i < nr:
                m = 0.5 * (mu_e[i - 1, j] + mu_e[i, j])
                sr[i, j] = E_CH * (m * max(ne[i - 1, j], ne[i, j]) + mu_i * max(ni[i - 1, j], ni[i, j]))
            else:
                sr[i, j] = E_CH * (mu_e[nr - 1, j] * ne[nr - 1, j] + mu_i * ni[nr - 1, j])
    for i in range(nr):
        for j in range(nz + 1):
            if j == 0:
                sz[i, j] = E_CH * (mu_e[i, 0] * ne[i, 0] + mu_i * ni[i, 0])
            elif j == nz:
                sz[i, j] = E_CH * (mu_e[i, nz - 1] * ne[i, nz - 1] + mu_i * ni[i, nz - 1])
            else:
                m = 0.5 * (mu_e[i, j - 1] + mu_e[i, j])
                sz[i, j] = E_CH * (m * max(ne[i, j - 1], ne[i, j]) + mu_i * max(ni[i, j - 1], ni[i, j]))
    return sr, sz


@njit(cache=True)
def k_sigma(mu_e, mu_i, Gre, Gze, Gri, Gzi, dt, lam_r, lam_z, have_lam):
    """Face sigma_f = e(mu_e G_e + mu_i G_i) and max_f dt sigma_f / (eps0 + Lambda_f)."""
    nr, nz = mu_e.shape
    sr = np.zeros((nr + 1, nz))
    sz = np.zeros((nr, nz + 1))
    worst = 0.0
    for i in range(1, nr + 1):
        for j in range(nz):
            m = 0.5 * (mu_e[i - 1, j] + mu_e[i, j]) if i < nr else mu_e[nr - 1, j]
            s = E_CH * (m * Gre[i, j] + mu_i * Gri[i, j])
            sr[i, j] = s
            if have_lam:
                worst = max(worst, dt * s / (EPS0 + lam_r[i, j]))
    for i in range(nr):
        for j in range(nz + 1):
            if j == 0:
                m = mu_e[i, 0]
            elif j == nz:
                m = mu_e[i, nz - 1]
            else:
                m = 0.5 * (mu_e[i, j - 1] + mu_e[i, j])
            s = E_CH * (m * Gze[i, j] + mu_i * Gzi[i, j])
            sz[i, j] = s
            if have_lam:
                worst = max(worst, dt * s / (EPS0 + lam_z[i, j]))
    return sr, sz, worst


@njit(cache=True)
def k_rhs(ne, ni, V, Ar, Az, cr, cz, Fre, Fze, Fri, Fzi, lam_r, lam_z, phi, dt, dEapp):
    """b = Q^n - dt div(J^n) - div(Lambda dEapp z) + A_Lambda phi^n (cell-integrated)."""
    nr, nz = ne.shape
    b = np.empty(nr * nz)
    for i in range(nr):
        for j in range(nz):
            q = E_CH * (ni[i, j] - ne[i, j]) * V[i, j]
            divJ = E_CH * (Ar[i + 1, j] * (Fri[i + 1, j] - Fre[i + 1, j]) - Ar[i, j] * (Fri[i, j] - Fre[i, j])
                           + Az[i, j + 1] * (Fzi[i, j + 1] - Fze[i, j + 1]) - Az[i, j] * (Fzi[i, j] - Fze[i, j]))
            app = (Az[i, j + 1] * lam_z[i, j + 1] - Az[i, j] * lam_z[i, j]) * dEapp
            p = phi[i, j]
            wr0 = lam_r[i, j] * cr[i, j]
            wr1 = lam_r[i + 1, j] * cr[i + 1, j]
            wz0 = lam_z[i, j] * cz[i, j]
            wz1 = lam_z[i, j + 1] * cz[i, j + 1]
            op = (wr0 + wr1 + wz0 + wz1) * p
            if i > 0:
                op -= wr0 * phi[i - 1, j]
            if i < nr - 1:
                op -= wr1 * phi[i + 1, j]
            if j > 0:
                op -= wz0 * phi[i, j - 1]
            if j < nz - 1:
                op -= wz1 * phi[i, j + 1]
            b[i * nz + j] = q - dt * divJ - app + op
    return b


@njit(cache=True)
def k_update(n, Fr, Fz, V, Ar, Az, src, dt):
    """n + dt[-div(F)/V + src], with src = nu_iz n_e - alpha n_e n_i (same for both species).

    Returns new density, its minimum and the boundary outflow rate (particles/s).
    """
    nr, nz = n.shape
    out = np.empty((nr, nz))
    mn = np.inf
    for i in range(nr):
        for j in range(nz):
            div = (Ar[i + 1, j] * Fr[i + 1, j] - Ar[i, j] * Fr[i, j]
                   + Az[i, j + 1] * Fz[i, j + 1] - Az[i, j] * Fz[i, j])
            v = n[i, j] - dt * div / V[i, j] + dt * src[i, j]
            out[i, j] = v
            mn = min(mn, v)
    flow = 0.0
    for j in range(nz):
        flow += Ar[nr, j] * Fr[nr, j]
    for i in range(nr):
        flow += Az[i, nz] * Fz[i, nz] - Az[i, 0] * Fz[i, 0]
    return out, mn, flow


# --------------------------------------------------------------------------- model
class Streamer2D:
    def __init__(self, cfg: Config, grid: Grid | None = None):
        self.cfg = cfg
        self.g = grid if grid is not None else Grid.from_config(cfg)
        g = self.g
        self.N = cfg.N
        self.swarm = Swarm(cfg.swarm_file, self.N)
        self.mu_i = cfg.K0_ion_m2_Vs * N_LOSCHMIDT / self.N
        self.D_i = self.mu_i * KB * cfg.T_gas_K / E_CH
        R2 = g.rc[:, None] ** 2
        Z2 = g.zc[None, :] ** 2
        self.ne = cfg.n0_m3 * np.exp(-R2 / cfg.a_m ** 2 - Z2 / cfg.L_m ** 2) + cfg.n_bg_m3
        self.ni = self.ne.copy()
        self.t = 0.0
        self.step_count = 0
        self.phi = np.zeros(g.shape)
        # vacuum Laplacian, factorised once per grid
        self.A0 = assemble(g, EPS0 * g.cr, EPS0 * g.cz)
        self.lu0 = factor(self.A0)
        self.lu = None
        self.lam_r = None
        self.lam_z = None
        self.last_refactor_step = -10**9
        self.n_refactor = 0
        self.n_forced_refactor = 0
        self.out_e = 0.0   # cumulative particles lost through boundaries
        self.out_i = 0.0
        self.Q0 = self.total_charge()
        # electron-mobility envelope for the frozen semi-implicit weight (see k_sigma_env)
        mu_cap = float(self.swarm.mu.max()) if cfg.si_mu_envelope == "max" else 0.0
        self._mu_env_value = mu_cap
        self.phi = self.solve_vacuum((self.ni - self.ne) * E_CH * g.V)
        self._update_fields(self.phi, self.t)

    # ---------------- helpers
    def E_app(self, t):
        return self.cfg.E0_V_m * math.sin(2 * math.pi * self.cfg.freq_Hz * t)

    def total_charge(self):
        return float(E_CH * np.sum((self.ni - self.ne) * self.g.V))

    def charge_residual(self):
        """Q(t) + charge carried out - Q(0), in C."""
        return self.total_charge() + E_CH * (self.out_i - self.out_e) - self.Q0

    def solve_vacuum(self, rhs_cells):
        return self.lu0.solve(rhs_cells.ravel()).reshape(self.g.shape)

    def _update_fields(self, phi, t):
        sw = self.swarm
        (self.Er, self.Ez, self.erc, self.ezc, self.Eabs, self.mu_e, self.D_e, nu) = k_fields(
            phi, self.E_app(t), self.g.hr, self.g.hz, 1.0 / (self.N * TD), sw.x0, sw.x1,
            float(sw.xf[1] - sw.xf[0]), sw.mu, sw.D, sw.nu)
        self.nu_iz = nu if self.cfg.ionization else np.zeros_like(nu)

    @property
    def _mu_env(self):
        return np.maximum(self.mu_e, self._mu_env_value)

    def _refactor(self, lam_r, lam_z):
        g = self.g
        self.lam_r, self.lam_z = lam_r, lam_z
        self.lu = factor(assemble(g, (EPS0 + lam_r) * g.cr, (EPS0 + lam_z) * g.cz))
        self.n_refactor += 1
        self.last_refactor_step = self.step_count

    # ---------------- one time step
    def step(self, dt_max=np.inf):
        cfg, g = self.cfg, self.g
        if not hasattr(self, "_mu_i_c"):
            self._mu_i_c = np.full(g.shape, self.mu_i)
            self._D_i_c = np.full(g.shape, self.D_i)
        # fluxes at the current field (for the linearisation and the dt bound)
        Fre, Fze, Gre, Gze, oute = sg_fluxes(self.ne, self.mu_e, self.D_e, -1.0, self.Er, self.Ez,
                                             g.hr, g.hz, g.V, g.Ar, g.Az, True)
        Fri, Fzi, Gri, Gzi, outi = sg_fluxes(self.ni, self._mu_i_c, self._D_i_c, 1.0, self.Er, self.Ez,
                                             g.hr, g.hz, g.V, g.Ar, g.Az, True)
        loss = cfg.alpha_rec_m3_s * max(float(self.ni.max()), float(self.ne.max()))
        rate = max(float(oute.max()), float(outi.max()))
        dt = cfg.safety / rate
        src = float(self.nu_iz.max()) + loss
        if src > 0:
            dt = min(dt, cfg.src_dt_frac / src)
        dt = min(dt, 1.0 / (cfg.freq_Hz * cfg.steps_per_rf_period), dt_max)
        t1 = self.t + dt
        S = self.nu_iz * self.ne
        if cfg.alpha_rec_m3_s > 0:
            S = S - cfg.alpha_rec_m3_s * self.ne * self.ni
        sw = self.swarm
        for attempt in range(6):
            phi1 = self._poisson(dt, t1, Fre, Fze, Fri, Fzi, Gre, Gze, Gri, Gzi)
            new = k_fields(phi1, self.E_app(t1), g.hr, g.hz, 1.0 / (self.N * TD), sw.x0, sw.x1,
                           float(sw.xf[1] - sw.xf[0]), sw.mu, sw.D, sw.nu)
            Er1, Ez1 = new[0], new[1]
            # fluxes: predicted field E^{n+1}, transport coefficients and densities at t^n
            Fre1, Fze1, _, _, _ = sg_fluxes(self.ne, self.mu_e, self.D_e, -1.0, Er1, Ez1,
                                            g.hr, g.hz, g.V, g.Ar, g.Az, False)
            Fri1, Fzi1, _, _, _ = sg_fluxes(self.ni, self._mu_i_c, self._D_i_c, 1.0, Er1, Ez1,
                                            g.hr, g.hz, g.V, g.Ar, g.Az, False)
            ne1, mne, flow_e = k_update(self.ne, Fre1, Fze1, g.V, g.Ar, g.Az, S, dt)
            ni1, mni, flow_i = k_update(self.ni, Fri1, Fzi1, g.V, g.Ar, g.Az, S, dt)
            if mne >= 0.0 and mni >= 0.0:
                break
            # positivity violated with the predicted field: retry with a smaller step (no floor)
            self.n_dt_retries = getattr(self, "n_dt_retries", 0) + 1
            dt *= 0.5
            t1 = self.t + dt
        else:
            raise FloatingPointError(f"negative density at t={self.t:.4e} s after dt reductions")
        self.out_e += dt * flow_e
        self.out_i += dt * flow_i
        self.ne, self.ni, self.phi, self.t = ne1, ni1, phi1, t1
        self.step_count += 1
        (self.Er, self.Ez, self.erc, self.ezc, self.Eabs, self.mu_e, self.D_e, nu) = new
        self.nu_iz = nu if cfg.ionization else np.zeros_like(nu)
        return dt

    def _poisson(self, dt, t1, Fre, Fze, Fri, Fzi, Gre, Gze, Gri, Gzi):
        g, cfg = self.g, self.cfg
        if cfg.si_mode == "explicit":
            return self.solve_vacuum(E_CH * (self.ni - self.ne) * g.V)
        # sigma_f = d(charge flux)/d(face field) of the SG fluxes (>= 0)
        have = self.lam_r is not None
        sr, sz, worst = k_sigma(self.mu_e, self.mu_i, Gre, Gze, Gri, Gzi, dt,
                                self.lam_r if have else sr_dummy(g), self.lam_z if have else sz_dummy(g), have)
        if cfg.si_mode == "exact":
            self._refactor(dt * sr, dt * sz)
        else:
            violated = (not have) or worst > cfg.si_violation
            if violated or self.step_count - self.last_refactor_step >= cfg.si_refresh_steps:
                if violated and have:
                    self.n_forced_refactor += 1
                er, ez = k_sigma_env(self.ne, self.ni, self._mu_env, self.mu_i)
                self._refactor(cfg.si_headroom * dt * np.maximum(er, sr), cfg.si_headroom * dt * np.maximum(ez, sz))
        dEapp = self.E_app(t1) - self.E_app(self.t)
        # Gauss's law for the predicted charge Q^n - dt div J^n - div(Lambda (E^{n+1} - E^n))
        b = k_rhs(self.ne, self.ni, g.V, g.Ar, g.Az, g.cr, g.cz, Fre, Fze, Fri, Fzi,
                  self.lam_r, self.lam_z, self.phi, dt, dEapp)
        return self.lu.solve(b).reshape(g.shape)

    # ---------------- diagnostics
    def sigma(self):
        return E_CH * (self.mu_e * self.ne + self.mu_i * self.ni)

    def diagnostics(self) -> dict:
        g, cfg = self.g, self.cfg
        ne = self.ne
        j0 = int(np.argmin(np.abs(g.zc)))
        d = {"t_s": self.t, "E_app_V_m": self.E_app(self.t), "step": self.step_count}
        d["peak_ne_m3"] = float(ne.max())
        d["peak_E_V_m"] = float(self.Eabs.max())
        ip, jp = np.unravel_index(int(np.argmax(self.Eabs)), g.shape)
        d["peak_E_r_m"], d["peak_E_z_m"] = float(g.rc[ip]), float(g.zc[jp])
        d["r_1e20_mid_m"] = contour_radius(g.rc, ne[:, j0], 1e20)
        d["r_1e20_max_m"] = float(max(contour_radius(g.rc, ne[:, j], 1e20) for j in range(g.nz)))
        s = self.sigma()
        w = s[:, j0] * (g.V[:, j0])
        d["r_sigma_mid_m"] = float(math.sqrt(np.sum(w * g.rc ** 2) / np.sum(w)))
        d["line_conductance_mid_S_m"] = float(np.sum(s[:, j0] * g.V[:, j0] / g.dzc[j0]))
        wall = s * g.V
        d["r_sigma_global_m"] = float(math.sqrt(np.sum(wall * g.rc[:, None] ** 2) / np.sum(wall)))
        colmax = ne.max(axis=0)
        for thr in cfg.thresholds_m3:
            zm, zp = axial_extent(g.zc, colmax, thr)
            key = f"{thr:.0e}".replace("+", "")
            d[f"z_minus_{key}_m"], d[f"z_plus_{key}_m"] = zm, zp
            d[f"length_{key}_m"] = (zp - zm) if np.isfinite(zp) else 0.0
        d["N_e"] = float(np.sum(ne * g.V))
        d["charge_residual_C"] = self.charge_residual()
        d["Q_e_total_C"] = E_CH * d["N_e"]
        return d

    def snapshot(self, path: Path):
        g = self.g
        np.savez_compressed(path, t_s=self.t, rc=g.rc, zc=g.zc, rf=g.rf, zf=g.zf,
                            ne=self.ne.astype(np.float32), ni=self.ni.astype(np.float32),
                            phi=self.phi.astype(np.float32), Er=self.erc.astype(np.float32),
                            Ez=self.ezc.astype(np.float32), Eabs=self.Eabs.astype(np.float32),
                            E_app=self.E_app(self.t))


def contour_radius(rc, col, thr):
    idx = np.nonzero(col >= thr)[0]
    if idx.size == 0:
        return 0.0
    i = idx[-1]
    if i == len(rc) - 1:
        return float(rc[i])
    a, b = math.log(col[i]), math.log(max(col[i + 1], 1e-300))
    f = (a - math.log(thr)) / (a - b) if a != b else 0.0
    return float(rc[i] + f * (rc[i + 1] - rc[i]))


def axial_extent(zc, colmax, thr):
    idx = np.nonzero(colmax >= thr)[0]
    if idx.size == 0:
        return float("nan"), float("nan")

    def interp(i, k):
        if k < 0 or k >= len(zc):
            return float(zc[i])
        a, b = math.log(colmax[i]), math.log(max(colmax[k], 1e-300))
        f = (a - math.log(thr)) / (a - b) if a != b else 0.0
        return float(zc[i] + f * (zc[k] - zc[i]))

    return interp(idx[0], idx[0] - 1), interp(idx[-1], idx[-1] + 1)


# --------------------------------------------------------------------------- driver
def tip_velocity(hist: dict, t_min=1e-9, key="1e20"):
    """Linear-fit speeds of the two n_e >= 1e20 axial tips (m/s) for t >= t_min."""
    t = np.asarray(hist["t_s"])
    zp = np.asarray(hist[f"z_plus_{key}_m"])
    zm = np.asarray(hist[f"z_minus_{key}_m"])
    m = (t >= t_min) & np.isfinite(zp) & np.isfinite(zm)
    if m.sum() < 5:
        return {"v_plus_m_s": float("nan"), "v_minus_m_s": float("nan"), "v_mean_m_s": float("nan")}
    vp = float(np.polyfit(t[m], zp[m], 1)[0])
    vm = float(-np.polyfit(t[m], zm[m], 1)[0])
    return {"v_plus_m_s": vp, "v_minus_m_s": vm, "v_mean_m_s": 0.5 * (vp + vm),
            "fit_window_s": [float(t[m][0]), float(t[m][-1])]}


def run_case(cfg: Config, outdir: Path | None = None, verbose=True) -> dict:
    outdir = Path(outdir) if outdir else RESULTS / cfg.name
    outdir.mkdir(parents=True, exist_ok=True)
    wall0 = time.time()
    model = Streamer2D(cfg)
    g = model.g
    hist: dict = {}

    def record():
        for k, v in model.diagnostics().items():
            hist.setdefault(k, []).append(v)

    snaps = sorted(t for t in cfg.snapshot_times_s if t <= cfg.t_end_s + 1e-15)
    snap_files = []
    record()
    model.snapshot(outdir / "snap_t0000.0ps.npz"); snap_files.append("snap_t0000.0ps.npz")
    snaps = [t for t in snaps if t > 0]
    next_diag = cfg.diag_dt_s
    stop_reason = "t_end"
    dts = []
    if verbose:
        print(f"[{cfg.name}] grid {g.nr}x{g.nz} = {g.nr * g.nz} cells, N={model.N:.3e} m^-3, "
              f"E/N={cfg.EN_Td:.1f} Td, mu_i={model.mu_i:.3e} m^2/Vs", flush=True)
    while model.t < cfg.t_end_s * (1 - 1e-12):
        targets = [next_diag, cfg.t_end_s] + snaps
        dt_max = min(t for t in targets if t > model.t + 1e-18) - model.t
        dts.append(model.step(dt_max=dt_max))
        if model.t >= next_diag * (1 - 1e-9):
            record()
            next_diag += cfg.diag_dt_s
            if verbose and len(hist["t_s"]) % 25 == 0:
                print(f"  t={model.t*1e9:7.3f} ns step={model.step_count} dt={dts[-1]:.2e} "
                      f"ne_max={hist['peak_ne_m3'][-1]:.3e} Emax={hist['peak_E_V_m'][-1]:.3e} "
                      f"r20={hist['r_1e20_max_m'][-1]*1e6:.1f}um "
                      f"z20=[{hist['z_minus_1e20_m'][-1]*1e3:.3f},{hist['z_plus_1e20_m'][-1]*1e3:.3f}]mm "
                      f"refac={model.n_refactor}({model.n_forced_refactor} forced) "
                      f"Qres={hist['charge_residual_C'][-1]/hist['Q_e_total_C'][0]:.1e} "
                      f"wall={time.time()-wall0:.0f}s", flush=True)
            zp, zm = hist["z_plus_1e20_m"][-1], hist["z_minus_1e20_m"][-1]
            rr = hist["r_1e20_max_m"][-1]
            if (np.isfinite(zp) and (zp > g.zf[-1] - cfg.stop_margin_m or zm < g.zf[0] + cfg.stop_margin_m)):
                stop_reason = "channel_reached_z_boundary"
            elif rr > g.rf[-1] - cfg.stop_margin_m:
                stop_reason = "channel_reached_r_boundary"
            elif not np.isfinite(zp):
                stop_reason = "channel_decayed_below_1e20"
            elif time.time() - wall0 > cfg.wall_limit_s:
                stop_reason = "wall_clock_limit"
            if stop_reason != "t_end":
                break
        if snaps and model.t >= snaps[0] * (1 - 1e-9):
            fn = f"snap_t{model.t*1e12:06.1f}ps.npz"
            model.snapshot(outdir / fn); snap_files.append(fn)
            snaps.pop(0)
    if hist["t_s"][-1] < model.t:
        record()
    fn = f"snap_t{model.t*1e12:06.1f}ps.npz"
    if fn not in snap_files:
        model.snapshot(outdir / fn); snap_files.append(fn)
    wall = time.time() - wall0
    summary = build_summary(cfg, model, hist, snap_files, stop_reason, wall, np.array(dts))
    (outdir / "summary.json").write_text(json.dumps(summary, indent=1, default=_json_default), encoding="utf-8")
    np.savez_compressed(outdir / "history.npz", **{k: np.asarray(v, float) for k, v in hist.items()})
    if verbose:
        print(f"[{cfg.name}] done: {stop_reason} at t={model.t*1e9:.3f} ns, {model.step_count} steps, "
              f"{wall:.0f} s", flush=True)
    return summary


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def _at(hist, key, t):
    tt = np.asarray(hist["t_s"])
    if t > tt[-1] + 1e-15:
        return None
    return float(np.interp(t, tt, np.asarray(hist[key], float)))


def build_summary(cfg, model, hist, snap_files, stop_reason, wall, dts):
    g = model.g
    sw = model.swarm
    tv = tip_velocity(hist)
    Q0 = hist["Q_e_total_C"][0]
    res = np.abs(np.asarray(hist["charge_residual_C"])) / Q0
    return {
        "case": cfg.name,
        "config": asdict(cfg),
        "derived": {"N_m3": model.N, "E0_over_N_Td": cfg.EN_Td, "mu_i_m2_Vs": model.mu_i, "D_i_m2_s": model.D_i,
                    "kiz_zero_below_Td": sw.kiz_EN_min_Td,
                    "mu_e_at_E0_m2_Vs": float(sw.at(np.array([cfg.EN_Td]))[0][0]),
                    "nu_iz_at_E0_1_s": float(sw.at(np.array([cfg.EN_Td]))[2][0]),
                    "grid": {"nr": g.nr, "nz": g.nz, "dr_min_m": float(g.rf[1]), "dr_max_m": float(np.diff(g.rf).max()),
                             "dz_m": float(g.dzc[0]), "r_max_m": float(g.rf[-1]), "z_half_m": float(g.zf[-1])}},
        "status": {"stop_reason": stop_reason, "t_final_s": model.t, "steps": model.step_count,
                   "wall_s": wall, "dt_min_s": float(dts.min()) if dts.size else None,
                   "dt_median_s": float(np.median(dts)) if dts.size else None,
                   "lu_refactorisations": model.n_refactor,
                   "forced_lu_refactorisations": model.n_forced_refactor,
                   "max_rel_charge_residual": float(res.max())},
        "results": {
            "radius_1e20_max_m_at": {f"{t*1e9:g}ns": _at(hist, "r_1e20_max_m", t) for t in (0, 1e-9, 2e-9, 5e-9, 10e-9)},
            "radius_1e20_mid_m_at": {f"{t*1e9:g}ns": _at(hist, "r_1e20_mid_m", t) for t in (0, 1e-9, 2e-9, 5e-9, 10e-9)},
            "radius_sigma_mid_m_at": {f"{t*1e9:g}ns": _at(hist, "r_sigma_mid_m", t) for t in (0, 1e-9, 2e-9, 5e-9, 10e-9)},
            "length_1e20_m_at": {f"{t*1e9:g}ns": _at(hist, "length_1e20_m", t) for t in (0, 1e-9, 2e-9, 5e-9, 10e-9)},
            "tip_velocity_1e20": tv,
            "peak_ne_m3": float(np.max(hist["peak_ne_m3"])),
            "peak_E_V_m": float(np.max(hist["peak_E_V_m"])),
            "final_peak_ne_m3": hist["peak_ne_m3"][-1],
        },
        "time_series": {k: hist[k] for k in ("t_s", "E_app_V_m", "r_1e20_max_m", "r_1e20_mid_m", "r_sigma_mid_m",
                                              "r_sigma_global_m", "length_1e20_m", "z_plus_1e20_m", "z_minus_1e20_m",
                                              "length_1e19_m", "length_1e21_m", "peak_ne_m3", "peak_E_V_m",
                                              "line_conductance_mid_S_m", "charge_residual_C")},
        "snapshots": snap_files,
        "environment": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
                        "numba": numba.__version__, "platform": platform.platform()},
        "command": " ".join(sys.argv),
    }


# --------------------------------------------------------------------------- figures
def make_figures(cases=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    figdir = RESULTS / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    names = cases or [n for n in CASES if (RESULTS / n / "summary.json").exists()]
    sums = {n: json.loads((RESULTS / n / "summary.json").read_text()) for n in names}
    made = []
    for n, s in sums.items():
        files = sorted((RESULTS / n).glob("snap_t*.npz"))
        times = [float(np.load(f)["t_s"]) for f in files]
        want = [1e-9, 2e-9, 5e-9, 10e-9]
        pick = []
        for w in want:
            k = int(np.argmin([abs(t - w) for t in times]))
            if files[k] not in pick:
                pick.append(files[k])
        if len(pick) < 4:
            for f in files[::-1]:
                if f not in pick and len(pick) < 4:
                    pick.append(f)
            pick.sort()
        fig, axs = plt.subplots(len(pick), 1, figsize=(9, 2.1 * len(pick)), sharex=True, constrained_layout=True)
        axs = np.atleast_1d(axs)
        for ax, f in zip(axs, pick):
            d = np.load(f)
            rc, zc, ne = d["rc"], d["zc"], d["ne"].astype(float)
            rmax = 0.6e-3
            m = rc <= rmax
            R = np.concatenate([-rc[m][::-1], rc[m]])
            Z = np.concatenate([ne[m][::-1], ne[m]], axis=0)
            pc = ax.pcolormesh(zc * 1e3, R * 1e3, np.maximum(Z, 1e14), norm=LogNorm(1e16, 1e24),
                               cmap="inferno", shading="auto", rasterized=True)
            ax.contour(zc * 1e3, R * 1e3, Z, levels=[1e20], colors="c", linewidths=0.7)
            ax.set_ylabel("r (mm)")
            ax.set_title(f"t = {float(d['t_s'])*1e9:.2f} ns, E_app = {float(d['E_app'])/1e5:.1f} kV/cm", fontsize=9)
        axs[-1].set_xlabel("z (mm)")
        fig.colorbar(pc, ax=axs, label="n_e (m$^{-3}$); cyan: 1e20 contour", shrink=0.8)
        fig.suptitle(f"{n}: E0/N = {s['derived']['E0_over_N_Td']:.1f} Td (simulation result)", fontsize=10)
        p = figdir / f"ne_map_{n}.png"
        fig.savefig(p, dpi=130); plt.close(fig); made.append(p)

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    for n, s in sums.items():
        ts = s["time_series"]
        t = np.array(ts["t_s"]) * 1e9
        l, = axs[0].plot(t, np.array(ts["r_1e20_max_m"]) * 1e6, label=n, lw=1)
        axs[0].plot(t, np.array(ts["r_sigma_mid_m"]) * 1e6, ls=":", color=l.get_color(), lw=1)
        axs[1].plot(t, np.array(ts["length_1e20_m"]) * 1e3, color=l.get_color(), lw=1)
        axs[2].semilogy(t, ts["peak_ne_m3"], color=l.get_color(), lw=1)
    axs[0].set(xlabel="t (ns)", ylabel="radius (µm)", title="solid: n_e≥1e20 half-width (max over z); dotted: σ-weighted (z=0)")
    axs[1].set(xlabel="t (ns)", ylabel="axial extent n_e≥1e20 (mm)", title="channel length (domain ±4 mm)")
    axs[2].set(xlabel="t (ns)", ylabel="peak n_e (m$^{-3}$)", title="peak electron density")
    axs[0].title.set_fontsize(8)
    axs[0].legend(fontsize=7)
    p = figdir / "radius_length_vs_time.png"
    fig.savefig(p, dpi=130); plt.close(fig); made.append(p)

    fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
    for n, s in sums.items():
        v = s["results"]["tip_velocity_1e20"]
        if v.get("v_mean_m_s") is not None and np.isfinite(v["v_mean_m_s"]):
            ax.errorbar(s["derived"]["E0_over_N_Td"], v["v_mean_m_s"],
                        yerr=[[max(v["v_mean_m_s"] - min(v["v_plus_m_s"], v["v_minus_m_s"]), 0)],
                              [max(max(v["v_plus_m_s"], v["v_minus_m_s"]) - v["v_mean_m_s"], 0)]],
                        fmt="o", capsize=3, label=f"{n} ({s['config']['pressure_Pa']/1e6:g} MPa)")
    ax.axhline(5e5, color="k", ls="--", label="measured AIR filament, ~5e7 cm/s (Technion)")
    ax.set(xlabel="applied peak E0/N (Td)", ylabel="tip speed of n_e≥1e20 contour (m/s)",
           title="He model tip speed vs E/N (air value is NOT helium)")
    ax.set_yscale("symlog", linthresh=1e3)
    ax.legend(fontsize=7)
    p = figdir / "tip_velocity_vs_EN.png"
    fig.savefig(p, dpi=130); plt.close(fig); made.append(p)

    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    labels = []
    for k, (n, s) in enumerate(sums.items()):
        ts = s["time_series"]
        d20 = 2 * np.array(ts["r_1e20_max_m"]) * 1e3
        ds = 2 * np.array(ts["r_sigma_mid_m"]) * 1e3
        ax.plot([k - 0.15, k - 0.15], [d20.min(), d20.max()], color="C0", lw=6, alpha=0.5)
        ax.plot(k - 0.15, d20[-1], "o", color="C0")
        ax.plot([k + 0.15, k + 0.15], [ds.min(), ds.max()], color="C1", lw=6, alpha=0.5)
        ax.plot(k + 0.15, ds[-1], "s", color="C1")
        labels.append(n)
    ax.axhspan(0.2, 0.6, color="0.8", label="measured AIR filament diameter 0.2-0.6 mm (Technion)")
    ax.plot([], [], "o", color="C0", label="2 x r(n_e≥1e20), final (bar: range)")
    ax.plot([], [], "s", color="C1", label="2 x σ-weighted radius at z=0, final (bar: range)")
    ax.set_xticks(range(len(labels)), labels, rotation=30, ha="right", fontsize=7)
    ax.set_yscale("log")
    ax.set_ylabel("channel diameter (mm)")
    ax.set_title("Model channel diameter vs measured air value (not helium)")
    ax.legend(fontsize=7)
    p = figdir / "diameter_vs_measured.png"
    fig.savefig(p, dpi=130); plt.close(fig); made.append(p)
    return made


# --------------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser(prog="streamer2d", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    r = sub.add_parser("run")
    r.add_argument("cases", nargs="*")
    r.add_argument("--all", action="store_true")
    r.add_argument("--t-end", type=float, help="end time (s)")
    r.add_argument("--dz", type=float)
    r.add_argument("--dr-min", type=float)
    r.add_argument("--suffix", default="", help="appended to the case/output name")
    r.add_argument("--alpha", type=float, help="override recombination coefficient (m^3/s)")
    r.add_argument("--z-half", type=float, help="axial half-length of the domain (m), default 4e-3")
    f = sub.add_parser("figures")
    f.add_argument("cases", nargs="*")
    a = ap.parse_args(argv)
    if a.cmd == "list":
        for n, c in CASES.items():
            print(f"{n:28s} p={c.pressure_Pa/1e6:6.3f} MPa  E0={c.E0_V_m/1e5:7.2f} kV/cm  E0/N={c.EN_Td:6.1f} Td  "
                  f"n0={c.n0_m3:.0e} a={c.a_m*1e6:g} um bg={c.n_bg_m3:.0e}")
        return 0
    if a.cmd == "figures":
        for p in make_figures(a.cases or None):
            print(p)
        return 0
    names = list(CASES) if a.all else a.cases
    for n in names:
        cfg = CASES[n]
        kw = {}
        if a.t_end:
            kw["t_end_s"] = a.t_end
        if a.dz:
            kw["dz_m"] = a.dz
        if a.dr_min:
            kw["dr_min_m"] = a.dr_min
        if a.alpha is not None:
            kw["alpha_rec_m3_s"] = a.alpha
        if a.z_half:
            kw["z_half_m"] = a.z_half
        cfg = replace(cfg, name=n + a.suffix, **kw)
        run_case(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
