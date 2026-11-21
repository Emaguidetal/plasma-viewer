"""Uniform-grid 3D Yee FDTD with PEC masks, z-CPML and a time-varying Drude current.

Field layout (cells nx, ny, nz; spacing dx):
  Ex (nx, ny+1, nz+1) at ((i+1/2)dx, j dx, k dx)
  Ey (nx+1, ny, nz+1) at (i dx, (j+1/2)dx, k dx)
  Ez (nx+1, ny+1, nz) at (i dx, j dx, (k+1/2)dx)
  Hx (nx+1, ny, nz), Hy (nx, ny+1, nz), Hz (nx, ny, nz+1)

The electron current obeys dJ/dt + nu J = eps0 wp^2 E and is advanced with the
trapezoidal (Crank-Nicolson) auxiliary-differential-equation scheme. The local
plasma update is unconditionally stable for any density or collision rate, and
the work J.E is accumulated exactly as the discrete energy removed from the field.
"""

from dataclasses import dataclass, field
import math

import numpy as np
from numba import njit, prange
from scipy.constants import epsilon_0 as EPS0, mu_0 as MU0, c as C0


@njit(parallel=True, cache=True)
def _update_h(Ex, Ey, Ez, Hx, Hy, Hz, ch, psi_hx, psi_hy, bh, chz):
    nxp1, ny, nz = Hx.shape
    for i in prange(nxp1):
        for j in range(ny):
            for k in range(nz):
                dz = Ey[i, j, k + 1] - Ey[i, j, k]
                if bh[k] != 0.0:
                    psi_hx[i, j, k] = bh[k] * psi_hx[i, j, k] + chz[k] * dz
                    dz += psi_hx[i, j, k]
                Hx[i, j, k] -= ch * ((Ez[i, j + 1, k] - Ez[i, j, k]) - dz)
    nx, nyp1, nz = Hy.shape
    for i in prange(nx):
        for j in range(nyp1):
            for k in range(nz):
                dz = Ex[i, j, k + 1] - Ex[i, j, k]
                if bh[k] != 0.0:
                    psi_hy[i, j, k] = bh[k] * psi_hy[i, j, k] + chz[k] * dz
                    dz += psi_hy[i, j, k]
                Hy[i, j, k] -= ch * (dz - (Ez[i + 1, j, k] - Ez[i, j, k]))
    nx, ny, nzp1 = Hz.shape
    for i in prange(nx):
        for j in range(ny):
            for k in range(nzp1):
                Hz[i, j, k] -= ch * ((Ey[i + 1, j, k] - Ey[i, j, k]) - (Ex[i, j + 1, k] - Ex[i, j, k]))


@njit(cache=True, inline="always")
def _drude_step(e_old, curl, ce, a, g, j_old, half_dt_eps):
    # E^{n+1}(1+g) = E^n(1-g) + ce curl - dt/(2 eps0)(1+a) J^n ;  J^{n+1} = a J^n + b (E^{n+1}+E^n)/2
    e_new = (e_old * (1.0 - g) + ce * curl - half_dt_eps * (1.0 + a) * j_old) / (1.0 + g)
    return e_new


@njit(parallel=True, cache=True)
def _update_e(Ex, Ey, Ez, Hx, Hy, Hz, ce, dt, mEx, mEy, mEz,
              aEx, aEy, aEz, gEx, gEy, gEz, Jx, Jy, Jz, Px, Py, Pz,
              psi_ex, psi_ey, be, cez):
    half_dt_eps = dt / (2.0 * EPS0)
    b_scale = 4.0 * EPS0 / dt
    nx, nyp1, nzp1 = Ex.shape
    for i in prange(nx):
        for j in range(1, nyp1 - 1):
            for k in range(1, nzp1 - 1):
                if mEx[i, j, k] == 0:
                    continue
                dz = Hy[i, j, k] - Hy[i, j, k - 1]
                if be[k] != 0.0:
                    psi_ex[i, j, k] = be[k] * psi_ex[i, j, k] + cez[k] * dz
                    dz += psi_ex[i, j, k]
                curl = (Hz[i, j, k] - Hz[i, j - 1, k]) - dz
                eo = Ex[i, j, k]
                g = gEx[i, j, k]
                if g == 0.0 and Jx[i, j, k] == 0.0:
                    Ex[i, j, k] = eo + ce * curl
                else:
                    a = aEx[i, j, k]
                    en = _drude_step(eo, curl, ce, a, g, Jx[i, j, k], half_dt_eps)
                    jo = Jx[i, j, k]
                    jn = a * jo + b_scale * g * 0.5 * (en + eo)
                    Px[i, j, k] += 0.25 * (jn + jo) * (en + eo) * dt
                    Jx[i, j, k] = jn
                    Ex[i, j, k] = en
    nxp1, ny, nzp1 = Ey.shape
    for i in prange(1, nxp1 - 1):
        for j in range(ny):
            for k in range(1, nzp1 - 1):
                if mEy[i, j, k] == 0:
                    continue
                dz = Hx[i, j, k] - Hx[i, j, k - 1]
                if be[k] != 0.0:
                    psi_ey[i, j, k] = be[k] * psi_ey[i, j, k] + cez[k] * dz
                    dz += psi_ey[i, j, k]
                curl = dz - (Hz[i, j, k] - Hz[i - 1, j, k])
                eo = Ey[i, j, k]
                g = gEy[i, j, k]
                if g == 0.0 and Jy[i, j, k] == 0.0:
                    Ey[i, j, k] = eo + ce * curl
                else:
                    a = aEy[i, j, k]
                    en = _drude_step(eo, curl, ce, a, g, Jy[i, j, k], half_dt_eps)
                    jo = Jy[i, j, k]
                    jn = a * jo + b_scale * g * 0.5 * (en + eo)
                    Py[i, j, k] += 0.25 * (jn + jo) * (en + eo) * dt
                    Jy[i, j, k] = jn
                    Ey[i, j, k] = en
    nxp1, nyp1, nz = Ez.shape
    for i in prange(1, nxp1 - 1):
        for j in range(1, nyp1 - 1):
            for k in range(nz):
                if mEz[i, j, k] == 0:
                    continue
                curl = (Hy[i, j, k] - Hy[i - 1, j, k]) - (Hx[i, j, k] - Hx[i, j - 1, k])
                eo = Ez[i, j, k]
                g = gEz[i, j, k]
                if g == 0.0 and Jz[i, j, k] == 0.0:
                    Ez[i, j, k] = eo + ce * curl
                else:
                    a = aEz[i, j, k]
                    en = _drude_step(eo, curl, ce, a, g, Jz[i, j, k], half_dt_eps)
                    jo = Jz[i, j, k]
                    jn = a * jo + b_scale * g * 0.5 * (en + eo)
                    Pz[i, j, k] += 0.25 * (jn + jo) * (en + eo) * dt
                    Jz[i, j, k] = jn
                    Ez[i, j, k] = en


@njit(cache=True)
def _multipole_edges(Ey, e_before, i0, k0, js, ce_curl, dt, w, a, b, J, P):
    """Replace a vacuum update on selected Ey edges by a multi-ring Drude update.

    ce_curl holds (E_vacuum_new - E_old) i.e. ce*curl for each listed edge.
    w[e, r] is ring area / dx^2; J[e, r] is the physical current density in ring r.
    """
    for e in range(js.shape[0]):
        j = js[e]
        eo = e_before[e]
        gsum = 0.0
        jsum = 0.0
        for r in range(w.shape[1]):
            gsum += dt * w[e, r] * b[e, r] / (4.0 * EPS0)
            jsum += w[e, r] * (1.0 + a[e, r]) * J[e, r]
        en = (eo * (1.0 - gsum) + ce_curl[e] - dt / (2.0 * EPS0) * jsum) / (1.0 + gsum)
        em = 0.5 * (en + eo)
        for r in range(w.shape[1]):
            jo = J[e, r]
            jn = a[e, r] * jo + b[e, r] * em
            P[e, r] += 0.5 * (jn + jo) * em * dt
            J[e, r] = jn
        Ey[i0, j, k0] = en


@njit(parallel=True, cache=True)
def _edge_coefficients(wp2_node, nu_node, dt, aEx, aEy, aEz, gEx, gEy, gEz, mEx, mEy, mEz):
    nxp1, nyp1, nzp1 = wp2_node.shape
    for i in prange(nxp1 - 1):
        for j in range(nyp1):
            for k in range(nzp1):
                if mEx[i, j, k]:
                    w1 = wp2_node[i, j, k]; w2 = wp2_node[i + 1, j, k]
                    wp2 = 0.5 * (w1 + w2)
                    nu = (w1 * nu_node[i, j, k] + w2 * nu_node[i + 1, j, k]) / (w1 + w2) if wp2 > 0 else 0.0
                    d = 1.0 + 0.5 * nu * dt
                    aEx[i, j, k] = (1.0 - 0.5 * nu * dt) / d
                    gEx[i, j, k] = dt * dt * wp2 / (4.0 * d)
    for i in prange(nxp1):
        for j in range(nyp1 - 1):
            for k in range(nzp1):
                if mEy[i, j, k]:
                    w1 = wp2_node[i, j, k]; w2 = wp2_node[i, j + 1, k]
                    wp2 = 0.5 * (w1 + w2)
                    nu = (w1 * nu_node[i, j, k] + w2 * nu_node[i, j + 1, k]) / (w1 + w2) if wp2 > 0 else 0.0
                    d = 1.0 + 0.5 * nu * dt
                    aEy[i, j, k] = (1.0 - 0.5 * nu * dt) / d
                    gEy[i, j, k] = dt * dt * wp2 / (4.0 * d)
    for i in prange(nxp1):
        for j in range(nyp1):
            for k in range(nzp1 - 1):
                if mEz[i, j, k]:
                    w1 = wp2_node[i, j, k]; w2 = wp2_node[i, j, k + 1]
                    wp2 = 0.5 * (w1 + w2)
                    nu = (w1 * nu_node[i, j, k] + w2 * nu_node[i, j, k + 1]) / (w1 + w2) if wp2 > 0 else 0.0
                    d = 1.0 + 0.5 * nu * dt
                    aEz[i, j, k] = (1.0 - 0.5 * nu * dt) / d
                    gEz[i, j, k] = dt * dt * wp2 / (4.0 * d)


@njit(parallel=True, cache=True)
def _edge_power_to_nodes(Px, Py, Pz, node_energy):
    """Distribute accumulated edge work (J/m^3) half to each end node, then clear edges."""
    nxp1, nyp1, nzp1 = node_energy.shape
    for i in prange(nxp1):
        for j in range(nyp1):
            for k in range(nzp1):
                s = 0.0
                if i > 0:
                    s += Px[i - 1, j, k]
                if i < nxp1 - 1:
                    s += Px[i, j, k]
                if j > 0:
                    s += Py[i, j - 1, k]
                if j < nyp1 - 1:
                    s += Py[i, j, k]
                if k > 0:
                    s += Pz[i, j, k - 1]
                if k < nzp1 - 1:
                    s += Pz[i, j, k]
                # each edge has two nodes; a node owns half of six edges -> volume share dx^3
                node_energy[i, j, k] = 0.5 * s
    Px[:] = 0.0
    Py[:] = 0.0
    Pz[:] = 0.0


@njit(parallel=True, cache=True)
def _em_energy(Ex, Ey, Ez, Hx, Hy, Hz, k_lo, k_hi):
    """Electromagnetic energy density sum (J/m^3 * cells) for k_lo <= z-index < k_hi."""
    tot = 0.0
    nx = Ex.shape[0]
    for i in prange(nx + 1):
        s = 0.0
        for k in range(k_lo, k_hi):
            if i < nx:
                for j in range(Ex.shape[1]):
                    s += 0.5 * EPS0 * Ex[i, j, k] ** 2
                for j in range(Hy.shape[1]):
                    s += 0.5 * MU0 * Hy[i, j, k] ** 2
                for j in range(Hz.shape[1]):
                    s += 0.5 * MU0 * Hz[i, j, k] ** 2
            for j in range(Ey.shape[1]):
                s += 0.5 * EPS0 * Ey[i, j, k] ** 2
            for j in range(Ez.shape[1]):
                s += 0.5 * EPS0 * Ez[i, j, k] ** 2
            for j in range(Hx.shape[1]):
                s += 0.5 * MU0 * Hx[i, j, k] ** 2
        tot += s
    return tot


@njit(cache=True)
def _plane_projection(Ey, Hx, k, profile_ey, profile_hx):
    v = 0.0
    cur = 0.0
    for i in range(profile_ey.shape[0]):
        for j in range(profile_ey.shape[1]):
            v += Ey[i, j, k] * profile_ey[i, j]
    for i in range(profile_hx.shape[0]):
        for j in range(profile_hx.shape[1]):
            cur += 0.5 * (Hx[i, j, k - 1] + Hx[i, j, k]) * profile_hx[i, j]
    return v, cur


@njit(cache=True)
def _poynting_z(Ex, Ey, Hx, Hy, k):
    """Instantaneous z Poynting flux sum (W/m^2 * cells) through node plane k, cell-centred."""
    s = 0.0
    nx = Ex.shape[0]
    ny = Ey.shape[1]
    for i in range(nx):
        for j in range(ny):
            ex = 0.5 * (Ex[i, j, k] + Ex[i, j + 1, k])
            ey = 0.5 * (Ey[i, j, k] + Ey[i + 1, j, k])
            hy = 0.25 * (Hy[i, j, k - 1] + Hy[i, j + 1, k - 1] + Hy[i, j, k] + Hy[i, j + 1, k])
            hx = 0.25 * (Hx[i, j, k - 1] + Hx[i + 1, j, k - 1] + Hx[i, j, k] + Hx[i + 1, j, k])
            s += ex * hy - ey * hx
    return s


def cpml_profiles(nz: int, n_lo: int, n_hi: int, dx: float, dt: float, f_ref: float,
                  order: int = 3, reflection: float = 1e-6):
    """CPML recursive-convolution coefficients along z for E (node) and H (half) positions."""
    eta0 = math.sqrt(MU0 / EPS0)
    def coeffs(positions):
        b = np.zeros(len(positions))
        cc = np.zeros(len(positions))
        for idx, zpos in enumerate(positions):
            depth = 0.0
            if n_lo > 0 and zpos < n_lo:
                depth = (n_lo - zpos) / n_lo
            if n_hi > 0 and zpos > nz - n_hi:
                depth = (zpos - (nz - n_hi)) / n_hi
            if depth <= 0:
                continue
            thickness = (n_lo if zpos < n_lo else n_hi) * dx
            sig_max = -(order + 1) * math.log(reflection) / (2 * eta0 * thickness)
            sig = sig_max * depth**order
            alpha = 2 * math.pi * f_ref * EPS0 * 0.1 * (1 - depth)
            bb = math.exp(-(sig + alpha) * dt / EPS0)
            b[idx] = bb
            cc[idx] = sig / (sig + alpha) * (bb - 1) if sig + alpha > 0 else 0.0
        return b, cc
    be, ce = coeffs(np.arange(nz + 1, dtype=float))
    bh, ch = coeffs(np.arange(nz, dtype=float) + 0.5)
    return be, ce, bh, ch


@dataclass
class Grid:
    nx: int
    ny: int
    nz: int
    dx: float
    courant: float = 0.98
    dt: float = field(init=False)

    def __post_init__(self):
        for name in ("nx", "ny", "nz"):
            if getattr(self, name) < 2:
                raise ValueError(f"{name} must be at least 2")
        if not self.dx > 0:
            raise ValueError("dx must be positive")
        self.dt = self.courant * self.dx / (C0 * math.sqrt(3.0))

    @property
    def cells(self) -> int:
        return self.nx * self.ny * self.nz


class FDTD:
    """Field state and stepping. Geometry is given as a node-inside boolean array."""

    def __init__(self, grid: Grid, node_inside: np.ndarray, pml_lo: int = 0, pml_hi: int = 0,
                 f_ref: float = 3e9):
        g = grid
        self.grid = g
        nx, ny, nz = g.nx, g.ny, g.nz
        if node_inside.shape != (nx + 1, ny + 1, nz + 1):
            raise ValueError("node_inside must have shape (nx+1, ny+1, nz+1)")
        f = np.float64
        self.Ex = np.zeros((nx, ny + 1, nz + 1), f)
        self.Ey = np.zeros((nx + 1, ny, nz + 1), f)
        self.Ez = np.zeros((nx + 1, ny + 1, nz), f)
        self.Hx = np.zeros((nx + 1, ny, nz), f)
        self.Hy = np.zeros((nx, ny + 1, nz), f)
        self.Hz = np.zeros((nx, ny, nz + 1), f)
        ins = node_inside.astype(bool)
        self.node_inside = ins
        self.mEx = (ins[:-1, :, :] & ins[1:, :, :]).astype(np.uint8)
        self.mEy = (ins[:, :-1, :] & ins[:, 1:, :]).astype(np.uint8)
        self.mEz = (ins[:, :, :-1] & ins[:, :, 1:]).astype(np.uint8)
        # outer box is PEC for tangential E
        self.mEx[:, 0, :] = 0; self.mEx[:, -1, :] = 0; self.mEx[:, :, 0] = 0; self.mEx[:, :, -1] = 0
        self.mEy[0, :, :] = 0; self.mEy[-1, :, :] = 0; self.mEy[:, :, 0] = 0; self.mEy[:, :, -1] = 0
        self.mEz[0, :, :] = 0; self.mEz[-1, :, :] = 0; self.mEz[:, 0, :] = 0; self.mEz[:, -1, :] = 0
        self.aEx = np.ones_like(self.Ex); self.aEy = np.ones_like(self.Ey); self.aEz = np.ones_like(self.Ez)
        self.gEx = np.zeros_like(self.Ex); self.gEy = np.zeros_like(self.Ey); self.gEz = np.zeros_like(self.Ez)
        self.Jx = np.zeros_like(self.Ex); self.Jy = np.zeros_like(self.Ey); self.Jz = np.zeros_like(self.Ez)
        self.Px = np.zeros_like(self.Ex); self.Py = np.zeros_like(self.Ey); self.Pz = np.zeros_like(self.Ez)
        self.pml_lo, self.pml_hi = pml_lo, pml_hi
        self.be, self.cez, self.bh, self.chz = cpml_profiles(nz, pml_lo, pml_hi, g.dx, g.dt, f_ref)
        self.psi_ex = np.zeros_like(self.Ex); self.psi_ey = np.zeros_like(self.Ey)
        self.psi_hx = np.zeros_like(self.Hx); self.psi_hy = np.zeros_like(self.Hy)
        self.ce = g.dt / (EPS0 * g.dx)
        self.ch = g.dt / (MU0 * g.dx)
        self.step_count = 0
        self.multipole = None  # dict set by attach_multipole

    def attach_multipole(self, i0: int, k0: int, js: np.ndarray, rings: int):
        js = np.asarray(js, dtype=np.int64)
        if np.any(self.mEy[i0, js, k0] == 0):
            raise ValueError("multipole edges must be interior Ey edges")
        n = len(js)
        self.multipole = dict(
            i0=int(i0), k0=int(k0), js=js,
            w=np.zeros((n, rings)), a=np.ones((n, rings)), b=np.zeros((n, rings)),
            J=np.zeros((n, rings)), P=np.zeros((n, rings)),
        )
        # these edges are excluded from the bulk Drude update
        self.gEy[i0, js, k0] = 0.0
        self.Jy[i0, js, k0] = 0.0

    def set_bulk_plasma(self, wp2_node: np.ndarray, nu_node: np.ndarray):
        _edge_coefficients(wp2_node, nu_node, self.grid.dt, self.aEx, self.aEy, self.aEz,
                           self.gEx, self.gEy, self.gEz, self.mEx, self.mEy, self.mEz)
        if self.multipole is not None:
            mp = self.multipole
            self.gEy[mp["i0"], mp["js"], mp["k0"]] = 0.0

    def step(self, source=None):
        """Advance one time step. source(fdtd) may add soft currents to E after the update."""
        _update_h(self.Ex, self.Ey, self.Ez, self.Hx, self.Hy, self.Hz, self.ch,
                  self.psi_hx, self.psi_hy, self.bh, self.chz)
        mp = self.multipole
        if mp is not None:
            e_before = self.Ey[mp["i0"], mp["js"], mp["k0"]].copy()
        _update_e(self.Ex, self.Ey, self.Ez, self.Hx, self.Hy, self.Hz, self.ce, self.grid.dt,
                  self.mEx, self.mEy, self.mEz, self.aEx, self.aEy, self.aEz,
                  self.gEx, self.gEy, self.gEz, self.Jx, self.Jy, self.Jz,
                  self.Px, self.Py, self.Pz, self.psi_ex, self.psi_ey, self.be, self.cez)
        if mp is not None:
            ce_curl = self.Ey[mp["i0"], mp["js"], mp["k0"]] - e_before
            _multipole_edges(self.Ey, e_before, mp["i0"], mp["k0"], mp["js"], ce_curl, self.grid.dt,
                             mp["w"], mp["a"], mp["b"], mp["J"], mp["P"])
        if source is not None:
            source(self)
        self.step_count += 1

    def energy(self, k_lo: int = 0, k_hi: int | None = None) -> float:
        k_hi = self.grid.nz if k_hi is None else k_hi
        return _em_energy(self.Ex, self.Ey, self.Ez, self.Hx, self.Hy, self.Hz, k_lo, k_hi) * self.grid.dx**3

    def collect_bulk_work(self, node_energy: np.ndarray):
        _edge_power_to_nodes(self.Px, self.Py, self.Pz, node_energy)

    def state(self) -> dict:
        names = ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "Jx", "Jy", "Jz", "psi_ex", "psi_ey", "psi_hx", "psi_hy"]
        return {n: getattr(self, n).copy() for n in names} | {"step_count": self.step_count}

    def load_state(self, st: dict):
        for n, v in st.items():
            if n == "step_count":
                self.step_count = int(v)
            else:
                getattr(self, n)[...] = v


def te10_discrete_beta(frequency_hz: float, width_m: float, dx: float, dt: float) -> float:
    """Propagation constant of TE10 on the Yee grid (numerical dispersion relation)."""
    s_t = math.sin(math.pi * frequency_hz * dt) / (C0 * dt)
    s_x = math.sin(math.pi / width_m * dx / 2) / dx
    rem = s_t**2 - s_x**2
    if rem <= 0:
        raise ValueError("TE10 is below cutoff on this grid")
    return 2 / dx * math.asin(dx * math.sqrt(rem))
