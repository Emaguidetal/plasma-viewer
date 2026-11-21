"""3D rod growth along E in a WR-284 cavity gap: a conducting-channel (tree) model.

Geometry: x across the broad wall (0..a), y across the gap along E (0..b), z along the guide (0..d).
The laser crosses the cavity along x at (y, z) = (b/2, d/2), perpendicular to E; the TE101 field at
the laser line is E_y = E0(t) sin(pi x / a) sin(omega t + phase).

The discharge is a graph of thin conducting line elements (nodes) joined by ohmic edges, following
the tree models of Luque & Ebert (New J. Phys. 16, 013039, 2014):
* laser-line nodes along x carry the laser seed and avalanche in the local field before inception;
* when the local electron density on the line reaches n_inc, a vertical rod (3 nodes, spacing h)
  is born there; each rod tip then advances along +-y with a local tip law;
* node charges obey dQ/dt = -sum(I), I_e = G_e (phi_i - phi_j + emf_e), solved implicitly
  with phi = P Q; P is the thin-wire potential matrix including image charges of all six cavity
  walls (Dirichlet box), so rods feel the broad walls, the side walls and each other in 3D;
* the edge line conductance is G = e mu_e(E/N) Lambda / l, and the electron line density Lambda
  of every node grows with nu_iz(E/N) of the local field along the channel (LFA, BOLOS He table)
  and decays by recombination;
* tip law: Delta U = space-charge potential of the tip node (tip potential above the unperturbed
  potential), r_tip = max(r_min, |Delta U| / E_m), E_tip = |Delta U| / r_tip,
  v = c_v r_tip nu_iz(E_tip / N); c_v and E_m are calibrated against resolved axisymmetric
  streamer2d runs (see calibrate_tip_law); new nodes get n_ch = (eps0/e) int_0^E_tip alpha dE;
* the cavity mode amplitude obeys dE0/dt = (E_set - E0)/tau_amp - P_ext / (dW/dE0) with
  W = eps0 E0^2 a b d / 8 and P_ext = sum(I_e emf_e), so the growing rods load the cavity.

Quasi-static: valid while the rods are short compared with lambda/4 and before the bridged rod
current becomes inductance limited; the run therefore stops shortly after the first bridge and the
bridged rod is continued by rod_drive.py.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
from numba import njit, prange

from .breakdown_probability import Setup, seed_per_length
from .streamer2d import SWARM_FILE, Swarm

EPS0 = 8.8541878128e-12
QE = 1.602176634e-19
KB = 1.380649e-23
TD = 1e-21
PSI = 6894.757
ATM = 101325.0
A_WG, B_WG, D_CAV = 0.072136, 0.034036, 0.06927
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "rod3d"
K_COUL = 1.0 / (4.0 * math.pi * EPS0)


@dataclass
class RodConfig:
    name: str = "test"
    p_pa: float = 200 * PSI + ATM
    t_k: float = 300.0
    f_hz: float = 2.95e9
    e0_v_m: float = 2.0e7          # empty-cavity peak field at x = a/2 (drive set point)
    rf_phase: float = 0.0
    x_focus_m: float = 0.0          # laser focus along x relative to the field maximum a/2
    energy_j: float = 0.75e-3
    stochastic: bool = True         # Poisson number of seed electrons per node and Furry gain
    rng_seed: int = 1
    dx_c: float = 0.25e-3           # spacing of laser-line nodes (candidate rod sites)
    x_half: float = 5e-3            # laser-line half span around the focus
    h: float = 0.5e-3               # rod element length
    r_s: float = 6e-6               # seed line radius
    n_inc: float = 1e20             # line density at which a vertical rod is born (m^-3)
    e_m_td: float = 550.0           # tip-field closure E_m / N
    c_v: float = 1.0e-2             # tip-law constant (calibrated)
    r_min: float = 10e-6
    r_max: float = 0.3e-3           # cap on the radius of new channel elements (streamer2d: 0.05-0.3 mm)
    alpha_rec: float = 1e-14        # unverified, as in streamer2d
    n_sat: float = 1.0e24           # channel density cap (streamer2d peak 1e23-7e23 m^-3 at 1.48 MPa)
    max_ion_step: float = 0.5       # max nu*dt per step (adaptive dt)
    tau_amp_s: float = 35.7e-9      # 2 Q_L / omega, Q_L = 331
    cavity_feedback: bool = True
    drive_off_s: float = math.inf   # klystron drive switched off (cavity then rings down with tau_amp)
    walls: bool = True
    n_img: tuple = (1, 3, 1)        # image orders in x, y, z
    t_end_s: float = 20e-9
    steps_per_period: int = 64
    frame_dt_s: float = 50e-12
    stop_after_bridge_s: float = 1.0e-9

    @property
    def n_gas(self):
        return self.p_pa / (KB * self.t_k)


# ----------------------------------------------------------------------------- kernels
@njit(cache=True)
def _seg(px, py, pz, cx, cy, cz, axis, hl, rmin):
    """Potential and E_y at p of a uniform unit line density (1 C/m) on a segment along `axis`."""
    if axis == 0:
        s = px - cx; d1 = py - cy; d2 = pz - cz
    elif axis == 1:
        s = py - cy; d1 = px - cx; d2 = pz - cz
    else:
        s = pz - cz; d1 = px - cx; d2 = py - cy
    rho_t = math.sqrt(d1 * d1 + d2 * d2)
    rho = max(rho_t, rmin)
    l1 = -hl - s
    l2 = hl - s
    phi = math.asinh(l2 / rho) - math.asinh(l1 / rho)
    r1 = math.sqrt(rho * rho + l1 * l1)
    r2 = math.sqrt(rho * rho + l2 * l2)
    if axis == 1:
        ey = 1.0 / r2 - 1.0 / r1
    else:
        ey = 0.0
        if rho_t > 1e-15:
            e_perp = (l2 / r2 - l1 / r1) / rho
            dy = (py - cy) if axis == 0 else 0.0
            if axis == 2:
                dy = py - cy
            ey = e_perp * dy / rho_t
    return phi, ey


@njit(cache=True)
def _seg_images(px, py, pz, cx, cy, cz, axis, hl, rmin, walls, nx, ny, nz, a, b, d):
    if not walls:
        return _seg(px, py, pz, cx, cy, cz, axis, hl, rmin)
    phi = 0.0
    ey = 0.0
    for ix in range(-nx, nx + 1):
        for sx in (1, -1):
            xi = 2 * ix * a + sx * cx
            for iy in range(-ny, ny + 1):
                for sy in (1, -1):
                    yi = 2 * iy * b + sy * cy
                    for iz in range(-nz, nz + 1):
                        for sz in (1, -1):
                            zi = 2 * iz * d + sz * cz
                            sg = sx * sy * sz
                            # thin-wire reduced kernel for every term (collinear y-images of rods)
                            p, e = _seg(px, py, pz, xi, yi, zi, axis, hl, rmin)
                            phi += sg * p
                            ey += sg * e
    return phi, ey


@njit(parallel=True, cache=True)
def _block(rows_pos, rows_r, cols_pos, cols_axis, cols_hl, walls, nx, ny, nz, a, b, d, out_p, out_e):
    """out_p[i, j] = potential at row point i per unit charge on column element j; out_e = E_y."""
    for i in prange(rows_pos.shape[0]):
        for j in range(cols_pos.shape[0]):
            p, e = _seg_images(rows_pos[i, 0], rows_pos[i, 1], rows_pos[i, 2], cols_pos[j, 0], cols_pos[j, 1],
                               cols_pos[j, 2], cols_axis[j], cols_hl[j], rows_r[i], walls, nx, ny, nz, a, b, d)
            out_p[i, j] = K_COUL * p / (2 * cols_hl[j])
            out_e[i, j] = K_COUL * e / (2 * cols_hl[j])


# ----------------------------------------------------------------------------- swarm helpers
class TipPhysics:
    def __init__(self, n_gas, e_m_td, c_v, r_min, r_max=np.inf):
        self.sw = Swarm(SWARM_FILE, n_gas)
        self.N = n_gas
        self.e_m = e_m_td * TD * n_gas
        self.c_v = c_v
        self.r_min = r_min
        self.r_max = r_max
        en = np.exp(self.sw.xf)
        e = en * TD * n_gas
        mu, _, nu = self.sw.at(en)
        alpha = nu / (mu * e)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (alpha[1:] + alpha[:-1]) * np.diff(e))])
        self._e_tab, self._nch_tab = e, EPS0 / QE * cum

    def nu(self, e):
        return self.sw.at(np.abs(e) / (self.N * TD))[2]

    def mu(self, e):
        return self.sw.at(np.abs(e) / (self.N * TD))[0]

    def tip(self, du, e_a=0.0, r_ch=0.0):
        """Tip radius, field and speed: r = min(r_max, max(r_ch, r_min, |dU|/E_m));
        E_tip = |E_a| + min(|dU|/r, E_m) (applied field plus a space-charge excess that saturates at the
        calibrated E_m once the head radius reaches r_max); v = c_v r nu(E_tip)."""
        du = np.abs(du)
        r = np.minimum(np.maximum(np.maximum(self.r_min, r_ch), du / self.e_m), self.r_max)
        et = np.abs(e_a) + np.minimum(du / r, self.e_m)
        return r, et, self.c_v * r * self.nu(et)

    def n_channel(self, e_tip):
        return np.interp(np.abs(e_tip), self._e_tab, self._nch_tab)


# ----------------------------------------------------------------------------- model
class Rod3D:
    def __init__(self, cfg: RodConfig):
        self.cfg = cfg
        self.tp = TipPhysics(cfg.n_gas, cfg.e_m_td, cfg.c_v, cfg.r_min, cfg.r_max)
        self.omega = 2 * math.pi * cfg.f_hz
        self.y0, self.z0 = B_WG / 2, D_CAV / 2
        cap = 256
        self.pos = np.zeros((cap, 3)); self.axis = np.zeros(cap, np.int64); self.hl = np.zeros(cap)
        self.rad = np.zeros(cap); self.lam = np.zeros(cap); self.q = np.zeros(0)
        self.kind = np.zeros(cap, np.int64); self.rod = np.full(cap, -1, np.int64)
        self.born_t = np.zeros(cap)
        self.M = 0
        self.P = np.zeros((0, 0))
        self.edges = []           # [i, j(-1 wall), length, dy]
        self.t = 0.0
        self.e0 = cfg.e0_v_m
        self.wall_q = 0.0
        self.w_ext = self.w_ohm = 0.0
        self.u0 = 0.0
        self.rods = []            # dict per rod
        # laser line
        xs = np.arange(-cfg.x_half, cfg.x_half + 1e-12, cfg.dx_c) + cfg.x_focus_m
        st = Setup(p_pa=cfg.p_pa, axis="x", focus_m=cfg.x_focus_m, energy_j=cfg.energy_j)
        zz = np.linspace(xs[0] - cfg.dx_c / 2, xs[-1] + cfg.dx_c / 2, 40 * len(xs) + 1)
        lam_s = seed_per_length(st, cfg.energy_j, zz)
        counts = np.add.reduceat(0.5 * (lam_s[1:] + lam_s[:-1]) * np.diff(zz), np.arange(0, len(zz) - 1, 40))
        rng = np.random.default_rng(cfg.rng_seed)
        if cfg.stochastic:
            k = rng.poisson(counts)
            counts = np.where(k > 0, rng.gamma(np.maximum(k, 1), 1.0), 0.0)
        self.seed_electrons = counts
        self.laser_nodes = []
        for x, c in zip(xs, counts):
            self.laser_nodes.append(self._add_node((A_WG / 2 + x, self.y0, self.z0), 0, cfg.dx_c / 2, cfg.r_s,
                                                   c / cfg.dx_c, kind=0))
        for a_, b_ in zip(self.laser_nodes[:-1], self.laser_nodes[1:]):
            self.edges.append([a_, b_, cfg.dx_c, 0.0])
        self.born = np.zeros(len(xs), bool)
        self.F = np.zeros((len(xs), 0))
        self._rebuild()

    # ---- bookkeeping
    def _grow(self):
        cap = 2 * len(self.lam)
        for name in ("pos", "axis", "hl", "rad", "lam", "kind", "rod", "born_t"):
            old = getattr(self, name)
            new = np.zeros((cap,) + old.shape[1:], old.dtype)
            if name == "rod":
                new[:] = -1
            new[: len(old)] = old
            setattr(self, name, new)

    def _add_node(self, p, axis, hl, r, lam, kind=1, rod=-1):
        if self.M >= len(self.lam):
            self._grow()
        i = self.M
        self.pos[i] = p; self.axis[i] = axis; self.hl[i] = hl; self.rad[i] = r; self.lam[i] = lam
        self.kind[i] = kind; self.rod[i] = rod; self.born_t[i] = self.t
        self.M += 1
        return i

    def _kernel(self, rows, cols):
        c = self.cfg
        nx, ny, nz = c.n_img
        rp = np.ascontiguousarray(self.pos[rows]); rr = np.ascontiguousarray(self.rad[rows])
        cp = np.ascontiguousarray(self.pos[cols]); ca = np.ascontiguousarray(self.axis[cols])
        ch = np.ascontiguousarray(self.hl[cols])
        op = np.zeros((len(rows), len(cols))); oe = np.zeros((len(rows), len(cols)))
        _block(rp, rr, cp, ca, ch, c.walls, nx, ny, nz, A_WG, B_WG, D_CAV, op, oe)
        return op, oe

    def _rebuild(self):
        idx = np.arange(self.M)
        self.P, _ = self._kernel(idx, idx)
        self.P = 0.5 * (self.P + self.P.T)
        lp = np.array(self.laser_nodes)
        _, self.F = self._kernel(lp, idx)
        self.q = np.zeros(self.M)

    def _extend(self, n_new):
        m0 = self.M - n_new
        old, new = np.arange(m0), np.arange(m0, self.M)
        allx = np.arange(self.M)
        p_new, _ = self._kernel(allx, new)
        p_row, _ = self._kernel(new, old)
        P = np.zeros((self.M, self.M))
        P[:m0, :m0] = self.P
        P[:, m0:] = p_new
        P[m0:, :m0] = p_row
        # symmetrise the new blocks (thin-wire reciprocity holds only approximately)
        P[m0:, :m0] = 0.5 * (P[m0:, :m0] + P[:m0, m0:].T)
        P[:m0, m0:] = P[m0:, :m0].T
        P[m0:, m0:] = 0.5 * (P[m0:, m0:] + P[m0:, m0:].T)
        self.P = P
        lp = np.array(self.laser_nodes)
        _, f_new = self._kernel(lp, new)
        self.F = np.hstack([self.F, f_new])
        self.q = np.concatenate([self.q, np.zeros(n_new)])

    # ---- physics
    def e_applied(self, x, t):
        return self.e0 * np.sin(math.pi * x / A_WG) * math.sin(self.omega * t + self.cfg.rf_phase)

    def _birth(self, k):
        c = self.cfg
        li = self.laser_nodes[k]
        x, y, z = self.pos[li]
        lam = c.n_inc * math.pi * c.r_s ** 2
        rid = len(self.rods)
        # the two root elements span [y0, y0+h] and [y0-h, y0]; no element coincides with the laser node
        ui = self._add_node((x, y + c.h / 2, z), 1, c.h / 2, c.r_s, lam, rod=rid)
        di = self._add_node((x, y - c.h / 2, z), 1, c.h / 2, c.r_s, lam, rod=rid)
        self.edges += [[li, ui, c.h / 2, c.h / 2], [di, li, c.h / 2, c.h / 2]]
        self.rods.append(dict(x=float(x - A_WG / 2), laser=int(li), born_s=self.t, up=ui, down=di, s_up=0.0,
                              s_down=0.0, wall_up=False, wall_down=False, bridged_s=None))
        self.born[k] = True
        self._extend(2)

    def _advance_tip(self, rod, side, v, dt, r_tip, e_tip):
        c = self.cfg
        key, wkey = ("s_" + side), ("wall_" + side)
        if rod[wkey]:
            return
        rod[key] += v * dt
        sgn = 1.0 if side == "up" else -1.0
        added = 0
        while True:
            tip = rod[side]
            y = self.pos[tip, 1]
            gap = (B_WG - y) if side == "up" else y        # from element centre to the wall
            if rod[key] >= gap - c.h / 2:
                self.edges.append([tip, -1, gap, sgn * gap])
                rod[wkey] = True
                if rod["wall_up"] and rod["wall_down"]:
                    rod["bridged_s"] = self.t
                break
            if rod[key] < c.h:
                break
            rod[key] -= c.h
            lam = min(self.tp.n_channel(e_tip), c.n_sat) * math.pi * r_tip ** 2
            ni = self._add_node((self.pos[tip, 0], y + sgn * c.h, self.pos[tip, 2]), 1, c.h / 2,
                                min(max(r_tip, c.r_s), c.r_max), lam, rod=self.rods.index(rod))
            if side == "up":
                self.edges.append([tip, ni, c.h, c.h])
            else:
                self.edges.append([ni, tip, c.h, c.h])
            rod[side] = ni
            added += 1
        if added:
            self._extend(added)

    def step(self, dt):
        c, tp = self.cfg, self.tp
        t1 = self.t + dt
        m = self.M
        ed = np.array(self.edges, float)
        ei = ed[:, 0].astype(np.int64); ej = ed[:, 1].astype(np.int64)
        length, dy = ed[:, 2], ed[:, 3]
        wall = ej < 0
        jj = np.where(wall, 0, ej)
        x_e = self.pos[ei, 0]
        emf = self.e_applied(x_e, t1) * dy
        lam_e = np.where(wall, self.lam[ei], 2 * self.lam[ei] * self.lam[jj] /
                         np.maximum(self.lam[ei] + self.lam[jj], 1e-300))
        phi0 = self.P @ self.q
        e_edge0 = np.abs(phi0[ei] - np.where(wall, 0.0, phi0[jj]) + emf) / length
        g = QE * tp.mu(e_edge0) * lam_e / length
        # implicit: (I + dt K P) q1 = q0 - dt A (g emf)
        D = self.P[ei] - np.where(wall[:, None], 0.0, self.P[jj])
        S = np.eye(m)
        gd = (dt * g)[:, None] * D
        np.add.at(S, ei, gd)
        np.add.at(S, ej[~wall], -gd[~wall])
        rhs = self.q.copy()
        np.add.at(rhs, ei, -dt * g * emf)
        np.add.at(rhs, ej[~wall], dt * g[~wall] * emf[~wall])
        q1 = np.linalg.solve(S, rhs)
        phi = self.P @ q1
        cur = g * (phi[ei] - np.where(wall, 0.0, phi[jj]) + emf)
        self.wall_q += dt * cur[wall].sum()
        p_ext = float(cur @ emf)
        self.w_ext += dt * p_ext
        self.w_ohm += dt * float(np.sum(cur ** 2 / np.maximum(g, 1e-300)))
        self.q = q1
        # ionization along the channel and on the laser line
        e_tot = np.where(g > 0, (np.abs(phi[ei] - np.where(wall, 0.0, phi[jj]) + emf)) / length, 0.0)
        acc = np.zeros(m); cnt = np.zeros(m)
        vert = dy != 0
        np.add.at(acc, ei[vert], e_tot[vert]); np.add.at(cnt, ei[vert], 1)
        np.add.at(acc, ej[vert & ~wall], e_tot[vert & ~wall]); np.add.at(cnt, ej[vert & ~wall], 1)
        e_node = acc / np.maximum(cnt, 1)
        lp = np.array(self.laser_nodes)
        ey_line = self.e_applied(self.pos[lp, 0], t1) + self.F @ q1
        free = ~self.born
        e_node[lp[free]] = np.abs(ey_line[free])
        nu = tp.nu(e_node)
        n_cap = c.n_sat * math.pi * self.rad[:m] ** 2
        lam = np.minimum(self.lam[:m] * np.exp(nu * dt), np.maximum(n_cap, self.lam[:m]))
        self.lam[:m] = lam / (1 + c.alpha_rec * lam / (math.pi * self.rad[:m] ** 2) * dt)
        growing = self.lam[:m] < 0.99 * n_cap
        self.nu_max = float(nu[growing].max()) if growing.any() else 0.0
        # cavity loading
        if c.cavity_feedback:
            dwde = EPS0 * max(self.e0, 1.0) * A_WG * B_WG * D_CAV / 4
            e_set = c.e0_v_m if self.t < c.drive_off_s else 0.0
            self.e0 += dt * ((e_set - self.e0) / c.tau_amp_s - p_ext / dwde)
            self.e0 = max(self.e0, 0.0)
        self.t = t1
        # tips
        tip_info = []
        for rod in self.rods:
            for side in ("up", "down"):
                if rod["wall_" + side]:
                    continue
                du = phi[rod[side]]
                tip_node = rod[side]
                e_a = self.e_applied(self.pos[tip_node, 0], t1)
                r_tip, e_tip, v = tp.tip(du, e_a, self.rad[tip_node])
                # mean field across the remaining tip-wall gap (same convention as the wall edge)
                y = self.pos[tip_node, 1]
                gap = max((B_WG - y) if side == "up" else y, c.h / 2)
                e_gap = abs(du + e_a * (gap if side == "up" else -gap)) / gap
                if e_gap > e_tip:
                    e_tip, v = e_gap, tp.c_v * r_tip * float(tp.nu(e_gap))
                tip_info.append((rod, side, float(v), float(r_tip), float(e_tip), float(du)))
        for rod, side, v, r_tip, e_tip, _ in tip_info:
            self._advance_tip(rod, side, v, dt, r_tip, e_tip)
        # inception; sites within one channel diameter (2 r_max) of an existing rod merge into it
        n_line = self.lam[lp] / (math.pi * c.r_s ** 2)
        cand = np.nonzero((~self.born) & (n_line >= c.n_inc))[0]
        for k in cand[np.argsort(-n_line[cand])]:
            xk = self.pos[lp[k], 0] - A_WG / 2
            if any(abs(xk - r["x"]) < 2 * c.r_max for r in self.rods):
                continue
            self._birth(int(k))
        self.last = dict(p_ext=p_ext, phi=phi, tips=tip_info, ey_line=ey_line)
        vmax = max([ti[2] for ti in tip_info], default=0.0)
        return vmax

    def energy_residual(self):
        u = 0.5 * float(self.q @ self.P @ self.q)
        return (self.w_ext - self.w_ohm - (u - self.u0)) / max(abs(self.w_ext), 1e-300)


# ----------------------------------------------------------------------------- driver
def run(cfg: RodConfig, out_root: Path = RESULTS, save=True, verbose=True) -> dict:
    model = Rod3D(cfg)
    period = 1 / cfg.f_hz
    dt0 = period / cfg.steps_per_period
    frames, next_frame = [], 0.0
    hist = {k: [] for k in ("t_s", "e0_v_m", "p_ext_w", "n_rods", "max_len_m", "max_du_v", "line_ey_v_m",
                            "max_v_m_s")}
    t_bridge = None
    wall0 = time.time()
    vmax = 0.0
    while model.t < cfg.t_end_s:
        dt = dt0
        if vmax > 0:
            dt = min(dt, 0.25 * cfg.h / vmax)
        if getattr(model, "nu_max", 0.0) > 0:
            dt = min(dt, cfg.max_ion_step / model.nu_max)
        dt = max(dt, dt0 / 256)
        vmax = model.step(dt)
        lens = [model.pos[r["up"], 1] - model.pos[r["down"], 1] + cfg.h for r in model.rods]
        if model.t >= next_frame - 1e-15:
            next_frame += cfg.frame_dt_s
            hist["t_s"].append(model.t); hist["e0_v_m"].append(model.e0)
            hist["p_ext_w"].append(model.last["p_ext"]); hist["n_rods"].append(len(model.rods))
            hist["max_len_m"].append(max(lens, default=0.0))
            hist["max_du_v"].append(max([abs(ti[5]) for ti in model.last["tips"]], default=0.0))
            hist["line_ey_v_m"].append(float(np.abs(model.last["ey_line"]).max()))
            hist["max_v_m_s"].append(vmax)
            frames.append(dict(t=model.t, m=model.M, lam=model.lam[: model.M].copy(), q=model.q.copy()))
            if verbose and len(frames) % 20 == 0:
                print(f"  t={model.t*1e9:6.2f} ns  E0={model.e0/1e5:6.1f} kV/cm  rods={len(model.rods)}  "
                      f"Lmax={hist['max_len_m'][-1]*1e3:5.2f} mm  vmax={vmax:.2e}  M={model.M}  "
                      f"wall={time.time()-wall0:.0f}s", flush=True)
        bridged = [r["bridged_s"] for r in model.rods if r["bridged_s"] is not None]
        if bridged and t_bridge is None:
            t_bridge = min(bridged)
        if t_bridge is not None and model.t > t_bridge + cfg.stop_after_bridge_s:
            break
    summary = summarize(cfg, model, hist, t_bridge, time.time() - wall0)
    if save:
        out = out_root / cfg.name
        out.mkdir(parents=True, exist_ok=True)
        m = model.M
        lam_frames = np.zeros((len(frames), m), np.float32)
        for k, f in enumerate(frames):
            lam_frames[k, : f["m"]] = f["lam"]
        ed = np.array(model.edges, float)
        np.savez_compressed(out / "rod3d.npz", t_s=np.array([f["t"] for f in frames]),
                            m_frame=np.array([f["m"] for f in frames]), lam=lam_frames,
                            pos=model.pos[:m], axis=model.axis[:m], rad=model.rad[:m], kind=model.kind[:m],
                            rod=model.rod[:m], born_t=model.born_t[:m], edges=ed,
                            seed_electrons=model.seed_electrons,
                            **{"hist_" + k: np.array(v) for k, v in hist.items()})
        (out / "summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    if verbose:
        print(json.dumps({k: summary[k] for k in ("name", "t_first_birth_ns", "n_rods", "t_bridge_ns",
                                                  "n_long_rods", "e0_min_kv_cm")}), flush=True)
    return summary


def summarize(cfg, model, hist, t_bridge, wall):
    c = cfg
    rods = []
    for r in model.rods:
        length = model.pos[r["up"], 1] - model.pos[r["down"], 1] + c.h
        rods.append(dict(x_mm=1e3 * r["x"], born_ns=1e9 * r["born_s"], length_mm=1e3 * float(length),
                         bridged_ns=None if r["bridged_s"] is None else 1e9 * r["bridged_s"],
                         top_wall=r["wall_up"], bottom_wall=r["wall_down"]))
    long_rods = [r for r in rods if r["length_mm"] >= 0.25 * B_WG * 1e3]
    return dict(name=c.name, config=asdict(c), wall_s=wall, nodes=int(model.M),
                n_seed_electrons=float(model.seed_electrons.sum()),
                t_first_birth_ns=min([r["born_ns"] for r in rods], default=None),
                n_rods=len(rods), n_long_rods=len(long_rods),
                long_rod_x_mm=[r["x_mm"] for r in long_rods],
                t_bridge_ns=None if t_bridge is None else 1e9 * t_bridge,
                n_bridged=sum(r["bridged_ns"] is not None for r in rods),
                e0_min_kv_cm=float(min(hist["e0_v_m"], default=c.e0_v_m)) / 1e5,
                energy_from_cavity_j=model.w_ext, ohmic_j=model.w_ohm,
                energy_residual_rel=model.energy_residual(),
                charge_residual_c=float(model.q.sum() + model.wall_q),
                rods=rods)


# ----------------------------------------------------------------------------- calibration
def calibrate_tip_law(runs, e_m_grid=(1, 2, 3, 5, 8, 10, 20, 50, 100, 300, 1000), t_start_s=None):
    """Fit c_v and E_m of v = c_v r_tip nu(E_tip) to axisymmetric streamer2d runs.

    The 1e20 front displacement over each RF period is compared with the integral of
    r_tip nu(E_tip) over the same period (using the resolved tip potential, applied field
    and channel radius). Runs are weighted equally; returns the E_m minimising the rms
    log residual and the per-run bias factor."""
    from .streamer2d import CASES
    t_start_s = t_start_s or {}
    data = []
    for run_name in runs:
        h = np.load(ROOT / "results" / "streamer2d" / run_name / "history.npz")
        if "phi_tip_plus_V" not in h.files:
            continue
        cfg = CASES[run_name.replace("_z8", "")]
        t = h["t_s"]
        t0 = t_start_s.get(run_name, 0.3e-9)
        for side, sgn in (("plus", 1), ("minus", -1)):
            z = sgn * h[f"z_{side}_1e20_m"]
            du = h[f"phi_tip_{side}_V"]
            ok = np.isfinite(z) & np.isfinite(du) & np.isfinite(h["r_1e20_mid_m"]) & (t > t0)
            data.append((run_name, cfg.N, t[ok], z[ok], du[ok], h["E_app_V_m"][ok], h["r_1e20_mid_m"][ok],
                         1.0 / cfg.freq_Hz))
    best = None
    for em in e_m_grid:
        res = {}
        for run_name, n, t, z, du, ea, rch, period in data:
            if len(t) < 3:
                continue
            _, _, g = TipPhysics(n, em, 1.0, 10e-6, 0.3e-3).tip(du, ea, rch)
            gi = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(t))])
            edges = np.arange(t[0], t[-1] + 1e-15, period)
            dz = np.diff(np.interp(edges, t, z))
            dg = np.diff(np.interp(edges, t, gi))
            ok = (dz > 1e-6) & (dg > 0)
            res.setdefault(run_name, []).extend(np.log(dz[ok] / dg[ok]))
        res = {k: np.asarray(v) for k, v in res.items() if len(v)}
        if not res:
            continue
        lc = float(np.mean([v.mean() for v in res.values()]))
        rms = float(np.sqrt(np.mean([np.mean((v - lc) ** 2) for v in res.values()])))
        if best is None or rms < best["rms_log"]:
            best = dict(e_m_td=em, c_v=math.exp(lc), rms_log=rms,
                        per_run_ratio={k: float(np.exp(v.mean() - lc)) for k, v in res.items()},
                        n_points=int(sum(len(v) for v in res.values())))
    return best


# ----------------------------------------------------------------------------- static check
def static_dipole(L, r, h, e0=1e6):
    """Equipotential straight conductor of half-length L in a uniform field, free space.
    Returns (dipole moment, slender-cylinder analytic 4 pi eps0 L^3 E0 / (3 (ln(4L/r) - 7/3)))."""
    n = int(round(2 * L / h))
    ys = -L + h * (np.arange(n) + 0.5)
    pos = np.column_stack([np.zeros(n), ys, np.zeros(n)])
    op = np.zeros((n, n)); oe = np.zeros((n, n))
    _block(pos, np.full(n, r), pos, np.ones(n, np.int64), np.full(n, h / 2), False, 0, 0, 0, 1.0, 1.0, 1.0, op, oe)
    P = 0.5 * (op + op.T)
    # phi_sc + phi_ext = V (const), phi_ext = -E0 y, sum q = 0
    A = np.zeros((n + 1, n + 1))
    A[:n, :n] = P; A[:n, n] = -1; A[n, :n] = 1
    rhs = np.concatenate([e0 * ys, [0.0]])
    sol = np.linalg.solve(A, rhs)
    p = float(sol[:n] @ ys)
    analytic = 4 * math.pi * EPS0 * L ** 3 * e0 / (3 * (math.log(4 * L / r) - 7.0 / 3.0))
    return p, analytic


# ----------------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser(prog="rod3d")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("calibrate")
    c.add_argument("runs", nargs="+")
    r = sub.add_parser("run")
    r.add_argument("--name", required=True)
    r.add_argument("--e0-kv-cm", type=float, default=200.0)
    r.add_argument("--focus-mm", type=float, default=0.0)
    r.add_argument("--energy-mj", type=float, default=0.75)
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--h-mm", type=float, default=0.5)
    r.add_argument("--dx-mm", type=float, default=0.25)
    r.add_argument("--t-end-ns", type=float, default=20.0)
    r.add_argument("--phase", type=float, default=0.0)
    r.add_argument("--no-feedback", action="store_true")
    r.add_argument("--drive-off-ns", type=float, default=math.inf)
    r.add_argument("--frame-ps", type=float, default=50.0)
    r.add_argument("--deterministic", action="store_true")
    r.add_argument("--calibration", type=Path, default=RESULTS / "tip_calibration.json")
    a = ap.parse_args(argv)
    if a.cmd == "calibrate":
        best = calibrate_tip_law(a.runs)
        RESULTS.mkdir(parents=True, exist_ok=True)
        (RESULTS / "tip_calibration.json").write_text(json.dumps(dict(best, runs=a.runs), indent=1), encoding="utf-8")
        print(json.dumps(best, indent=1))
        return 0
    cal = json.loads(a.calibration.read_text(encoding="utf-8")) if a.calibration.exists() else {}
    cfg = RodConfig(name=a.name, e0_v_m=a.e0_kv_cm * 1e5, x_focus_m=a.focus_mm * 1e-3, energy_j=a.energy_mj * 1e-3,
                    rng_seed=a.seed, h=a.h_mm * 1e-3, dx_c=a.dx_mm * 1e-3, t_end_s=a.t_end_ns * 1e-9,
                    rf_phase=a.phase, cavity_feedback=not a.no_feedback, stochastic=not a.deterministic,
                    drive_off_s=a.drive_off_ns * 1e-9, frame_dt_s=a.frame_ps * 1e-12)
    for k in ("c_v", "e_m_td", "n_inc"):
        if k in cal:
            cfg = replace(cfg, **{k: cal[k]})
    run(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
