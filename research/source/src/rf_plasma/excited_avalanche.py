"""Radial reaction-diffusion avalanche with excited-state reservoirs (stepwise ionization, pooling, excimers).

Purpose: test whether ionization of stored electronic excitation (cumulative / stepwise ionization) can change the
He / Ar / N2 ordering of laser-seeded and RF-only breakdown fields at 50-300 psig, 2.95 GHz. It extends the
direct-ionization avalanche of gas_compare (same effective-field closure, field envelope and N_br criterion).

Model (per radial cell; sphere for a point electron, cylinder of fixed length L for a laser line seed):
    dn_e/dt = div(D_s grad n_e) + n_e N k_iz + sum_r n_e n_r k_si,r + sum_p y_p k_p n_a n_b - alpha n_e n_M+
    dn_A+/dt = n_e N k_iz + (atomic stepwise, pooling) - k_conv N^2 n_A+
    dn_M+/dt = k_conv N^2 n_A+ + (molecular stepwise, pooling) - alpha n_e n_M+
    dn_r/dt = n_e N sum_{g in r} k_g + (three-body / recombination feed) - n_r [A_r + k2 N + k3 N^2 + n_e(k_si + k_sup)]
              - (pooling consumption)
* Electron rates (k_iz, group excitation counts k_g, D_e = DN/N, mean energy) from the two-term BOLOS tables at
  the RF effective field E_eff = E_rms nu_m / sqrt(nu_m^2 + omega^2) (local-field approximation; the EEDF is not
  modified by superelastic or electron-electron collisions -- a stated limitation).
* Rates of excited species: Maxwellian fits A Te^b exp(-c/Te) at Te = (2/3)<eps> (stated limitation) or constants.
* Transport: electrons only; ions and excited species are immobile on the simulated time scale. Free-to-ambipolar
  transition D_s = D_a + (D_e - D_a) / (1 + (sigma/lambda_D)^2 / c_ar) per face (Allis-Rose-type interpolation;
  c_ar is a sensitivity parameter), with sigma the instantaneous rms cloud radius per dimension.
* Pooling convention: reaction rate R = k n_a n_b (a = b: R = k n_a^2, consuming 2 n_a per event).
Breakdown when the total electron number reaches N_br (default 1e12, as gas_compare). The early stochastic phase is
handled by scaling the initial electron number with the Yule factor W (exact for linear growth).
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp
from scipy.optimize._numdiff import approx_derivative, group_columns

from .swarm import SWARM_DIR, read_table
from . import swarm_excited as sx

ROOT = Path(__file__).resolve().parents[2]
KINETICS_FILE = ROOT / "cases" / "excited_kinetics.json"
KB = 1.380649e-23
QE = 1.602176634e-19
EPS0 = 8.8541878128e-12
TD = 1e-21
N_LOSCHMIDT = 2.686780111e25
PSI = 6894.757
ATM = 101325.0


# ----------------------------------------------------------------------------- electron swarm (LFA)
class Swarm:
    """Log-interpolated two-term BOLOS data for one gas: mean energy, k_m, k_iz, group excitation counts, N*D_e."""

    def __init__(self, gas: str):
        base = read_table(SWARM_DIR / f"{gas}_swarm.csv")
        grp = sx.load(gas)
        self.gas = gas
        self.groups = sx.GROUP_COLUMNS[gas]
        en = np.asarray(base["E_N_Td"])
        o = np.argsort(en)
        self.ln_en = np.log(en[o])
        self.eps = np.asarray(base["mean_energy_eV"])[o]
        self.km = np.asarray(base["km_m3_s"])[o]
        self.dn = np.asarray(base["DN_1_ms"])[o]
        g_en = np.asarray(grp["E_N_Td"])
        if not np.allclose(np.sort(g_en), en[o], rtol=1e-5):
            raise ValueError("group table and swarm table use different E/N grids")
        og = np.argsort(g_en)
        self.ln_kiz = np.log(np.maximum(np.asarray(grp["kiz_m3_s"])[og], 1e-80))
        self.ln_kai = np.log(np.maximum(np.asarray(grp["k_ai"])[og], 1e-80))
        self.ln_kg = {g: np.log(np.maximum(np.asarray(grp[g])[og], 1e-80)) for g in self.groups}

    def _x(self, en_td):
        return np.clip(np.log(max(en_td, 1e-6)), self.ln_en[0], self.ln_en[-1])

    def effective_en(self, e0_v_m: float, n: float, f_hz: float) -> float:
        e_rms = e0_v_m / math.sqrt(2.0)
        nu_m = n * float(np.interp(self._x(e_rms / n / TD), self.ln_en, self.km))
        w = 2 * math.pi * f_hz
        return e_rms * nu_m / math.sqrt(nu_m * nu_m + w * w) / n / TD

    def at(self, en_td: float) -> dict:
        x = self._x(en_td)
        out = dict(eps=float(np.interp(x, self.ln_en, self.eps)),
                   dn=float(np.interp(x, self.ln_en, self.dn)),
                   kiz=float(np.exp(np.interp(x, self.ln_en, self.ln_kiz))),
                   k_ai=float(np.exp(np.interp(x, self.ln_en, self.ln_kai))))
        for g in self.groups:
            out[g] = float(np.exp(np.interp(x, self.ln_en, self.ln_kg[g])))
        return out


# ----------------------------------------------------------------------------- kinetics specification
def rate_fit(spec, te_ev: float, te_k: float | None = None) -> float:
    """spec: None | number | {"A":, "b":, "c":} -> A Te^b exp(-c/Te) (Te in eV)
    | {"a300":, "n":} -> a300 (Te_K/300)^-n (recombination form)."""
    if spec is None:
        return 0.0
    if isinstance(spec, (int, float)):
        return float(spec)
    if "a300" in spec:
        tk = te_k if te_k is not None else np.asarray(te_ev) * QE / KB
        return float(spec["a300"]) * (tk / 300.0) ** (-float(spec["n"]))
    te = np.asarray(te_ev)
    return float(spec["A"]) * te ** float(spec.get("b", 0.0)) * np.exp(-float(spec.get("c", 0.0)) / te)


def load_kinetics(path: Path = KINETICS_FILE) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def scaled_kinetics(kin: dict, which: str) -> dict:
    """Return a copy of a gas spec with every uncertain coefficient multiplied by its uncertainty factor ('low' ->
    1/f for ionizing channels and f for loss channels; 'high' the reverse). which in {'nominal','low','high'}."""
    if which == "nominal":
        return json.loads(json.dumps(kin))
    k = json.loads(json.dumps(kin))
    up = which == "high"

    def sc(spec, f, ionizing):
        if spec is None or f in (None, 1.0):
            return spec
        m = f if (ionizing == up) else 1.0 / f
        if isinstance(spec, (int, float)):
            return spec * m
        spec = dict(spec)
        if "a300" in spec:
            spec["a300"] *= m
        else:
            spec["A"] *= m
        return spec

    for r in k["reservoirs"]:
        uf = r.get("uncertainty", {})
        r["k_si"] = sc(r.get("k_si"), uf.get("k_si"), True)
        for key in ("k_sup", "k2", "k3", "A_s"):
            if key in r:
                r[key] = sc(r[key], uf.get(key), False)
    for p in k.get("pooling", []):
        p["k"] = sc(p["k"], p.get("uncertainty"), True)
    return k


@dataclass
class Case:
    gas: str = "Ar"
    psig: float = 200.0
    t_k: float = 300.0
    f_hz: float = 2.95e9
    e0_v_m: float = 50e5             # peak RF field amplitude (E0 cos wt)
    geometry: str = "cylinder"       # 'cylinder' (laser line seed) | 'sphere' (single electron)
    n0: float = 376.0                # initial electrons (already including the Yule factor W)
    sigma0_m: float = 4.1e-6         # initial Gaussian rms radius per dimension
    length_m: float = 2.03e-3        # cylinder length (seed line)
    t_flat_s: float = 100e-9         # flat RF after t = 0
    tau_s: float = 35.7e-9           # ring-down amplitude time constant after t_flat
    n_tau_after: float = 10.0
    n_br: float = 1e12
    criterion: str = "n_br"          # 'n_br' (total electrons) | 'screen' (peak n_e reaches RF-screening density)
    excited: bool = True             # False -> direct ionization only (gas_compare limit)
    nonlinear: bool = True           # stepwise / pooling / e-transfer channels (quadratic in densities)
    hm_branch: float = 0.0           # Hornbeck-Molnar: fraction of AI-eligible excitations giving X2+ + e
    eta_gs: float = 0.0              # photoemission feedback: f_geo * P_survival (x reservoir photoyield Y)
    ee: bool | str = False           # k_iz at the local ionization degree: True (swarm_ee tables) | 'dense'
    c_ar: float = 10.0               # free-ambipolar transition parameter
    cells: int = 120
    rtol: float = 1e-5
    kinetics: dict | None = None     # per-gas spec (defaults to KINETICS_FILE[gas])

    @property
    def p_pa(self):
        return self.psig * PSI + ATM

    @property
    def n_gas(self):
        return self.p_pa / (KB * self.t_k)

    @property
    def t_end(self):
        return self.t_flat_s + self.n_tau_after * self.tau_s

    def envelope(self, t):
        return 1.0 if t <= self.t_flat_s else math.exp(-(t - self.t_flat_s) / max(self.tau_s, 1e-30))


# ----------------------------------------------------------------------------- grid
def make_grid(case: Case, d_max: float):
    t_size = case.t_flat_s + 3 * case.tau_s
    r_max = 7.0 * math.sqrt(case.sigma0_m ** 2 + 2 * d_max * t_size) + 3 * case.sigma0_m
    m = case.cells
    dr0 = case.sigma0_m / 6.0
    # geometric stretching with first cell ~ dr0
    if r_max / m <= dr0:
        faces = np.linspace(0.0, r_max, m + 1)
    else:
        lo, hi = 1e-9, 50.0
        for _ in range(200):
            beta = 0.5 * (lo + hi)
            first = r_max * math.expm1(beta / m) / math.expm1(beta)
            if first > dr0:
                lo = beta
            else:
                hi = beta
        faces = r_max * np.expm1(beta * np.arange(m + 1) / m) / math.expm1(beta)
    rc = 0.5 * (faces[1:] + faces[:-1])
    if case.geometry == "sphere":
        vol = 4.0 / 3.0 * math.pi * (faces[1:] ** 3 - faces[:-1] ** 3)
        area = 4.0 * math.pi * faces ** 2
    else:
        vol = math.pi * (faces[1:] ** 2 - faces[:-1] ** 2) * case.length_m
        area = 2.0 * math.pi * faces * case.length_m
    return faces, rc, vol, area




# ----------------------------------------------------------------------------- e-e (ionization-degree) tables
class EETable:
    """k_iz(E/N, x = ne/N) from the e-e two-term Boltzmann tables (rf_plasma.swarm_ee, Hagelaar & Pitchford 2005
    sec. 2.6). Rows where the nonlinear solve did not converge (high x, low E/N) are filled with the value of the
    next-lower x level (a conservative lower bound: Maxwellization raises the tail). N enters only via ln Lambda,
    so one table set (built at N_table) is used for all pressures."""

    def __init__(self, gas: str, n_table: float = 2.414e26):
        from .swarm_ee import X_LEVELS, table_path
        base = read_table(SWARM_DIR / f"{gas}_swarm.csv")
        self.ln_en = np.log(np.sort(np.asarray(base["E_N_Td"])))
        rows = []
        for x in X_LEVELS:
            d = read_table(table_path(gas, n_table, float(x)))
            en = np.asarray(d["E_N_Td"])
            k = np.asarray(d["kiz_m3_s"])
            ok = k > 0
            lk = np.full(self.ln_en.shape, np.nan)
            for i, le in enumerate(self.ln_en):
                j = np.nonzero(np.isclose(np.log(en[ok]), le, atol=1e-6))[0]
                if j.size:
                    lk[i] = math.log(k[ok][j[0]])
            rows.append(lk)
        tab = np.array(rows)
        tab[0] = np.where(np.isfinite(tab[0]), tab[0], math.log(1e-80))
        for i in range(1, len(tab)):
            tab[i] = np.where(np.isfinite(tab[i]), np.maximum(tab[i], tab[i - 1]), tab[i - 1])
        self.tab = tab
        self.lx = np.log(np.maximum(np.asarray(X_LEVELS, float), 1e-30))
        self.lx[0] = math.log(1e-8)                 # x = 0 level placed at 1e-8 for interpolation

    def kiz(self, en_td: float, x: np.ndarray) -> np.ndarray:
        le = np.clip(math.log(max(en_td, 1e-6)), self.ln_en[0], self.ln_en[-1])
        col = np.array([np.interp(le, self.ln_en, r) for r in self.tab])
        lx = np.clip(np.log(np.maximum(x, 1e-30)), self.lx[0], self.lx[-1])
        return np.exp(np.interp(lx, self.lx, col))
class DenseTable(EETable):
    """k_iz(E/N, x) from results/excited_avalanche/dense_regime.json (e-e + Coulomb e-ion momentum transfer,
    10 E/N points 3-50 Td, x = 0..1e-2). Missing (non-converged) points take the next-lower x value; outside
    3-50 Td the edge column is used."""

    def __init__(self, gas: str, path: Path | None = None):
        d = json.loads(Path(path or ROOT / "results" / "excited_avalanche" / "dense_regime.json").read_text())
        en = np.array(sorted(d["en_td"]))
        xs = list(d["x"])
        tab = np.full((len(xs), len(en)), np.nan)
        etab = np.full((len(xs), len(en)), np.nan)
        for r in d["results"]:
            if r["gas"] != gas or (r["x"] > 0 and not r["ei"]):
                continue
            i = xs.index(r["x"])
            for row in r["rows"]:
                if row["converged"] and row["kiz"] > 0:
                    j = int(np.argmin(abs(en - row["E_N_Td"])))
                    tab[i, j] = math.log(row["kiz"])
                    etab[i, j] = row["mean_energy_eV"]
        tab[0] = np.where(np.isfinite(tab[0]), tab[0], math.log(1e-80))
        etab[0] = np.where(np.isfinite(etab[0]), etab[0], np.nanmin(etab[0]))
        self.filled = ~np.isfinite(tab)
        for i in range(1, len(tab)):
            fill = ~np.isfinite(tab[i])
            tab[i] = np.where(fill, tab[i - 1], tab[i])
            etab[i] = np.where(fill | ~np.isfinite(etab[i]), etab[i - 1], etab[i])
        self.ln_en = np.log(en)
        self.tab = tab
        self.etab = etab
        self.lx = np.log(np.maximum(np.asarray(xs, float), 1e-30))
        self.lx[0] = math.log(1e-8)
        self.max_x = 0.0
        self.used_filled = False

    def kiz(self, en_td: float, x: np.ndarray) -> np.ndarray:
        xm = float(np.max(x)) if np.size(x) else 0.0
        if xm > self.max_x:
            self.max_x = xm
        if not self.used_filled and xm > 0:
            le = np.clip(math.log(max(en_td, 1e-6)), self.ln_en[0], self.ln_en[-1])
            j = int(np.searchsorted(self.ln_en, le))
            cols = [c for c in (j - 1, j) if 0 <= c < len(self.ln_en)]
            i = int(np.searchsorted(self.lx, math.log(max(xm, 1e-30))))
            rows = [r for r in (i - 1, i) if 0 <= r < len(self.lx)]
            if self.filled[np.ix_(rows, cols)].any():
                self.used_filled = True
        return super().kiz(en_td, x)

    def te(self, en_td: float, x: np.ndarray) -> np.ndarray:
        le = np.clip(math.log(max(en_td, 1e-6)), self.ln_en[0], self.ln_en[-1])
        col = np.array([np.interp(le, self.ln_en, r) for r in self.etab])
        lx = np.clip(np.log(np.maximum(x, 1e-30)), self.lx[0], self.lx[-1])
        return np.maximum(2.0 * np.interp(lx, self.lx, col) / 3.0, 0.03)


# ----------------------------------------------------------------------------- solver
class Avalanche:
    """Radial avalanche (cells) + a 0-D linear population of wall-photoemission avalanches (eta_gs > 0)."""

    COUNTERS = ["ion_direct", "ion_step", "ion_pool", "ion_hm", "recomb", "exc_made", "exc_lost_rad",
                "exc_lost_quench", "exc_lost_sup", "exc_to_ion", "outflow_e", "photo_e", "wall_ion"]

    def __init__(self, case: Case, swarm: Swarm | None = None):
        self.c = case
        self.sw = swarm or Swarm(case.gas)
        kin = case.kinetics if case.kinetics is not None else load_kinetics()[case.gas]
        self.kin = kin
        self.res = kin["reservoirs"] if case.excited else []
        self.pool = kin.get("pooling", []) if (case.excited and case.nonlinear) else []
        self.ion = kin["ion"]
        self.rnames = [r["name"] for r in self.res]
        self.nr = len(self.res)
        self.nsp = 3 + self.nr                     # e, A+, M+, reservoirs
        n = case.n_gas
        self.N = n
        en0 = self.sw.effective_en(case.e0_v_m, n, case.f_hz)
        d_max = self.sw.at(en0)["dn"] / n
        self.faces, self.rc, self.vol, self.area = make_grid(case, d_max)
        self.m = len(self.rc)
        self.dist = np.diff(self.rc)
        k0 = float(self.ion["K0_cm2_Vs"]) * 1e-4
        self.mu_i = k0 * N_LOSCHMIDT / n
        self.counters = self.COUNTERS
        self.iw = self.m * self.nsp                # start of wall block [Nw, Rw_0..]
        self.ic = self.iw + 1 + self.nr            # start of counters
        self.ee = (DenseTable(case.gas) if case.ee == "dense" else EETable(case.gas)) if case.ee else None
        self._sparsity = self._jac_sparsity()

    def scalars(self, t: float) -> dict:
        c = self.c
        e0 = c.e0_v_m * c.envelope(t)
        en = self.sw.effective_en(e0, self.N, c.f_hz) if e0 > 0 else 1e-3
        s = self.sw.at(en)
        te = max(2.0 * s["eps"] / 3.0, 0.03)
        s["en"], s["te"] = en, te
        s["de"] = s["dn"] / self.N
        s["da"] = self.mu_i * te
        k_tot = sum(s[g] for g in self.sw.groups)
        h = min(c.hm_branch * s["k_ai"] / max(k_tot, 1e-300), 1.0) if c.excited else 0.0
        s["h"] = h
        s["src"] = [self.N * (1.0 - h) * (sum(s[g] for g in r["groups"])
                                          + sum(f * s[g] for g, f in r.get("group_frac", {}).items()))
                    for r in self.res]
        s["nu_hm"] = self.N * c.hm_branch * s["k_ai"] if c.excited else 0.0
        return s

    def _jac_sparsity(self):
        m, ns, nr = self.m, self.nsp, self.nr
        n_tot = self.ic
        S = sp.lil_matrix((n_tot, n_tot), dtype=np.int8)
        for j in range(m):
            b = j * ns
            S[b:b + ns, b:b + ns] = 1
            for jj in (j - 1, j + 1):
                if 0 <= jj < m:
                    S[b, jj * ns] = 1
        S[self.iw:self.ic, self.iw:self.ic] = 1
        # The one-way coupling cloud reservoirs -> wall photon source is left out of the finite-difference
        # sparsity: a dense row would force one column group per unknown, and an approximate Newton Jacobian
        # does not change the converged solution.
        return S.tocsr()

    def initial(self):
        c, m, ns = self.c, self.m, self.nsp
        y = np.zeros(self.ic)
        s2 = c.sigma0_m ** 2
        if c.geometry == "sphere":
            prof = np.exp(-self.rc ** 2 / (2 * s2))
        else:
            prof = np.exp(-self.rc ** 2 / (2 * s2))
        prof = prof / np.sum(prof * self.vol)
        y[0: m * ns: ns] = c.n0 * prof
        y[1: m * ns: ns] = c.n0 * prof
        return y

    def cloud_total(self, y):
        return float(np.sum(y[0: self.m * self.nsp: self.nsp] * self.vol))

    def totals(self, y):
        return self.cloud_total(y) + max(float(y[self.iw]), 0.0)

    def rms_radius(self, ne):
        w = ne * self.vol
        tot = w.sum()
        if tot <= 0:
            return self.c.sigma0_m
        dim = 3.0 if self.c.geometry == "sphere" else 2.0
        return math.sqrt(max(float(np.sum(w * self.rc ** 2) / tot) / dim, 1e-30))

    def _res_linear_loss(self, r):
        N = self.N
        return float(r.get("A_s", 0.0)) + float(r.get("k2", 0.0)) * N + float(r.get("k3", 0.0)) * N * N

    def rhs(self, t, y):
        return self._eval(t, y)[0]

    def counter_rates(self, t, y):
        return self._eval(t, y)[1]

    def _eval(self, t, y):
        c, m, ns, N = self.c, self.m, self.nsp, self.N
        s = self.scalars(t)
        te = s["te"]
        Y = y[: m * ns].reshape(m, ns)
        ne = np.maximum(Y[:, 0], 0.0)
        na = np.maximum(Y[:, 1], 0.0)
        nm = np.maximum(Y[:, 2], 0.0)
        nr = np.maximum(Y[:, 3:], 0.0)
        d = np.zeros_like(Y)
        cnt = np.zeros(len(self.counters))
        # transport of electrons
        sig = self.rms_radius(ne)
        nf = np.maximum(0.5 * (ne[1:] + ne[:-1]), 1.0)
        x = sig * sig / (EPS0 * te / (QE * nf))
        dsf = s["da"] + (s["de"] - s["da"]) / (1.0 + x / c.c_ar)
        flux = -dsf * (ne[1:] - ne[:-1]) / self.dist * self.area[1:-1]
        div = np.zeros(m)
        div[:-1] -= flux
        div[1:] += flux
        out_flux = s["de"] * ne[-1] / (self.faces[-1] - self.rc[-1]) * self.area[-1]
        div[-1] -= out_flux
        d[:, 0] += div / self.vol
        cnt[10] += out_flux
        # direct + Hornbeck-Molnar ionization
        kiz_cell = self.ee.kiz(s["en"], ne / N) if self.ee is not None else s["kiz"]
        te = self.ee.te(s["en"], ne / N) if isinstance(self.ee, DenseTable) else te   # per-cell Te (dense EEDF)
        r_iz = ne * N * kiz_cell
        r_hm = ne * s["nu_hm"]
        d[:, 0] += r_iz + r_hm
        d[:, 1] += r_iz
        d[:, 2] += r_hm
        cnt[0] += np.sum(r_iz * self.vol)
        cnt[3] += np.sum(r_hm * self.vol)
        conv = rate_fit(self.ion.get("k_conv"), te) * N * N * na
        d[:, 1] -= conv
        d[:, 2] += conv
        rec = rate_fit(self.ion.get("alpha_dr"), te) * ne * nm
        d[:, 0] -= rec
        d[:, 2] -= rec
        cnt[4] += np.sum(rec * self.vol)
        dr_to = self.ion.get("dr_product")
        if dr_to in self.rnames:
            d[:, 3 + self.rnames.index(dr_to)] += rec
        # electron-electron-ion (three-body) recombination into excited states, both ion pools
        a3 = rate_fit(self.ion.get("alpha_3b"), te) if self.ion.get("alpha_3b") else 0.0
        if np.any(a3):
            rec3a, rec3m = a3 * ne * ne * na, a3 * ne * ne * nm
            d[:, 0] -= rec3a + rec3m
            d[:, 1] -= rec3a
            d[:, 2] -= rec3m
            cnt[4] += np.sum((rec3a + rec3m) * self.vol)
            if dr_to in self.rnames:
                d[:, 3 + self.rnames.index(dr_to)] += rec3a + rec3m
        # reservoirs (cloud)
        photons = 0.0
        for i, r in enumerate(self.res):
            n_i = nr[:, i]
            src = ne * s["src"][i]
            d[:, 3 + i] += src
            cnt[5] += np.sum(src * self.vol)
            a_s = float(r.get("A_s", 0.0))
            r_rad = a_s * n_i
            r_q = float(r.get("k2", 0.0)) * N * n_i
            r_3 = float(r.get("k3", 0.0)) * N * N * n_i
            d[:, 3 + i] -= r_rad + r_q + r_3
            cnt[6] += np.sum(r_rad * self.vol)
            photons += float(r.get("photo", 0.0)) * float(r.get("Y", 0.0)) * np.sum(r_rad * self.vol)
            for tgt, frac in r.get("k3_to", []):
                if tgt in self.rnames:
                    d[:, 3 + self.rnames.index(tgt)] += frac * r_3
            for tgt, frac in r.get("k2_to", []):
                if tgt in self.rnames:
                    d[:, 3 + self.rnames.index(tgt)] += frac * r_q
                else:
                    cnt[7] += np.sum(frac * r_q * self.vol)
            if not r.get("k2_to"):
                cnt[7] += np.sum(r_q * self.vol)
            if c.nonlinear:
                r_si = rate_fit(r.get("k_si"), te) * ne * n_i
                r_sup = rate_fit(r.get("k_sup"), te) * ne * n_i
                d[:, 3 + i] -= r_si + r_sup
                d[:, 0] += r_si
                d[:, 1 if r.get("si_ion", "atomic") == "atomic" else 2] += r_si
                cnt[1] += np.sum(r_si * self.vol)
                cnt[8] += np.sum(r_sup * self.vol)
                cnt[9] += np.sum(r_si * self.vol)
                for tgt, k in r.get("e_transfer", []):
                    if tgt in self.rnames:
                        tr = float(k) * ne * n_i
                        d[:, 3 + i] -= tr
                        d[:, 3 + self.rnames.index(tgt)] += tr
        for p in self.pool:
            ia, ib = self.rnames.index(p["a"]), self.rnames.index(p["b"])
            rr = rate_fit(p["k"], te) * nr[:, ia] * nr[:, ib]
            d[:, 3 + ia] -= rr
            d[:, 3 + ib] -= rr
            y_e = float(p.get("e_yield", 1.0))
            d[:, 0] += y_e * rr
            d[:, 1 if p.get("ion", "atomic") == "atomic" else 2] += y_e * rr
            for tgt, frac in p.get("products", []):
                if tgt in self.rnames:
                    d[:, 3 + self.rnames.index(tgt)] += frac * rr
            cnt[2] += np.sum(y_e * rr * self.vol)
            cnt[9] += np.sum(2 * rr * self.vol)
        # wall-photoemission population (0-D, linear chemistry only)
        dw = np.zeros(1 + self.nr)
        nw = max(float(y[self.iw]), 0.0)
        rw = np.maximum(y[self.iw + 1: self.ic], 0.0)
        dw[0] += nw * (N * s["kiz"] + s["nu_hm"])
        cnt[12] += nw * (N * s["kiz"] + s["nu_hm"])
        for i, r in enumerate(self.res):
            dw[1 + i] += nw * s["src"][i] - self._res_linear_loss(r) * rw[i]
            a_s = float(r.get("A_s", 0.0))
            photons += float(r.get("photo", 0.0)) * float(r.get("Y", 0.0)) * a_s * rw[i]
            r_3 = float(r.get("k3", 0.0)) * N * N * rw[i]
            r_q = float(r.get("k2", 0.0)) * N * rw[i]
            for tgt, frac in r.get("k3_to", []):
                if tgt in self.rnames:
                    dw[1 + self.rnames.index(tgt)] += frac * r_3
            for tgt, frac in r.get("k2_to", []):
                if tgt in self.rnames:
                    dw[1 + self.rnames.index(tgt)] += frac * r_q
        pe = c.eta_gs * photons
        dw[0] += pe
        cnt[11] += pe
        return np.concatenate([d.ravel(), dw]), cnt

    def run(self, max_step_frac: float = 0.02):
        c = self.c
        y0 = self.initial()
        ln_br = math.log(c.n_br)

        def ev(t, y):
            return math.log(max(self.totals(y), 1e-300)) - ln_br
        ev.terminal = True
        ev.direction = 1
        n_scr = self.n_screen()

        def ev_scr(t, y):
            return math.log(max(float(y[0: self.m * self.nsp: self.nsp].max()), 1e-300) / n_scr)
        ev_scr.terminal = c.criterion == "screen"
        ev_scr.direction = 1
        atol = np.full_like(y0, 1e-12 * max(float(y0.max()), 1.0))
        atol[self.iw:] = 1e-6
        S = self._sparsity
        groups = group_columns(S)

        def jac(t, y):
            # fixed-step sparse finite differences: scipy's adaptive num_jac overflows its step factor on columns
            # whose derivative is exactly zero (e.g. molecular ions in cells where n_e underflows to 0)
            return approx_derivative(lambda yy: self.rhs(t, yy), y, method="2-point", sparsity=(S, groups))
        sol = solve_ivp(self.rhs, (0.0, c.t_end), y0, method="BDF", jac=jac,
                        rtol=c.rtol, atol=atol, events=[ev, ev_scr], max_step=max(c.t_end * max_step_frac, 1e-12))
        return self._result(sol)

    def n_screen(self, factor: float | None = None) -> float:
        """Electron density at which omega_p^2 / (omega nu_m) = factor (RF depolarization of a compact cloud,
        nu_m >> omega); beyond it the uniform-field radial model is not valid. factor = 1/L: 3 for a sphere,
        2 for a cylinder perpendicular to E (the laser line here is perpendicular to the TE10 field)."""
        c = self.c
        if factor is None:
            factor = 3.0 if c.geometry == "sphere" else 2.0
        e_rms = c.e0_v_m / math.sqrt(2.0)
        nu_m = self.N * float(np.interp(self.sw._x(e_rms / self.N / TD), self.sw.ln_en, self.sw.km))
        w = 2 * math.pi * c.f_hz
        return factor * EPS0 * 9.1093837015e-31 * w * nu_m / QE ** 2

    def _result(self, sol) -> dict:
        m, ns = self.m, self.nsp
        ys = sol.y
        n_cloud = np.array([self.cloud_total(ys[:, k]) for k in range(ys.shape[1])])
        n_wall = np.maximum(ys[self.iw], 0.0)
        hit = sol.t_events[0].size > 0
        yf = ys[:, -1]
        Y = yf[: m * ns].reshape(m, ns)
        rates = np.array([self.counter_rates(sol.t[k], ys[:, k]) for k in range(ys.shape[1])])
        cnt = dict(zip(self.counters, (float(v) for v in _logmean_integral(rates, sol.t))))
        ion_cloud = float(np.sum((Y[:, 1] + Y[:, 2]) * self.vol))
        made = cnt["ion_direct"] + cnt["ion_step"] + cnt["ion_pool"] + cnt["ion_hm"] - cnt["recomb"]
        e_tot = n_cloud[-1] + n_wall[-1] + cnt["outflow_e"]
        sig = self.rms_radius(np.maximum(Y[:, 0], 0))
        s = self.scalars(0.0)
        scr = sol.t_events[1]
        t_scr = float(scr[0]) if scr.size else None
        n_at_scr = float(self.totals(sol.y_events[1][0])) if scr.size else None
        if self.c.criterion == "screen" and scr.size:
            hit = True
        t_bd = (float(sol.t_events[0][0]) if sol.t_events[0].size else t_scr) if hit else None
        return dict(breakdown=bool(hit), t_bd_s=t_bd, status=int(sol.status),
                    t_screen_s=t_scr, n_at_screen=n_at_scr, n_screen_m3=self.n_screen(),
                    message=str(sol.message), t=sol.t, n_cloud=n_cloud, n_wall=n_wall,
                    final_reservoirs={nm: float(np.sum(Y[:, 3 + i] * self.vol)) for i, nm in enumerate(self.rnames)},
                    counters=cnt, rms_radius_m=sig, ne_peak_m3=float(Y[:, 0].max()),
                    electron_balance_rel=(e_tot - self.c.n0 - made - cnt["wall_ion"] - cnt["photo_e"]) / max(e_tot, 1.0),
                    charge_residual_rel=(n_cloud[-1] + cnt["outflow_e"] - ion_cloud) / max(e_tot, 1.0),
                    wall_fraction=float(n_wall[-1] / max(e_tot, 1.0)),
                    x_peak=float(Y[:, 0].max() / self.N),
                    dense_max_x=getattr(self.ee, "max_x", None),
                    dense_used_filled=getattr(self.ee, "used_filled", None),
                    en_eff_Td=s["en"], te_eV=s["te"], nfev=int(sol.nfev), r_max_m=float(self.faces[-1]))


def run_case(case: Case, swarm: Swarm | None = None) -> dict:
    return Avalanche(case, swarm).run()


def _logmean_integral(rates: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Integral of non-negative rates over solver steps with the log-mean rule (exact for exponential growth)."""
    r0, r1 = np.maximum(rates[:-1], 0.0), np.maximum(rates[1:], 0.0)
    dt = np.diff(t)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        lm = (r1 - r0) / np.log(r1 / r0)
    ok = (r0 > 0) & (r1 > 0) & (np.abs(r1 - r0) > 1e-9 * np.maximum(r0, r1))
    seg = np.where(ok, lm, 0.5 * (r0 + r1)) * dt
    return seg.sum(axis=0)


def threshold(case: Case, swarm=None, rel: float = 0.003, hi_v_m: float | None = None) -> float:
    """Lowest peak field for which the deterministic avalanche reaches N_br within the field window."""
    sw = swarm or Swarm(case.gas)
    if hi_v_m is None:
        hi_v_m = direct_threshold_0d(case, sw) * 1.05
    for _ in range(4):
        if run_case(replace(case, e0_v_m=hi_v_m), sw)["breakdown"]:
            break
        hi_v_m *= 1.5
    else:
        return math.inf
    lo_v_m = hi_v_m / 30.0
    while hi_v_m / lo_v_m > 1 + rel:
        mid = math.sqrt(lo_v_m * hi_v_m)
        if run_case(replace(case, e0_v_m=mid), sw)["breakdown"]:
            hi_v_m = mid
        else:
            lo_v_m = mid
    return math.sqrt(lo_v_m * hi_v_m)


def direct_threshold_0d(case: Case, sw: Swarm, dt: float = 0.2e-9) -> float:
    """Direct-ionization envelope threshold (no losses): N0 exp(int N k_iz dt) = N_br (gas_compare definition)."""
    t = np.arange(0.0, case.t_end + dt, dt)
    env = np.array([case.envelope(x) for x in t])
    g_req = math.log(case.n_br / case.n0)

    def gain(e0):
        nu = np.array([case.n_gas * sw.at(sw.effective_en(e0 * v, case.n_gas, case.f_hz) if v > 0 else 1e-3)["kiz"]
                       for v in env])
        return float(np.trapezoid(nu, t))
    lo, hi = 1e5, 1e8
    if gain(hi) < g_req:
        return math.inf
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if gain(mid) >= g_req:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)
