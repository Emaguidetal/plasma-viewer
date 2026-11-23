"""Two-term Boltzmann rate tables including electron-electron collisions (ionization-degree tables).

Replaces the interpolation k = (1-w) k_Boltzmann + w k_Maxwell used by the base kernels by a solved
EEDF. The e-e operator is the two-term form of Hagelaar & Pitchford, Plasma Sources Sci. Technol. 14,
722 (2005), Sec. 2.6 (eqs. 29-32, BOLSIG+): with a = e^2 gamma / (24 pi eps0^2) (ne/N) lnL,
  W_ee = 3 a A1,  D_ee = 2 a (A2 + eps^(3/2) A3),
  A1 = int_0^eps u^(1/2) F0 du,  A2 = int_0^eps u^(3/2) F0 du,  A3 = int_eps^inf F0 du,
added to the drift (-W) and diffusion (D) coefficients of the energy flux; the problem becomes
nonlinear and is solved by fixed-point iteration of BOLOS's linear solver.
Tables: for each ionization degree x = ne/N on a log grid, a sweep in E/N gives rows (km, kiz, kexc,
kel) vs mean energy, same definitions as swarm.py. The Coulomb logarithm uses the NRL form at the
gas density given (weak dependence; a table per pressure is written).
"""

from pathlib import Path
import math

import numpy as np
import scipy.constants as co
from scipy import integrate

from bolos import solver as bsolver

from .swarm import SWARM_DIR, _row, make_solver

GAMMA = math.sqrt(2 * co.e / co.m_e)
X_LEVELS = np.array([0.0, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2])


def coulomb_log(ne, te_ev):
    te = max(te_ev, 0.05)
    return max(23.0 - math.log(math.sqrt(ne * 1e-6) * te**-1.5), 2.0)


class EESolver:
    """Wraps a BOLOS solver; adds the e-e drift/diffusion terms on the cell boundaries.
    ei=True also adds electron-ion momentum transfer (Coulomb, sigma_ei = pi (e^2/(4 pi eps0 eps))^2 lnL) as an
    elastic target 'ION' of density x * lnL (cross section tabulated for lnL = 1), mass ratio m_e/M_ion."""

    def __init__(self, gas, n_gas, e_max=120.0, cells=400, tg_k=300.0, ei=False):
        self.b = make_solver(gas, e_max=e_max, cells=cells, tg_k=tg_k)
        self.n_gas = n_gas
        self.x = 0.0
        self.ei = ei
        b = self.b
        if ei:
            eps = np.geomspace(1e-4, e_max * 1.2, 400)
            k_c = (co.e / (4 * math.pi * co.epsilon_0)) ** 2      # (e^2/(4 pi eps0))^2 / e^2 -> eV^2 m^2 below
            sig = math.pi * k_c / eps ** 2                          # eps in eV: (e/(4 pi eps0 eps[V]))^2 pi
            mass = {"He": 4.0026, "Ar": 39.948, "N2": 2 * 28.014}.get(gas, 40.0) * 1.66054e-27
            b.load_collisions([dict(kind="ELASTIC", target="ION", mass_ratio=co.m_e / mass,
                                    data=np.c_[np.r_[0.0, eps], np.r_[sig[0], sig]].tolist(),
                                    comment="Coulomb e-ion momentum transfer per unit lnL")])
            b.target["ION"].density = 0.0
            b.init()
        self._W0 = None
        orig = b._scharf_gummel

        def sg(sigma_tilde, G=0, _self=self):
            if _self.x <= 0.0 or _self._F is None:
                return orig(sigma_tilde, G)
            W, DB = b.W.copy(), b.DB.copy()
            wee, dee = _self._ee_terms(_self._F)
            b.W = W - wee
            b.DB = DB + dee
            try:
                return orig(sigma_tilde, G)
            finally:
                b.W, b.DB = W, DB

        b._scharf_gummel = sg
        orig_it = b.iterate

        def it(f0, delta=1e14, _self=self):
            # e-e coefficients lagged by one inner iteration: BOLOS's adaptive converge then
            # solves the nonlinear problem directly
            if _self.x > 0.0:
                _self._F = f0
            return orig_it(f0, delta=delta)

        b.iterate = it
        self._F = None

    def _ee_terms(self, F):
        b = self.b
        ce, be = b.cenergy, b.benergy
        de = np.diff(be)
        c1 = np.r_[0.0, np.cumsum(np.sqrt(ce) * F * de)]
        c2 = np.r_[0.0, np.cumsum(ce**1.5 * F * de)]
        c3 = np.r_[np.cumsum((F * de)[::-1])[::-1], 0.0]
        eps_m = float(b.mean_energy(F))
        lnl = coulomb_log(self.x * self.n_gas, 2.0 * eps_m / 3.0)
        a = co.e**2 * GAMMA / (24 * math.pi * co.epsilon_0**2) * self.x * lnl
        return 3 * a * c1, 2 * a * (c2 + be**1.5 * c3)

    def solve(self, en_td, x, f0, tol=1e-7):
        b = self.b
        if self.ei:
            b.target["ION"].density = 0.0
        b.EN = en_td * bsolver.TOWNSEND
        b.init()
        self.x = 0.0
        self._F = None
        f = b.converge(f0, maxn=400, rtol=1e-6)
        if x <= 0.0:
            return f
        self.x = x
        try:
            if self.ei:
                for _ in range(3):           # lnL from the current EEDF (weak dependence; 3 passes)
                    te = 2.0 * float(b.mean_energy(f)) / 3.0
                    b.target["ION"].density = x * coulomb_log(x * self.n_gas, te)
                    b.init()
                    f = b.converge(f, maxn=3000, rtol=tol, delta0=1e12)
            else:
                f = b.converge(f, maxn=3000, rtol=tol, delta0=1e12)
        finally:
            self._F = None
        return f

    def ee_power_residual(self, F):
        """Net e-e energy flux into the grid (should vanish: e-e conserves energy)."""
        if self.x <= 0.0:
            return 0.0
        wee, dee = self._ee_terms(F)
        b = self.b
        fb = np.interp(b.benergy, b.cenergy, F)
        dfb = np.gradient(fb, b.benergy)
        flux = wee * fb - dee * dfb  # energy-space flux due to e-e (sign: + towards higher energy)
        return float(integrate.simpson(flux, x=b.benergy))


def ee_table(gas, n_gas, x, en_td=None):
    en_td = np.geomspace(0.03, 3000.0, 90) if en_td is None else en_td
    s = EESolver(gas, n_gas)
    f0 = s.b.maxwell(0.5)
    rows = []
    for en in en_td:
        try:
            f1 = s.solve(float(en), x, f0)
        except bsolver.ConvergenceError:
            continue
        s.b.EN = en * bsolver.TOWNSEND
        s.b.init()
        rows.append(_row(s.b, f1, float(en)))
        f0 = f1
    return np.array(rows)


def table_path(gas, n_gas, x):
    return SWARM_DIR / "ee" / f"{gas}_N{n_gas:.3e}_x{x:.0e}.csv"


def nearest_complete(gas, n_gas, max_ratio=10.0, levels=X_LEVELS):
    """N of an already-built complete table set closest (in log N) to n_gas, else n_gas.
    N enters only through ln(Lambda) at fixed ne/N, so a factor-of-few mismatch shifts ln(Lambda) by < 10 %."""
    best, best_d = n_gas, math.inf
    for p in (SWARM_DIR / "ee").glob(f"{gas}_N*_x0e+00.csv"):
        n = float(p.name.split("_N")[1].split("_x")[0])
        d = abs(math.log(n / n_gas))
        if d < best_d and d <= math.log(max_ratio) and all(table_path(gas, n, x).exists() for x in levels):
            best, best_d = n, d
    return best


def build(gas, n_gas, levels=X_LEVELS):
    from .swarm import write_table
    out = []
    for x in levels:
        p = table_path(gas, n_gas, x)
        if not p.exists():
            write_table(p, ee_table(gas, n_gas, float(x)), gas, f"boltzmann + e-e, x = ne/N = {x:g}, N = {n_gas:.4e}")
        out.append(p)
    return out


if __name__ == "__main__":
    import sys
    gas = sys.argv[1] if len(sys.argv) > 1 else "He"
    n = float(sys.argv[2]) if len(sys.argv) > 2 else 2.414e26
    for p in build(gas, n):
        print(p)
