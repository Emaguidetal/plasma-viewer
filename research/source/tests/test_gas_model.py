import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rf_plasma import plasma3d as k3  # noqa: E402
from rf_plasma.gas_model import load_gas_model  # noqa: E402
from rf_plasma.swarm import SWARM_DIR, read_table  # noqa: E402

N_GAS = 1.909e26  # 8 bar, 300 K
TD = 1e-21


def _dc_equilibrium(model, en_td, chunks=400, dt=2e-10):
    """Drive one cell with ohmic DC heating sigma(eps) E^2 until eps settles."""
    prm = model.prm.copy()
    for p in (k3.P_B1_0, k3.P_B1N_0, k3.P_B2_0, k3.P_B2N_0, k3.P_K3B0, k3.P_KCONV0, k3.P_KDISS0, k3.P_KDEX0):
        prm[p] = 0.0  # isolate electron energetics
    e_field = en_td * TD * N_GAS
    st = np.array([1e12, 0.0, 1e12 * 0.1])
    acc = np.zeros(9)
    for _ in range(chunks):
        ne = st[0] + st[1]
        eps = st[2] / ne
        nu = N_GAS * k3._lookup(model.tabs_bulk[0], model.ln_e0, model.h, max(eps, 1e-3))
        p_ext = ne * k3.QE * e_field**2 / (k3.ME * nu)  # eV m^-3 s^-1
        k3.kinetic_cell(st, dt, N_GAS, 300.0, p_ext, 0.0, prm, model.tabs_bulk, model.ln_e0, model.h,
                        model.ln_i0, model.hi, model.ppt_tab, acc)
    return st[2] / (st[0] + st[1])


class GasModelTests(unittest.TestCase):
    def test_tables_finite_nonnegative(self):
        for g in ("Ar", "He", "N2"):
            m = load_gas_model(g, 532e-9, allow_unverified=True)
            for t in (m.tabs_bulk, m.tabs_filament):
                self.assertTrue(np.isfinite(t).all() and (t >= 0).all(), g)

    def test_unverified_rejected_by_default(self):
        with self.assertRaises(ValueError):
            load_gas_model("N2", 532e-9)

    def test_dc_local_field_reproduces_boltzmann_mean_energy(self):
        """kinetic_cell + tables must recover the BOLOS eps(E/N) (local-field consistency)."""
        cases = {"Ar": 10.0, "He": 5.0, "N2": 60.0}
        for g, en in cases.items():
            m = load_gas_model(g, 532e-9, allow_unverified=True)
            ref = read_table(SWARM_DIR / f"{g}_boltzmann.csv")
            eps_ref = float(np.exp(np.interp(math.log(en), np.log(ref["E_N_Td"]), np.log(ref["mean_energy_eV"]))))
            eps = _dc_equilibrium(m, en)
            self.assertAlmostEqual(eps / eps_ref, 1.0, delta=0.03, msg=f"{g}: {eps:.3f} vs {eps_ref:.3f} eV")


if __name__ == "__main__":
    unittest.main()
