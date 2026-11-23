"""Focused tests for rf_plasma.excited_avalanche (radial avalanche with excited-state reservoirs)."""
import math
import unittest
from dataclasses import replace

import numpy as np

from rf_plasma import excited_avalanche as xa
from rf_plasma import gas_compare as gc
from rf_plasma import swarm_excited as sx


class TestTables(unittest.TestCase):
    def test_group_counts_sum_to_count_rate_and_ai_subset(self):
        for gas in ("He", "Ar", "N2"):
            d = sx.load(gas)
            sel = np.asarray(d["E_N_Td"]) > 5.0          # below ~1 Td all rates are round-off (~1e-32)
            tot = sum(np.asarray(d[g]) for g in sx.GROUP_COLUMNS[gas])[sel]
            self.assertTrue(np.all(np.asarray(d["k_ai"])[sel] <= tot * (1 + 1e-9)))
        ar = sx.load("Ar")
        i = int(np.argmin(abs(np.log(np.asarray(ar["E_N_Td"]) / 20.0))))
        ratio = sum(ar[g][i] for g in sx.GROUP_COLUMNS["Ar"]) / ar["kiz_m3_s"][i]
        self.assertGreater(ratio, 500)          # Ar excitation >> ionization near 20 Td


class TestSolver(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sw = {g: xa.Swarm(g) for g in ("He", "Ar")}

    def test_direct_limit_matches_gas_compare(self):
        # gas_compare has no recombination; with alpha_dr = 0 the radial direct limit must reproduce it
        for gas in ("He", "Ar", "N2"):
            kin = xa.load_kinetics()[gas]
            kin["ion"]["alpha_dr"] = None
            case = xa.Case(gas=gas, psig=200.0, excited=False, kinetics=kin)
            sw = xa.Swarm(gas)
            e0d = xa.direct_threshold_0d(case, sw)
            old = gc.threshold_envelope(gc.rates(gas), gc.setup(200.0), math.log(1e12 / 376.0))
            self.assertAlmostEqual(e0d / old, 1.0, delta=0.01)
            self.assertTrue(xa.run_case(replace(case, e0_v_m=1.02 * e0d), sw)["breakdown"])
            self.assertFalse(xa.run_case(replace(case, e0_v_m=0.98 * e0d), sw)["breakdown"])

    def test_free_diffusion_of_single_electron(self):
        # weak field: negligible ionization; one electron: lambda_D >> cloud -> free diffusion sigma^2 = s0^2 + 2Dt
        case = xa.Case(gas="He", psig=50.0, e0_v_m=1e3, excited=False, geometry="sphere", n0=1.0,
                       sigma0_m=5e-6, t_flat_s=20e-9, n_tau_after=0.0, cells=200)
        av = xa.Avalanche(case, self.sw["He"])
        r = av.run()
        s = av.scalars(0.0)
        expect = math.sqrt(case.sigma0_m ** 2 + 2 * s["de"] * case.t_flat_s)
        self.assertAlmostEqual(r["rms_radius_m"] / expect, 1.0, delta=0.02)

    def test_balances_and_channels(self):
        case = xa.Case(gas="Ar", psig=200.0)
        e0 = xa.direct_threshold_0d(case, self.sw["Ar"])
        lin = xa.run_case(replace(case, e0_v_m=e0, nonlinear=False), self.sw["Ar"])
        full = xa.run_case(replace(case, e0_v_m=e0), self.sw["Ar"])
        fb = xa.run_case(replace(case, e0_v_m=e0, eta_gs=1e-2), self.sw["Ar"])
        for r in (lin, full, fb):
            self.assertLess(abs(r["electron_balance_rel"]), 1e-2)
            self.assertLess(abs(r["charge_residual_rel"]), 1e-2)
        self.assertEqual(lin["counters"]["ion_step"], 0.0)
        self.assertEqual(lin["wall_fraction"], 0.0)
        self.assertGreater(full["counters"]["ion_step"] + full["counters"]["ion_pool"], 0.0)
        self.assertGreater(fb["counters"]["photo_e"], 0.0)
        # extra electron sources can only advance breakdown
        t = lambda r: r["t_bd_s"] if r["t_bd_s"] is not None else math.inf
        self.assertLessEqual(t(full), t(lin) * 1.001)
        self.assertLessEqual(t(fb), t(full) * 1.001)

    def test_uncertainty_sets_direction(self):
        kin = xa.load_kinetics()["Ar"]
        lo, hi = xa.scaled_kinetics(kin, "low"), xa.scaled_kinetics(kin, "high")
        r0 = kin["reservoirs"][0]
        self.assertLess(lo["reservoirs"][0]["k_si"]["A"], r0["k_si"]["A"])
        self.assertGreater(hi["reservoirs"][0]["k_si"]["A"], r0["k_si"]["A"])
        self.assertGreater(lo["reservoirs"][0]["k2"], r0["k2"])     # losses up in the low-ionization set
        self.assertGreater(hi["pooling"][0]["k"], kin["pooling"][0]["k"])


class TestDenseRegime(unittest.TestCase):
    def test_dense_table_dilute_limit_and_monotone_fill(self):
        for gas in ("He", "Ar"):
            dt = xa.DenseTable(gas)
            sw = xa.Swarm(gas)
            for en in (10.0, 20.0, 50.0):
                k0 = float(dt.kiz(en, np.array([0.0]))[0])
                self.assertAlmostEqual(k0 / sw.at(en)["kiz"], 1.0, delta=0.05)
            k = dt.kiz(6.0, np.array([0.0, 1e-6, 1e-5, 1e-4, 1e-3]))
            self.assertTrue(np.all(np.isfinite(k)) and np.all(k > 0))
            self.assertTrue(dt.filled[1:].any())            # non-converged points exist and are tracked
            self.assertTrue(np.all(dt.te(20.0, np.array([0.0, 1e-3])) > 0.5))

    def test_dense_seed_run_conserves(self):
        case = xa.Case(gas="Ar", psig=200.0, e0_v_m=80e5, n0=1.1e10, sigma0_m=5e-6, ee="dense", n_br=1e13)
        r = xa.run_case(case)
        self.assertTrue(r["breakdown"])
        self.assertLess(abs(r["electron_balance_rel"]), 1e-2)


if __name__ == "__main__":
    unittest.main()
