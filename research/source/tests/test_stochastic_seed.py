import math
import unittest

import numpy as np

from rf_plasma.stochastic_seed import Rates, Shot, simulate


class StochasticSeedTests(unittest.TestCase):
    def test_single_electron_yule_delay_is_gumbel(self):
        """One seed electron, no loss, constant field: t = (ln N_br - ln W)/nu with W ~ Exp(1)."""
        r = Rates("He")
        sh = Shot(laser=True, n_laser=1.0, n_br=1e10, t_laser_s=0.0, t_rf_s=0.5e-6, dt_s=0.01e-9)
        nu = float(r.nu_iz(sh.e0_v_m, sh.n_gas, sh.f_hz))
        tb, _ = simulate(sh, 6000, seed=3, rates=r, sources=("yule",))
        self.assertTrue(np.isfinite(tb).all())
        np.testing.assert_allclose(tb.mean(), (math.log(1e10) + 0.5772157) / nu, rtol=0.01)
        np.testing.assert_allclose(tb.std(), math.pi / math.sqrt(6) / nu, rtol=0.05)

    def test_deterministic_limit_and_loss_survival(self):
        r = Rates("He")
        sh = Shot(n_laser=1e6, n_br=1e12, t_laser_s=0.0, t_rf_s=0.5e-6, dt_s=0.01e-9)
        nu = float(r.nu_iz(sh.e0_v_m, sh.n_gas, sh.f_hz))
        tb, _ = simulate(sh, 3, rates=r, sources=())
        np.testing.assert_allclose(tb, math.log(1e6) / nu, rtol=1e-3)
        # with loss nu/2 only ~half of single-electron lineages survive
        sh1 = Shot(n_laser=1.0, n_br=1e8, t_laser_s=0.0, t_rf_s=1e-6, nu_loss_s=0.5 * nu)
        tb1, _ = simulate(sh1, 3000, seed=5, rates=r, sources=("yule",))
        self.assertAlmostEqual(np.isfinite(tb1).mean(), 0.5, delta=0.03)


if __name__ == "__main__":
    unittest.main()
