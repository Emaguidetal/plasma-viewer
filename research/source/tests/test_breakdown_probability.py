import math
import unittest

import numpy as np

from rf_plasma.breakdown_probability import Setup, seed_per_length, simulate
from rf_plasma.stochastic_seed import Rates


class BreakdownProbabilityTests(unittest.TestCase):
    def test_seed_scales_as_energy_power_and_saturates(self):
        st = Setup()
        z = np.linspace(-0.03, 0.03, 6001)
        n1 = seed_per_length(st, 1e-3, z).sum()
        n05 = seed_per_length(st, 0.5e-3, z).sum()
        self.assertAlmostEqual(n1 / n05, 2 ** st.k_order, delta=0.01 * 2 ** st.k_order)
        # moving the focus only moves the seed while the beam stays inside the window
        shifted = seed_per_length(Setup(focus_m=5e-3), 1e-3, z).sum()
        self.assertAlmostEqual(shifted / n1, 1.0, delta=1e-3)
        # saturation caps the on-axis density at n_sat
        sat = Setup(n_ref_m3=1e22, n_sat_m3=1e18)
        self.assertLess(seed_per_length(sat, 1e-3, z).max(), 1e18 * math.pi * sat.w0_m ** 2 * 10)

    def test_deterministic_delay_matches_growth_time(self):
        r = Rates("He")
        st = Setup(e_peak_v_m=150e5, n_ref_m3=1e22, sig_p=0.0, sig_el=0.0, n_br=1e14)
        d, n0 = simulate(st, 2, rates=r)
        nu = float(r.nu_iz(st.e_peak_v_m, st.n_gas, st.f_hz))
        np.testing.assert_allclose(d, math.log(st.n_br / n0) / nu, rtol=0.02)

    def test_ring_down_limits_breakdown_near_threshold(self):
        r = Rates("He")
        kw = dict(e_peak_v_m=65e5, sig_p=0.0, sig_el=0.0, n_ref_m3=1e20)
        long_p = np.isfinite(simulate(Setup(tau_amp_s=1e-6, **kw), 3, rates=r)[0]).mean()
        short_p = np.isfinite(simulate(Setup(tau_amp_s=1e-9, **kw), 3, rates=r)[0]).mean()
        self.assertEqual((long_p, short_p), (1.0, 0.0))


if __name__ == "__main__":
    unittest.main()
