import unittest

import numpy as np

from rf_plasma.helium_cr import (
    E_P, E_S, STATE_ENERGIES, excitation_rates, ionization_rates,
    local_collisions, radiation_rhs, resonance_transport,
)
from rf_plasma.sheet_cavity import C0, laser_front_fraction, laser_intensity, sheet_fraction


class HeliumCRTests(unittest.TestCase):
    def test_sequential_ionization_conserves_charge_nuclei_and_energy(self):
        population = np.array([0.75, 0.02, 0.03, 0.18, 0.02])
        ne = population[3] + 2 * population[4]
        dp, du = local_collisions(population, 1.5 * ne * 80, 1e20)
        self.assertGreater(dp[4], 0)
        self.assertAlmostEqual(float(dp.sum()), 0, delta=1e-8 * np.abs(dp).max())
        self.assertAlmostEqual(float(STATE_ENERGIES @ dp + du), 0, delta=1e-8 * abs(du))
        with self.assertRaises(ValueError):
            ionization_rates(2.0)

    def test_superelastic_pair_has_boltzmann_equilibrium(self):
        te = np.array([3.0, 10.0, 50.0])
        up, down = excitation_rates(te)
        expected = np.array([1.0, 3.0])[:, None] * np.exp(-np.array([E_S, E_P])[:, None] / te)
        np.testing.assert_allclose(up / down, expected, rtol=1e-13)
        # A cold electron population receives energy from overpopulated singlets.
        pop = np.array([0.0, 0.5, 0.4, 0.1, 0.0])
        dp, du = local_collisions(pop, 1.5 * 0.1 * 3, 1e20)
        self.assertGreater(du, 0)
        self.assertLess(dp[1], 0)
        self.assertLess(dp[2], 0)

    def test_radiation_transport_is_nonlocal_and_conservative(self):
        off, escape, outgoing = resonance_transport(7, 0.002, 1e20)
        self.assertTrue(np.all(off >= 0))
        self.assertTrue(np.all(escape >= 0))
        self.assertTrue(np.all(outgoing <= 1 + 1e-10))
        self.assertGreater(off[3, 2], 0)
        excited = np.array([0, 0, 0, 1e15, 0, 0, 0], dtype=float)
        dp, lost = radiation_rhs(excited, off, escape, outgoing)
        self.assertGreater(dp[2], 0)
        self.assertAlmostEqual(float(E_P * dp.sum() + lost.sum()), 0,
                               delta=1e-10 * float(lost.sum()))
        thin = resonance_transport(3, 0.002, 0)
        np.testing.assert_allclose(thin[0], 0, atol=1e-15)
        np.testing.assert_allclose(thin[1], 1, atol=1e-13)
        thick = resonance_transport(7, 0.002, 1e24)
        self.assertLess(thick[1][3], escape[3])

    def test_prescribed_trigger_is_not_active_before_its_time(self):
        f = sheet_fraction(np.array([0, 1, 1.5, 2, 4]), 1, 1)
        np.testing.assert_allclose(f, [0, 0, 0.5, 1, 1], atol=1e-15)

    def test_laser_front_ionizes_only_where_the_pulse_has_arrived(self):
        x = np.array([0.0, 0.036, 0.072])
        t = 1e-9 + 0.036 / C0
        f = laser_front_fraction(t, x, 1e-9, 0.6e-9)
        self.assertGreater(f[0], 0)
        self.assertEqual(f[1], 0)
        self.assertEqual(f[2], 0)
        self.assertAlmostEqual(float(laser_intensity(t, 0.036, 1e-9, 30e-12)), 1.0)


if __name__ == "__main__":
    unittest.main()
