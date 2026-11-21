import unittest

import numpy as np

from rf_plasma.kinetics import ReactionRates, single_ion_benchmark


class SingleIonBenchmarkTests(unittest.TestCase):
    def test_recombination_matches_analytic_decay_and_conserves_particles(self):
        total, initial, alpha = 1e20, 1e18, 2e-13
        time = np.linspace(0, 20e-6, 51)
        result = single_ion_benchmark(total, initial, ReactionRates(0, 0, alpha), time)
        expected = initial / (1 + alpha * initial * time)
        np.testing.assert_allclose(result["electron_density_m3"], expected, rtol=1e-7)
        np.testing.assert_allclose(result["ion_density_m3"], result["electron_density_m3"])
        np.testing.assert_allclose(result["ion_density_m3"] + result["neutral_density_m3"], total, rtol=1e-14)

    def test_impact_ionization_matches_logistic_limit(self):
        total, initial, rate = 1e20, 1e17, 1e-14
        time = np.linspace(0, 4e-6, 41)
        result = single_ion_benchmark(total, initial, ReactionRates(0, rate, 0), time)
        expected = total / (1 + (total / initial - 1) * np.exp(-rate * total * time))
        np.testing.assert_allclose(result["electron_density_m3"], expected, rtol=1e-7)

    def test_constant_photoionization_matches_neutral_depletion(self):
        total, rate = 1e20, 1e6
        time = np.linspace(0, 4e-6, 41)
        result = single_ion_benchmark(total, 0, ReactionRates(rate, 0, 0), time)
        expected = total * (1 - np.exp(-rate * time))
        np.testing.assert_allclose(result["electron_density_m3"], expected, rtol=1e-7, atol=1e6)

    def test_invalid_rates_and_seeds_are_rejected(self):
        with self.assertRaises(ValueError):
            ReactionRates(0, -1, 0)
        with self.assertRaises(ValueError):
            single_ion_benchmark(1e20, 2e20, ReactionRates(0, 0, 0), np.array([0, 1]))


if __name__ == "__main__":
    unittest.main()
