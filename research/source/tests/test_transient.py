from pathlib import Path
import unittest

import numpy as np

from rf_plasma.config import load_reference_case
from rf_plasma.drude import peak_phasor_heating_w_m3
from rf_plasma.transient import density_and_growth, simulate

CASE_FILE = Path(__file__).resolve().parents[1] / "cases" / "reference.json"


class PrescribedTransientTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        case = load_reference_case(CASE_FILE)
        cls.control = simulate(case, 300e-12, 1e15, 1e15, 2e11, 7.5e-12)
        cls.coarse = simulate(case, 300e-12, 1e22, 1e15, 2e11, 7.5e-12)
        cls.fine = simulate(case, 300e-12, 1e22, 1e15, 2e11, 3.75e-12)
        cls.case = case

    def test_prescribed_birth_integrates_to_density_change(self):
        time = np.linspace(-2e-9, 2e-9, 5001)
        density, growth = density_and_growth(time, 1e15, 1e22, 300e-12)
        self.assertAlmostEqual(float(np.trapezoid(growth, time) / (density[-1] - density[0])), 1, delta=1e-5)
        self.assertTrue(np.all(growth >= 0))
        with self.assertRaises(ValueError):
            density_and_growth(0, 1e16, 1e15, 300e-12)

    def test_seed_only_recovers_steady_drude_heating(self):
        control, case = self.control, self.case
        t = control["time_s"]
        end = float(t[-1])
        start = end - 5 / case.rf.frequency_hz
        sample = np.linspace(start, end, 1001)
        mean_power = np.trapezoid(np.interp(sample, t, control["collisional_power_w"]), sample) / (end - start)
        expected = (
            peak_phasor_heating_w_m3(1e15, 2e11, case.rf.frequency_hz, case.rf.empty_focus_field_peak_v_m)
            * control["metrics"]["plasma_volume_m3_illustrative"]
        )
        self.assertAlmostEqual(float(mean_power / expected), 1, delta=0.005)
        self.assertEqual(control["metrics"]["final_birth_mixing_energy_j"], 0)
        self.assertLess(control["metrics"]["max_energy_residual_relative"], 1e-6)

    def test_dynamic_energy_balance_and_timestep_convergence(self):
        coarse, fine = self.coarse["metrics"], self.fine["metrics"]
        self.assertLess(coarse["max_energy_residual_relative"], 1e-6)
        self.assertLess(fine["max_energy_residual_relative"], 1e-6)
        self.assertGreater(fine["final_collisional_energy_j"], 0)
        self.assertGreater(fine["final_birth_mixing_energy_j"], 0)
        relative_change = abs(
            coarse["final_collisional_energy_j"] / fine["final_collisional_energy_j"] - 1
        )
        self.assertLess(relative_change, 0.05)


if __name__ == "__main__":
    unittest.main()
