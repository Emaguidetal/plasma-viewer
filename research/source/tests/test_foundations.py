import json
import math
from pathlib import Path
import unittest

import numpy as np
from scipy.constants import c, epsilon_0, mu_0

from rf_plasma.circular_rf import tm01_impedance_ohm
from rf_plasma.cold_rf import GIB, _resonance, estimated_memory_bytes, memory_budget_bytes
from rf_plasma.config import Gas, gauge_psi_to_absolute_pa, load_reference_case
from rf_plasma.drude import conductivity_s_m, peak_phasor_heating_w_m3, relative_permittivity, te10_slab_scattering
from rf_plasma.laser_study import calculate as calculate_laser_study
from rf_plasma.optics import (
    beam_radius_m,
    energy_before_window_j,
    gaussian_intensity_w_m2,
    ideal_focal_length_m,
    rayleigh_length_m,
)
from rf_plasma.transition_rf import transition_polygon
from rf_plasma.waveguides import (
    circular_cutoff_hz,
    isolated_te101_copper_q,
    rectangular_cutoff_hz,
    skin_depth_m,
    surface_resistance_ohm,
    te10_copper_power_attenuation_per_m,
    te10_power_for_peak_electric_field_w,
    te10_propagation_constant_per_m,
    te101_frequency_hz,
    te101_length_for_frequency_m,
)

CASE_FILE = Path(__file__).resolve().parents[1] / "cases" / "reference.json"


class ReferenceCaseTests(unittest.TestCase):
    def setUp(self):
        self.case = load_reference_case(CASE_FILE)

    def test_pressure_and_reduced_field_use_absolute_pressure(self):
        expected = gauge_psi_to_absolute_pa(100, 101325)
        self.assertAlmostEqual(self.case.gas.absolute_pressure_pa, expected, delta=0.01)
        self.assertAlmostEqual(self.case.gas.ideal_neutral_density_m3 / 1e26, 1.9092, delta=1e-4)
        self.assertAlmostEqual(
            self.case.rf.empty_focus_field_peak_v_m / self.case.gas.ideal_neutral_density_m3 / 1e-21,
            52.38,
            delta=0.05,
        )

    def test_invalid_gas_and_pressure_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "species"):
            Gas("air", 1, 101325, 300)
        with self.assertRaisesRegex(ValueError, "absolute pressure"):
            Gas("He", -200000, 101325, 300)

    def test_case_is_explicitly_illustrative(self):
        self.assertEqual(self.case.status, "illustrative")
        with CASE_FILE.open(encoding="utf-8") as handle:
            raw = json.load(handle)
        self.assertEqual(raw["gas"]["species"], "Ar")


class GaussianBeamTests(unittest.TestCase):
    def setUp(self):
        self.laser = load_reference_case(CASE_FILE).laser

    def test_gaussian_integrates_to_focal_energy(self):
        times = np.linspace(-3 * self.laser.intensity_fwhm_s, 3 * self.laser.intensity_fwhm_s, 1001)
        radii = np.linspace(0, 3 * self.laser.waist_radius_m, 1001)
        intensity = gaussian_intensity_w_m2(self.laser, radii[:, None], 0, times[None, :])
        spatial_power = np.trapezoid(2 * math.pi * radii[:, None] * intensity, radii, axis=0)
        energy = np.trapezoid(spatial_power, times)
        self.assertLess(abs(energy / self.laser.energy_at_focus_j - 1), 0.005)

    def test_fwhm_retardation_and_focus_geometry(self):
        peak = gaussian_intensity_w_m2(self.laser, 0, 0, 0)
        self.assertAlmostEqual(
            float(gaussian_intensity_w_m2(self.laser, 0, 0, self.laser.intensity_fwhm_s / 2) / peak),
            0.5,
        )
        zr = rayleigh_length_m(self.laser)
        self.assertAlmostEqual(zr * 1000, 2.362, delta=0.002)
        self.assertAlmostEqual(float(beam_radius_m(self.laser, zr) / self.laser.waist_radius_m), math.sqrt(2))
        delayed_peak = gaussian_intensity_w_m2(self.laser, 0, zr, zr / c)
        self.assertAlmostEqual(float(delayed_peak / peak), 0.5)
        self.assertAlmostEqual(float(peak / 1e4 / 1e11), 4.984, delta=0.005)
        self.assertAlmostEqual(ideal_focal_length_m(self.laser, 1e-3), 0.1181, delta=0.0002)

    def test_window_loss_is_not_silently_assumed(self):
        self.assertAlmostEqual(energy_before_window_j(self.laser, 0.8), 0.00125)
        with self.assertRaises(ValueError):
            energy_before_window_j(self.laser, 0)

    def test_sensitivity_scan_covers_pulse_and_focus_endpoints(self):
        records = calculate_laser_study(load_reference_case(CASE_FILE))
        self.assertEqual(len(records), 6)
        self.assertEqual({row["waist_radius_m"] for row in records}, {10e-6, 20e-6, 50e-6})
        self.assertEqual({row["intensity_fwhm_s"] for row in records}, {300e-12, 3e-9})
        peaks = {row["intensity_fwhm_s"]: row["peak_intensity_w_m2"]
                 for row in records if row["waist_radius_m"] == 20e-6}
        self.assertAlmostEqual(peaks[300e-12] / peaks[3e-9], 10)


class ElectromagneticBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.case = load_reference_case(CASE_FILE)

    def test_guide_and_closed_cavity_frequencies(self):
        geom = self.case.rectangular
        self.assertAlmostEqual(rectangular_cutoff_hz(geom.width_m, geom.height_m, 1, 0) / 1e9, 2.078, delta=0.002)
        self.assertAlmostEqual(te101_frequency_hz(geom) / 1e9, 3, delta=0.001)
        self.assertAlmostEqual(te101_length_for_frequency_m(geom.width_m, 3e9), geom.closed_cavity_length_m, delta=1e-5)
        radius = self.case.circular.diameter_m / 2
        self.assertAlmostEqual(circular_cutoff_hz(radius, "TE", 1, 1) / 1e9, 1.952, delta=0.002)
        self.assertAlmostEqual(circular_cutoff_hz(radius, "TM", 0, 1) / 1e9, 2.550, delta=0.002)
        with self.assertRaises(ValueError):
            te101_length_for_frequency_m(geom.width_m, 1e9)

    def test_copper_skin_depth_and_q_are_positive(self):
        copper = self.case.copper
        depth = skin_depth_m(3e9, copper)
        self.assertAlmostEqual(depth / 1e-6, 1.21, delta=0.02)
        self.assertGreater(isolated_te101_copper_q(self.case.rectangular, copper), 0)

    def test_copper_q_matches_independent_field_integrals(self):
        geom, copper = self.case.rectangular, self.case.copper
        a, b, d = geom.width_m, geom.height_m, geom.closed_cavity_length_m
        omega = 2 * math.pi * te101_frequency_hz(geom)
        x, z = np.linspace(0, a, 101), np.linspace(0, d, 101)
        xx, zz = np.meshgrid(x, z, indexing="ij")
        ey = np.sin(math.pi * xx / a) * np.sin(math.pi * zz / d)
        hx = math.pi / (d * omega * mu_0) * np.sin(math.pi * xx / a) * np.cos(math.pi * zz / d)
        hz = math.pi / (a * omega * mu_0) * np.cos(math.pi * xx / a) * np.sin(math.pi * zz / d)
        integral = lambda field: np.trapezoid(np.trapezoid(field, z, axis=1), x)
        energy = epsilon_0 * b * integral(ey**2) / 2
        wall_h2 = (
            b * (np.trapezoid(hz[0] ** 2, z) + np.trapezoid(hz[-1] ** 2, z))
            + b * (np.trapezoid(hx[:, 0] ** 2, x) + np.trapezoid(hx[:, -1] ** 2, x))
            + 2 * integral(hx**2 + hz**2)
        )
        wall_power = surface_resistance_ohm(omega / (2 * math.pi), copper) * wall_h2 / 2
        self.assertAlmostEqual(isolated_te101_copper_q(geom, copper), omega * energy / wall_power, delta=0.1)

    def test_te10_copper_attenuation_matches_wall_power_integral(self):
        geom, copper, frequency = self.case.rectangular, self.case.copper, 3e9
        a, b = geom.width_m, geom.height_m
        omega = 2 * math.pi * frequency
        beta = te10_propagation_constant_per_m(frequency, a).real
        x = np.linspace(0, a, 501)
        hx = beta / (omega * mu_0) * np.sin(math.pi * x / a)
        hz = math.pi / (a * omega * mu_0) * np.cos(math.pi * x / a)
        surface_h2_per_m = 2 * b * hz[0]**2 + 2 * np.trapezoid(hx**2 + hz**2, x)
        loss_per_m = surface_resistance_ohm(frequency, copper) * surface_h2_per_m / 2
        power = te10_power_for_peak_electric_field_w(frequency, a, b, 1)
        self.assertAlmostEqual(te10_copper_power_attenuation_per_m(frequency, geom, copper), loss_per_m / power)

    def test_guide_peak_field_power_and_cutoff(self):
        geom = self.case.rectangular
        f, e = 3e9, 1e5
        beta = te10_propagation_constant_per_m(f, geom.width_m).real
        impedance = 2 * math.pi * f * mu_0 / beta
        expected = geom.width_m * geom.height_m / 4 * e**2 / impedance
        self.assertAlmostEqual(
            te10_power_for_peak_electric_field_w(f, geom.width_m, geom.height_m, e), expected
        )
        self.assertEqual(te10_power_for_peak_electric_field_w(f, geom.width_m, geom.height_m, 0), 0)
        with self.assertRaisesRegex(ValueError, "cutoff"):
            te10_power_for_peak_electric_field_w(1e9, geom.width_m, geom.height_m, e)

    def test_tm01_impedance_is_passive_and_requires_propagation(self):
        radius = self.case.circular.diameter_m / 2
        impedances = tm01_impedance_ohm(np.array([2.7e9, 3.0e9]), radius)
        self.assertGreater(impedances[0], 0)
        self.assertGreater(impedances[1], impedances[0])
        self.assertLess(impedances[1], math.sqrt(mu_0 / epsilon_0))
        with self.assertRaisesRegex(ValueError, "cutoff"):
            tm01_impedance_ohm(np.array([2.4e9]), radius)

    def test_taper_cross_sections_match_circle_and_rectangle_endpoints(self):
        radius = self.case.circular.diameter_m / 2
        a, b = self.case.rectangular.width_m, self.case.rectangular.height_m
        circle = transition_polygon(0, radius, a, b)
        rectangle = transition_polygon(1, radius, a, b)
        np.testing.assert_allclose(np.hypot(*circle), radius, rtol=1e-12)
        self.assertAlmostEqual(float(np.max(abs(rectangle[0]))), a / 2)
        self.assertAlmostEqual(float(np.max(abs(rectangle[1]))), b / 2)
        self.assertTrue(np.all((abs(rectangle[0]) <= a / 2 + 1e-14)))
        self.assertTrue(np.all((abs(rectangle[1]) <= b / 2 + 1e-14)))

    def test_drude_passivity_zero_density_and_phasor_heating(self):
        frequency = 3e9
        density = 1e18
        collision_rate = 2e11
        sigma = conductivity_s_m(density, collision_rate, frequency)
        eps = relative_permittivity(density, collision_rate, frequency)
        self.assertGreater(sigma.real, 0)
        self.assertGreater(eps.imag, 0)
        self.assertAlmostEqual(eps.imag, sigma.real / (epsilon_0 * 2 * math.pi * frequency))
        self.assertEqual(relative_permittivity(0, collision_rate, frequency), 1)
        self.assertAlmostEqual(
            peak_phasor_heating_w_m3(density, collision_rate, frequency, 1e5),
            0.5 * sigma.real * 1e10,
        )
        with self.assertRaises(ValueError):
            conductivity_s_m(-1, collision_rate, frequency)

    def test_te10_drude_slab_zero_density_loss_and_passivity(self):
        a = self.case.rectangular.width_m
        r0, t0 = te10_slab_scattering(3e9, a, 0.012, 0, 2e10)
        self.assertAlmostEqual(abs(r0), 0)
        self.assertAlmostEqual(abs(t0), 1)
        r, t = te10_slab_scattering(3e9, a, 0.012, 1e17, 2e10)
        self.assertGreater(abs(r), 0)
        self.assertGreater(1 - abs(r)**2 - abs(t)**2, 0)
        self.assertLessEqual(abs(r)**2 + abs(t)**2, 1)
        rc, tc = te10_slab_scattering(3e9, a, 0.012, 1e17, 0)
        self.assertAlmostEqual(abs(rc)**2 + abs(tc)**2, 1, delta=1e-12)


class ResourceBudgetTests(unittest.TestCase):
    def test_budget_is_quarter_of_available_capped_at_eight_gib(self):
        self.assertEqual(memory_budget_bytes(12 * GIB), 3 * GIB)
        self.assertEqual(memory_budget_bytes(40 * GIB), 8 * GIB)
        self.assertGreater(estimated_memory_bytes(100_000), 100_000 * 512)
        with self.assertRaises(ValueError):
            memory_budget_bytes(-1)

    def test_resonance_peak_is_reported_without_claiming_unbracketed_q(self):
        result = _resonance(np.arange(5.0), np.array([0.9, 0.95, 1.0, 0.9, 0.8]))
        self.assertEqual(result["frequency_hz"], 2.0)
        self.assertIsNone(result["loaded_q_half_power_approximate"])


if __name__ == "__main__":
    unittest.main()
