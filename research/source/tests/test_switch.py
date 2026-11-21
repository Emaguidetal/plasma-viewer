import math
import unittest

import numpy as np

from rf_plasma import rod3d, rod_drive, seed_sources


class Rod3DTests(unittest.TestCase):
    def test_static_dipole_matches_slender_body(self):
        p, analytic = rod3d.static_dipole(L=5e-3, r=50e-6, h=0.25e-3)
        self.assertAlmostEqual(p / analytic, 1.0, delta=0.05)

    def test_tip_law_bounded_beyond_calibration_range(self):
        tp = rod3d.TipPhysics(3.574e26, 5.0, 0.85, 10e-6, 0.3e-3)
        r, et, v = tp.tip(np.array([1e3, 3e5]), 2e7, 50e-6)
        self.assertTrue(np.all(r <= 0.3e-3))
        self.assertTrue(np.all(et <= 2e7 + tp.e_m * (1 + 1e-12)))
        self.assertLess(v.max(), 1e8)

    def test_overlapping_sites_merge_and_p_stays_positive_definite(self):
        cfg = rod3d.RodConfig(e0_v_m=2e7, n_inc=1e20, c_v=0.85, e_m_td=5.0, x_half=1.5e-3, stochastic=False)
        m = rod3d.Rod3D(cfg)
        while m.t < 1.6e-9:
            m.step(5e-12)
        xs = sorted(r["x"] for r in m.rods)
        self.assertGreater(len(xs), 0)
        self.assertTrue(all(b - a >= 2 * cfg.r_max - 1e-12 for a, b in zip(xs, xs[1:])))
        self.assertGreater(np.linalg.eigvalsh(m.P[: m.M, : m.M]).min(), 0.0)


class SeedSourceTests(unittest.TestCase):
    def test_impurity_mpi_dominates_helium_and_scales_with_photon_order(self):
        b = seed_sources.Budget()
        fns = seed_sources.per_source(b)
        self.assertGreater(fns["MPI H2O"](b), 1e6 * fns["MPI He"](b))
        k = seed_sources.energy_slope(b, fns["MPI H2O"])
        self.assertGreater(k, 4.0)
        self.assertLessEqual(k, 6.5)


class RodDriveTests(unittest.TestCase):
    def test_line_capacitance_equals_cavity_field_capacitance(self):
        cc = rod_drive.cavity_constants(rod_drive.DriveConfig())
        self.assertAlmostEqual(cc["cap"] / cc["c_eps"], 1.0, delta=1e-3)

    def test_short_drive_conserves_energy(self):
        c = rod_drive.DriveConfig(name="unittest", t_end_s=3e-9, nr=60, samples=31)
        s = rod_drive.run(c, save=False, verbose=False)
        self.assertLess(abs(s["circuit_energy_residual_rel"]), 0.02)
        self.assertLess(abs(s["plasma_energy_residual_rel"]), 1e-9)


if __name__ == "__main__":
    unittest.main()
