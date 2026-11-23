import math
import unittest

import numpy as np

from rf_plasma import plasma_story3d, rod3d, rod_drive, seed_sources


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


class StoryViewerTests(unittest.TestCase):
    def test_unbridged_case_reports_evolution_without_claiming_a_bridge(self):
        import json, tempfile
        from pathlib import Path
        from unittest.mock import patch

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = root / "rod3d" / "example"
            case.mkdir(parents=True)
            summary = {
                "config": {"e0_v_m": 1e7, "dx_c": 1e-3, "r_s": 1e-5, "x_focus_m": 0.0},
                "t_first_birth_ns": 1.5, "t_bridge_ns": None, "n_rods": 1,
            }
            (case / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
            np.savez(
                case / "rod3d.npz",
                pos=np.array([[0, 0.017, 0], [0, 0.018, 0]]),
                rad=np.array([1e-5, 1e-5]), kind=np.array([0, 1]), rod=np.array([-1, 0]),
                lam=np.full((3, 2), 1e8), t_s=np.array([0.0, 1e-9, 2e-9]),
                seed_electrons=np.array([100.0]), m_frame=np.array([2, 2, 2]),
                hist_e0_v_m=np.full(3, 1e7),
            )
            with patch.object(plasma_story3d, "ROOT", root):
                result = plasma_story3d.case_data("example")
            phase = next(frame["phase"] for frame in result["frames"] if frame["t"] == 2.0)
        self.assertIn("no wall-to-wall bridge", phase)
        self.assertIn("display-only", result["frames"][0]["phase"])
        self.assertIn("t=0", result["frames"][0]["phase"])
        self.assertIn("Pre-zero seed buildup is display-only", plasma_story3d.HTML)


class RodDriveTests(unittest.TestCase):
    def test_line_capacitance_equals_cavity_field_capacitance(self):
        cc = rod_drive.cavity_constants(rod_drive.DriveConfig())
        self.assertAlmostEqual(cc["cap"] / cc["c_eps"], 1.0, delta=1e-3)

    def test_short_drive_conserves_energy(self):
        c = rod_drive.DriveConfig(name="unittest", t_end_s=3e-9, nr=60, samples=31, e_start_v_m=1.5e7)
        s = rod_drive.run(c, save=False, verbose=False)
        self.assertLess(abs(s["circuit_energy_residual_rel"]), 1e-9)
        self.assertLess(abs(s["plasma_energy_residual_rel"]), 1e-9)

    def test_from_rod3d_hands_off_bridged_rod_at_bridge_time(self):
        import json, tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "case"
            d.mkdir()
            (d / "summary.json").write_text(json.dumps({"t_bridge_ns": 4.0, "config": {"e0_v_m": 2e7}}))
            # node 0: laser line; rod 0 spans 1 mm, rod 1 (the bridge) spans 3 mm along y
            pos = np.array([[0, 17e-3, 0], [0, 16.5e-3, 0], [0, 17.5e-3, 0], [1e-3, 15.5e-3, 0], [1e-3, 18.5e-3, 0]])
            rad = np.array([20e-6, 1e-4, 1e-4, 2e-4, 4e-4])
            ne = np.array([1e20, 1e21, 1e21, 1e22, 1e24])
            lam = ne * math.pi * rad ** 2
            t = np.array([0.0, 2e-9, 4e-9, 6e-9])
            np.savez(d / "rod3d.npz", t_s=t, m_frame=np.full(4, 5), kind=np.array([0, 1, 1, 1, 1]),
                     rod=np.array([-1, 0, 0, 1, 1]), lam=np.tile(lam, (4, 1)), rad=rad, pos=pos,
                     hist_t_s=t, hist_e0_v_m=np.array([2e7, 1.9e7, 1.5e7, 1.0e7]))
            h = rod_drive.from_rod3d("case", t_off_after_laser_s=100e-9, root=Path(tmp))
        self.assertAlmostEqual(h["e_start_v_m"], 1.5e7)
        lam_h = 2 / (1 / lam[3] + 1 / lam[4])
        self.assertAlmostEqual(h["n0_m3"] * math.pi * h["r0_m"] ** 2 / lam_h, 1.0, places=9)
        self.assertAlmostEqual(h["drive_off_s"], 96e-9)
        self.assertGreater(h["r0_m"], 2e-4)
        self.assertLess(h["r0_m"], 4e-4)


if __name__ == "__main__":
    unittest.main()
