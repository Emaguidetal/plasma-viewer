import tempfile
import unittest
from pathlib import Path

import numpy as np

from rf_plasma.coupled3d import Case, run_case


class LaserTriggerFramesTest(unittest.TestCase):
    def test_guide_frames_match_timeseries_and_laser_seeds_focus(self):
        case = Case(geometry="guide", gas="He", laser=True, rf=False, fwhm_s=100e-12, energy_j=10e-3,
                    t_laser_s=60e-12, t_end_s=0.15e-9, pressure_pa=1e6, frame_every_blocks=2,
                    snapshot_s=[], tag="unittest")
        with tempfile.TemporaryDirectory() as tmp:
            s = run_case(case, threads=2, out_root=Path(tmp))
            out = Path(tmp) / s["name"]
            fr = dict(np.load(out / "frames.npz"))
            ts = dict(np.load(out / "timeseries.npz"))
            nb = len(ts["t"])
            blocks = np.arange(0, nb, 2)
            self.assertEqual(len(fr["t"]), len(blocks))
            np.testing.assert_allclose(fr["t"], ts["t"][blocks])
            np.testing.assert_allclose(fr["fil_ne_axis"], ts["fil_ne_axis"][blocks], rtol=1e-6)
            self.assertEqual(fr["fil_ne_r_focus"].shape[1] + 1, fr["fil_rb_focus"].shape[1])
            self.assertEqual(float(np.abs(fr["ey_xz"]).max()), 0.0)
        e = int(np.argmin(np.abs(fr["fil_y"] / float(fr["dx"]) - 0.5 - int(fr["focus"][1]))))
        self.assertGreater(fr["fil_ne_axis"][-1, e], 10 * case.seed_m3)
        self.assertGreater(fr["fil_ne_axis"][-1, e], 10 * fr["fil_ne_axis"][-1, 0])


if __name__ == "__main__":
    unittest.main()
