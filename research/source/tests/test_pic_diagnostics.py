import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pic"))
from pic_dynamics_report import nrmse, window_average  # noqa: E402


class PhaseWindowDriftTest(unittest.TestCase):
    def test_window_average_recovers_drift_buried_in_thermal_noise(self):
        # 2000-macro mean v_z: SE ~2e4 m/s against a ~5e4 m/s RF drift amplitude.
        freq, dt = 2.95e9, 4.e-14
        t = np.arange(0., 2/freq, dt)
        drift = 5.e4*np.sin(2*np.pi*freq*t)
        noisy = drift + np.random.default_rng(1).normal(0., 2.e4, t.size)
        tw, vw, se = window_average(t, noisy, freq, t[-1])
        _, exact, _ = window_average(t, drift, freq, t[-1])
        last = tw >= t[-1]-1/freq
        self.assertGreater(nrmse(noisy[t >= t[-1]-1/freq], drift[t >= t[-1]-1/freq]), .3)
        self.assertLess(nrmse(vw[last], exact[last]), .05)
        self.assertTrue(np.allclose(np.median(se), 2.e4/np.sqrt(265), rtol=.2))


if __name__ == "__main__":
    unittest.main()
