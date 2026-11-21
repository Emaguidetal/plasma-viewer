"""Validation of the custom 3D FDTD core: TE10 port bookkeeping and a Drude slab."""

import math
import unittest

import numba
import numpy as np

from rf_plasma.drude import te10_slab_scattering
from rf_plasma.fdtd3d import FDTD, _plane_projection, _poynting_z
from rf_plasma.geometry3d import TE10Port, raised_cosine_ramp, straight_guide

numba.set_num_threads(2)
F0 = 3e9
A_WR284, B_WR284 = 72.136e-3, 34.036e-3


def _run_guide(slab=None, steps=3200, length=0.15):
    lay = straight_guide(1e-3, A_WR284, B_WR284, length)
    g = lay.grid
    fd = FDTD(g, lay.node_inside, lay.pml_lo, lay.pml_hi, F0)
    if slab is not None:
        density, nu, k0, nk = slab
        wp2 = np.zeros(lay.node_inside.shape)
        nun = np.zeros(lay.node_inside.shape)
        wp2[:, :, k0:k0 + nk] = density * 1.602176634e-19**2 / (8.8541878128e-12 * 9.1093837015e-31)
        nun[:, :, k0:k0 + nk] = nu
        fd.set_bulk_plasma(wp2, nun)
    port = TE10Port(lay, F0)
    port.amplitude = 1.0
    port.envelope = raised_cosine_ramp(3e-9)
    k_out = lay.trans_k
    rec = {n: np.zeros(steps) for n in ("v_in", "i_in", "v_out", "i_out", "s_out")}
    for n in range(steps):
        fd.step(port.inject)
        rec["v_in"][n], rec["i_in"][n] = _plane_projection(fd.Ey, fd.Hx, lay.meas_k, port.profile, port.profile)
        rec["v_out"][n], rec["i_out"][n] = _plane_projection(fd.Ey, fd.Hx, k_out, port.profile, port.profile)
        rec["s_out"][n] = _poynting_z(fd.Ex, fd.Ey, fd.Hx, fd.Hy, k_out) * g.dx**2
    fin, bin_ = port.modal_amplitudes(rec["v_in"], rec["i_in"], g.dt)
    fout, bout = port.modal_amplitudes(rec["v_out"], rec["i_out"], g.dt)
    period = int(round(1 / (F0 * g.dt)))
    flux = np.convolve(rec["s_out"], np.ones(period) / period, mode="valid")
    tail = slice(-period, None)
    return dict(port=port, fin=fin[tail].mean(), bin=bin_[tail].mean(), fout=fout[tail].mean(),
                bout=bout[tail].mean(), flux=flux[tail].mean())


class FDTDPortTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.empty = _run_guide()

    def test_matched_guide_reflection_and_flux(self):
        r = self.empty
        self.assertLess(abs(r["bin"]) / abs(r["fin"]), 2e-3)
        modal = float(r["port"].power(np.array([r["fout"]]))[0])
        self.assertAlmostEqual(r["flux"] / modal, 1.0, delta=0.01)

    def test_drude_slab_matches_analytic(self):
        density, nu, thick = 1e17, 2e10, 12
        lay = straight_guide(1e-3, A_WR284, B_WR284, 0.15)
        k0 = (lay.meas_k + lay.trans_k) // 2 - thick // 2
        s = _run_guide((density, nu, k0, thick))
        s11 = abs(s["bin"]) / abs(self.empty["fin"])
        s21 = abs(s["fout"]) / abs(self.empty["fout"])
        r, t = te10_slab_scattering(F0, 72e-3, thick * 1e-3, density, nu)
        # openEMS at the same 1 mm mesh gave |S11| 0.2666, |S21| 0.7695 (analytic 0.2504 / 0.7862)
        self.assertAlmostEqual(s11, abs(r), delta=0.02)
        self.assertAlmostEqual(s21, abs(t), delta=0.02)
        self.assertLessEqual(s11**2 + s21**2, 1.0)
        print(f"\nslab |S11| {s11:.5f} (analytic {abs(r):.5f}), |S21| {s21:.5f} (analytic {abs(t):.5f})")


if __name__ == "__main__":
    unittest.main()
