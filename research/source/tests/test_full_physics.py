"""Tests for the extended kernels (plasma3d_cr): CR chemistry, bulk transport/heating, propagation."""

import math
import unittest

import numpy as np

from rf_plasma import plasma3d as k3
from rf_plasma import plasma3d_cr as kc
from rf_plasma.gas_model import load_gas_model

LAM = 532e-9
N0 = 1e6 / (k3.KB * 300.0)


def _gm(cr=True, transport=True):
    return load_gas_model("He", LAM, allow_unverified=True, extended=True, cr=cr, transport=transport)


def _closure(out, w0, w1):
    lhs = out[0] + out[1]
    rhs = (w1 - w0) * k3.QE + out[6] + out[7] + out[8]
    return lhs, rhs


class CRKineticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.on = _gm(cr=True)
        cls.off = _gm(cr=False)

    def _run(self, gm, st, p_ext, dt=2e-9, inten=0.0):
        out = np.zeros(9)
        kc.kinetic_cell_cr(st, dt, N0, 300.0, p_ext, inten, gm.prm, gm.tabs, gm.ln_e0, gm.h,
                           gm.ln_i0, gm.hi, gm.ppt_tab, out, 0.0)
        return out

    def test_cr_off_reproduces_base_kernel(self):
        gm = self.off
        p_ext = 5e7 / k3.QE * 1e3  # 5e10 W/m^3
        s3 = np.array([1e18, 0.0, 1e18 * 1.0])
        o3 = np.zeros(9)
        k3.kinetic_cell(s3, 2e-9, N0, 300.0, p_ext, 0.0, gm.prm[:k3.N_PARAMS], gm.tabs, gm.ln_e0, gm.h,
                        gm.ln_i0, gm.hi, gm.ppt_tab, o3, 0.0)
        s6 = np.array([1e18, 0.0, 1e18, 0.0, 0.0, 0.0])
        o6 = self._run(gm, s6, p_ext)
        np.testing.assert_allclose(s6[:3], s3, rtol=1e-9)
        np.testing.assert_allclose(o6, o3, rtol=1e-9, atol=1e-12 * abs(o3).max())
        self.assertTrue(np.all(s6[3:] == 0.0))

    def test_energy_closure_and_positivity_with_cr(self):
        gm = self.on
        for p_ext in (0.0, 5e10 / k3.QE, 5e12 / k3.QE):
            st = np.array([1e21, 1e20, 1e21 * 2.0, 1e21, 1e20, 0.0])
            w0 = st[2]
            out = self._run(gm, st, p_ext)
            lhs, rhs = _closure(out, w0, st[2])
            scale = max(abs(lhs), abs(out[6]), abs(out[7]), abs(out[8]), 1e-6)
            self.assertLess(abs(lhs - rhs) / scale, 1e-9, (p_ext, lhs, rhs))
            self.assertTrue(np.all(st >= 0.0) and np.all(np.isfinite(st)))
        self.assertGreater(out[8], 0.0)  # strong heating stores ionization/excitation energy

    def test_rf_excitation_builds_metastables(self):
        st = np.array([1e20, 0.0, 1e20 * 1.0, 0.0, 0.0, 0.0])
        self._run(self.on, st, 1e12 / k3.QE, dt=1e-9)
        self.assertGreater(st[3] + st[4], 1e17)

    def test_empty_cell_creates_no_electrons(self):
        """The numerical floor must not inject a seed (it would avalanche into a spurious volume breakdown)."""
        st = np.zeros(6)
        self._run(self.on, st, 0.0, dt=1e-9)
        self.assertLessEqual(st[0] + st[1], 10 * k3.NE_FLOOR)

    def test_laser_seed_line_density(self):
        from rf_plasma.coupled3d import Case, _impurity_seed
        w = np.array([20e-6, 40e-6])
        rb = np.tile(np.concatenate([[0.0], np.geomspace(1e-6, 1e-3, 400)]), (2, 1))
        c = Case(gas="He", wavelength_m=LAM, laser_seed_m3=1e18, impurity_frac=0.0)
        cnt, exc = _impurity_seed(c, rb, w, 1e3)
        np.testing.assert_allclose(cnt.sum(axis=1), 1e18 * math.pi * w ** 2 / 2, rtol=2e-3)
        self.assertTrue(np.all(exc > 0))
        # effective order k: reference pulse at w_ref gives n_ref*pi w^2/(2k); doubling energy scales by 2^k
        ck = Case(gas="He", wavelength_m=LAM, laser_seed_m3=1e18, laser_seed_k=6.0, laser_seed_ref=(1e-3, 1e-9, 20e-6),
                  fwhm_s=1e-9, energy_j=1e-3)
        p_pk = 1e-3 * math.sqrt(4 * math.log(2) / math.pi) / 1e-9
        c1 = _impurity_seed(ck, rb[:1], w[:1], p_pk)[0].sum()
        c2 = _impurity_seed(ck, rb[:1], w[:1], 2 * p_pk)[0].sum()
        self.assertAlmostEqual(c1 / (1e18 * math.pi * w[0] ** 2 / 12), 1.0, delta=5e-3)
        self.assertAlmostEqual(c2 / c1, 64.0, delta=1e-6 * 64)

    def test_superelastic_heating_of_afterglow(self):
        """Without any input, a He* population heats cold electrons above the CR-off afterglow."""
        cold = 0.05
        a = np.array([1e20, 0.0, 1e20 * cold, 1e22, 0.0, 0.0])
        b = np.array([1e20, 0.0, 1e20 * cold, 0.0, 0.0, 0.0])
        self._run(self.on, a, 0.0, dt=0.5e-9)
        self._run(self.on, b, 0.0, dt=0.5e-9)
        ne_a = a[0] + a[1] + 2 * a[5]
        ne_b = b[0] + b[1]
        self.assertGreater(a[2] / ne_a, 3.0 * b[2] / ne_b)


class TransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gm = _gm()

    def _box(self, n=21, ne_peak=1e12, eps=1.0):
        shape = (n, n, n)
        x = (np.arange(n) - n // 2)
        r2 = x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2
        ne = 1e6 + ne_peak * np.exp(-r2 / 4.0)
        st = np.zeros((6,) + shape)
        st[0] = ne
        st[2] = ne * eps
        return st, np.full(shape, 300.0), np.ones(shape, np.bool_), r2

    def test_conservation_and_free_diffusion_rate(self):
        gm = self.gm
        st, tg, act, r2 = self._box(ne_peak=1e10)
        dx = 1e-3
        coef = np.zeros_like(st)
        delta = np.zeros((5,) + st.shape[1:])
        n_tot0, w_tot0 = st[0].sum(), st[2].sum()
        var0 = ((st[0] - 1e6) * r2).sum() / (st[0] - 1e6).sum()
        dt = 1e-6
        dmax = kc.bulk_transport(st, tg, act, dt, dx, N0, gm.prm, gm.tabs, gm.ln_e0, gm.h, coef, delta, True)
        d_e = coef[2, 10, 10, 10]
        self.assertGreater(d_e, 1e-3)
        var1 = ((st[0] - 1e6) * r2).sum() / (st[0] - 1e6).sum()
        # uniform eps and lambda_D >> dx: discrete Laplacian gives exactly d<r^2>/dt = 6 D_e
        np.testing.assert_allclose((var1 - var0) * dx * dx, 6 * d_e * dt, rtol=2e-3)
        self.assertGreaterEqual(dmax, d_e)
        for _ in range(20):
            kc.bulk_transport(st, tg, act, dt, dx, N0, gm.prm, gm.tabs, gm.ln_e0, gm.h, coef, delta, True)
        self.assertAlmostEqual(st[0].sum() / n_tot0, 1.0, places=12)
        self.assertAlmostEqual(st[2].sum() / w_tot0, 1.0, places=12)

    def test_ambipolar_limit_slows_dense_plasma(self):
        gm = self.gm
        st, tg, act, _ = self._box(ne_peak=1e20)
        coef = np.zeros_like(st)
        delta = np.zeros((5,) + st.shape[1:])
        kc.bulk_transport(st, tg, act, 1e-12, 1e-3, N0, gm.prm, gm.tabs, gm.ln_e0, gm.h, coef, delta, False)
        te, ti = 1.0 / 1.5, 300.0 / k3.EV_K
        d_a = gm.prm[k3.P_MU_XP] * k3.N_LOSCHMIDT / N0 * (te + ti)
        np.testing.assert_allclose(coef[0, 10, 10, 10], d_a, rtol=1e-3)
        self.assertLess(coef[0, 10, 10, 10], 0.02 * coef[2, 10, 10, 10])

    def test_gas_heating_matches_ledger(self):
        gm = self.gm
        shape = (3, 3, 3)
        st = np.zeros((6,) + shape)
        st[0] = 1e21
        st[2] = 1e21 * 1.0
        tg = np.full(shape, 300.0)
        act = np.ones(shape, np.bool_)
        work = np.full(shape, 1e3)  # J/m^3 over the step
        led = np.zeros((3, 9))
        wp2 = np.zeros(shape); nu = np.zeros(shape)
        kc.bulk_update_cr(st, tg, act, work, 1e-9, N0, gm.prm, gm.tabs, gm.ln_e0, gm.h, gm.ln_i0, gm.hi,
                          gm.ppt_tab, wp2, nu, led, True)
        heat = (0.5 * gm.prm[k3.P_DOF] * k3.KB * N0 * (tg - 300.0)).sum()
        np.testing.assert_allclose(heat, led[:, 6].sum(), rtol=1e-3)
        self.assertGreater(tg.max(), 300.0)
        np.testing.assert_allclose(led[:, 0].sum(), work.sum(), rtol=1e-9)


class PropagationTests(unittest.TestCase):
    def test_cn_vacuum_gaussian(self):
        lam, w0, dr = LAM, 20e-6, 1e-6
        k = 2 * math.pi / lam
        zr = math.pi * w0 ** 2 / lam
        npg = 800
        r = (np.arange(npg) + 0.5) * dr
        a = np.exp(-(r / w0) ** 2).astype(np.complex128)
        dA = math.pi * (((np.arange(npg) + 1) * dr) ** 2 - (np.arange(npg) * dr) ** 2)
        p0 = (abs(a) ** 2 * dA).sum()
        ws = [np.empty(npg, np.complex128) for _ in range(5)]
        nst = 80
        for _ in range(nst):
            kc._cn_radial(a, dr, zr / nst, k, *ws)
        p1 = (abs(a) ** 2 * dA).sum()
        self.assertAlmostEqual(p1 / p0, 1.0, places=10)
        w2 = 2 * (abs(a) ** 2 * r ** 2 * dA).sum() / p1  # second-moment width^2
        np.testing.assert_allclose(math.sqrt(w2), math.sqrt(2) * w0, rtol=0.01)

    def test_filament_diffusion_conserves_energy_with_te_gradient(self):
        gm = _gm()
        rings = 24
        edges = np.concatenate([np.linspace(0, 60e-6, 13), np.geomspace(60e-6, 1e-3 / math.sqrt(math.pi), 13)[1:]])
        rb = edges[None, :].copy()
        area = math.pi * np.diff(rb ** 2, axis=1)
        rc = 0.5 * (rb[:, 1:] + rb[:, :-1])
        cs6 = np.zeros((6, 1, rings))
        cs6[0] = (1e9 + 1e24 * np.exp(-(rc / 15e-6) ** 2)) * area
        cs6[5] = 0.1 * cs6[0]
        te = 0.5 + 3.0 * np.exp(-(rc / 30e-6) ** 2)
        cs6[2] = 1.5 * te * (cs6[0] + 2 * cs6[5])
        u0, n0 = cs6[2].sum(), (cs6[0] + 2 * cs6[5]).sum()
        kc._fil_diffuse_cr(0, rb, N0 * area, cs6, 5e-12, gm.prm, gm.tabs, gm.ln_e0, gm.h)
        self.assertAlmostEqual((cs6[0] + 2 * cs6[5]).sum() / n0, 1.0, delta=1e-13)
        self.assertAlmostEqual(cs6[2].sum() / u0, 1.0, delta=1e-13)
        self.assertLess(cs6[0, 0, 0], (1e9 + 1e24 * np.exp(-(rc[0, 0] / 15e-6) ** 2)) * area[0, 0])

    def test_filament_free_transport_front(self):
        """Allis-Rose filament transport: conservative, ambipolar in a dense core, free in a dilute one."""
        edges = np.concatenate([np.linspace(0, 60e-6, 13), np.geomspace(60e-6, 1e-3 / math.sqrt(math.pi), 13)[1:]])
        rb = edges[None, :].copy()
        area = math.pi * np.diff(rb ** 2, axis=1)
        rc = 0.5 * (rb[:, 1:] + rb[:, :-1])
        out = {}
        for peak in (1e16, 1e24):
            for free in (0.0, 1.0):
                gm = _gm()
                gm.prm[kc.P_FIL_FREE] = free
                cs6 = np.zeros((6, 1, 24))
                cs6[0] = (1e9 + peak * np.exp(-(rc / 15e-6) ** 2)) * area
                cs6[2] = 1.5 * 2.0 * cs6[0]
                n0, u0 = cs6[0].sum(), cs6[2].sum()
                for _ in range(20):
                    kc._fil_diffuse_cr(0, rb, N0 * area, cs6, 1e-10, gm.prm, gm.tabs, gm.ln_e0, gm.h)
                self.assertAlmostEqual(cs6[0].sum() / n0, 1.0, delta=1e-12)
                self.assertAlmostEqual(cs6[2].sum() / u0, 1.0, delta=1e-12)
                out[peak, free] = cs6[0, 0, rc[0] > 30e-6].sum() / n0  # fraction beyond 30 um
        self.assertGreater(out[1e16, 1.0], 3.0 * out[1e16, 0.0])        # dilute: free diffusion spreads
        self.assertAlmostEqual(out[1e24, 1.0] / out[1e24, 0.0], 1.0, delta=0.05)  # dense: stays ambipolar

    def test_filament_diffusion_near_empty_rings(self):
        """Round-off counts (~1e-40 /m) with inconsistent energy must not blow up the free-electron D."""
        edges = np.concatenate([np.linspace(0, 60e-6, 13), np.geomspace(60e-6, 1e-3 / math.sqrt(math.pi), 13)[1:]])
        rb = edges[None, :].copy()
        area = math.pi * np.diff(rb ** 2, axis=1)
        gm = _gm()
        gm.prm[kc.P_FIL_FREE] = 1.0
        cs6 = np.zeros((6, 1, 24))
        cs6[0, 0] = np.geomspace(1e-37, 1e-109, 24)
        cs6[2, 0] = 1e-42  # eps = ue/cnt up to 1e67 eV
        cs6[0, 0, :4] = 1e24 * area[0, :4]
        eps_top = math.exp(gm.ln_e0 + gm.h * (gm.tabs.shape[1] - 1))
        cs6[2, 0, :4] = 3.0 * eps_top * cs6[0, 0, :4]  # a hot core above the table ceiling
        u0 = cs6[2].sum()
        for _ in range(5):
            kc._fil_diffuse_cr(0, rb, N0 * area, cs6, 1e-11, gm.prm, gm.tabs, gm.ln_e0, gm.h)
        self.assertTrue(np.all(np.isfinite(cs6)))
        self.assertTrue(np.all(cs6[0] >= 0.0))
        self.assertAlmostEqual(cs6[2].sum() / u0, 1.0, delta=1e-12)

    def test_ee_table_level_interpolation(self):
        """e-e closure: ionization degree maps to table levels (linear below 1e-7, log10 above)."""
        prm = np.zeros(kc.N_PARAMS_CR)
        nlev, rows = 4, 4
        prm[kc.P_EE_NLEV] = nlev
        prm[kc.P_EE_ROW0] = 8
        ln_e0, h = 0.0, 1.0
        tabs = np.zeros((8 + rows * nlev, 4))
        for j in range(nlev):
            tabs[8 + rows * j] = 10.0 ** j  # row 0 of level j, flat in energy
        # values are interpolated linearly between adjacent levels
        for x, want in ((0.0, 1.0), (0.5e-7, 5.5), (1e-7, 10.0), (10 ** -6.5, 55.0), (1e-5, 1e3), (1.0, 1e3)):
            w = kc._wsel(x * 1e25, 1e25, 1.0, tabs, ln_e0, h, prm)
            self.assertLess(w, 0.0)
            v = kc._rx(tabs, 0, w, ln_e0, h, 1.0, prm)
            self.assertAlmostEqual(v / want, 1.0, places=9)

    def test_filament_energy_accounting_and_defocusing(self):
        gm = _gm()
        res = {}
        for ne0 in (1e9, 3e24):
            res[ne0] = self._fil_run(gm, ne0)
        for ne0, (inp, absd, trans, lat, axis) in res.items():
            np.testing.assert_allclose(absd + trans + lat, inp, rtol=1e-6)
        self.assertGreater(res[3e24][1], 1e3 * res[1e9][1])            # dense plasma absorbs
        self.assertLess(res[3e24][4], 0.8 * res[1e9][4])               # and defocuses (lower axis intensity)

    @staticmethod
    def _fil_run(gm, ne0):
        n_e, rings, dx = 9, 24, 1e-3
        edges = np.concatenate([np.linspace(0, 60e-6, rings // 2 + 1),
                                np.geomspace(60e-6, dx / math.sqrt(math.pi), rings // 2 + 1)[1:]])
        rb = np.tile(edges, (n_e, 1))
        ub = np.zeros_like(rb)
        area = math.pi * np.diff(rb ** 2, axis=1)
        heavy = N0 * area
        eg = 1.5 * k3.KB * 300.0 * heavy
        cs6 = np.zeros((6, n_e, rings))
        rc = 0.5 * (rb[:, 1:] + rb[:, :-1])
        cs6[0] = (1e9 + ne0 * np.exp(-(rc / 15e-6) ** 2)) * area
        cs6[2] = cs6[0] * 1.0
        ring_work = np.zeros((n_e, rings))
        order = np.arange(n_e, dtype=np.int64)
        laser_in = np.zeros(n_e)
        led = np.zeros((n_e, 12))
        dr, npg = 1e-6, 800
        r = (np.arange(npg) + 0.5) * dr
        k = 2 * math.pi / LAM
        w0 = 20e-6
        zr = math.pi * w0 ** 2 / LAM
        q = -4.5e-3 - 1j * zr
        field = np.exp(1j * k * r ** 2 / (2 * q)).astype(np.complex128)
        dA = math.pi * (((np.arange(npg) + 1) * dr) ** 2 - (np.arange(npg) * dr) ** 2)
        p_in = 1e8
        field *= math.sqrt(p_in / (abs(field) ** 2 * dA).sum())
        mask = np.ones(npg)
        n_abs = int(0.2 * npg)
        mask[-n_abs:] = np.cos(0.5 * math.pi * np.linspace(0, 1, n_abs)) ** 0.05
        beam = np.zeros(2)
        dt = 1e-12
        n_crit = k3.EPS0 * k3.ME * (2 * math.pi * k3.C0 / LAM) ** 2 / k3.QE ** 2
        kc.filament_update_cr(rb, ub, heavy, eg, cs6, ring_work, order, laser_in, dt,
                              N0, 300.0, 1e6, gm.prm, gm.tabs, gm.ln_e0, gm.h, gm.ln_i0, gm.hi, gm.ppt_tab,
                              field, dr, k, 0.0, 0.0, n_crit, led, False, dx, True, mask, beam)
        absorbed = laser_in.sum() * dx
        axis = abs(field[0]) ** 2
        return p_in * dt, absorbed, beam[1] * dt, beam[0], axis


if __name__ == "__main__":
    unittest.main()
