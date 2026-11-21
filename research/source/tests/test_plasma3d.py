"""Conservation and limit tests for the 3D kinetics/filament kernels and the PPT rate."""

import math
import unittest

import numpy as np

from rf_plasma import plasma3d as k3
from rf_plasma.ppt import adk_rate_s, keldysh_gamma, photons_required, ppt_rate_s

N_GAS = 1.9e26
LN_E0, NPTS = math.log(0.01), 400
H = (math.log(60.0) - LN_E0) / (NPTS - 1)
EPS = np.exp(LN_E0 + H * np.arange(NPTS))


def _params(**kw):
    prm = np.zeros(k3.N_PARAMS)
    prm[k3.P_MASS] = 39.948 * 1.66053906660e-27
    prm[k3.P_IP] = 15.76
    prm[k3.P_DOF] = 3.0
    prm[k3.P_PHOTON_EV] = 2.33
    prm[k3.P_OMEGA_L] = 2 * math.pi * 299792458.0 / 532e-9
    prm[k3.P_NPHOT] = 7
    prm[k3.P_ETA] = 0.3
    for name, val in kw.items():
        prm[getattr(k3, name)] = val
    return prm


def _tables(km=1e-13, iz=True, exc=True):
    te = EPS / 1.5
    tabs = np.zeros((4, NPTS))
    tabs[0] = km
    tabs[3] = 2 * 9.1093837015e-31 / (39.948 * 1.66053906660e-27) * km
    if iz:
        tabs[1] = 2.3e-14 * te**0.59 * np.exp(-17.44 / te)
    if exc:
        tabs[2] = 11.5 * 5e-15 * np.exp(-12.64 / te)
    return tabs


PPT_LN_I = np.linspace(math.log(1e12), math.log(1e19), 50)
PPT_TAB = np.full(50, -700.0)


def _ledger_balance(acc, w0, w1):
    lhs = acc[0] + acc[1]
    rhs = (w1 - w0) * k3.QE + acc[6] + acc[7] + acc[8]
    return lhs, rhs


class KineticCellTests(unittest.TestCase):
    def test_energy_ledger_closes_with_rf_heating_and_ionization(self):
        prm = _params(P_KCONV0=2.25e-43, P_B2_0=2.08e-12, P_B2_EXP=-0.67, P_K3B0=8.75e-39, P_K3B_EXP=-4.5)
        tabs = _tables()
        st = np.array([1e18, 0.0, 1e18 * 1.0])
        w0 = st[2]
        acc = np.zeros(9)
        k3.kinetic_cell(st, 2e-9, N_GAS, 300.0, 1e9 / k3.QE * 1e3, 0.0, prm, tabs, LN_E0, H,
                        PPT_LN_I[0], PPT_LN_I[1] - PPT_LN_I[0], PPT_TAB, acc)
        lhs, rhs = _ledger_balance(acc, w0, st[2])
        self.assertGreater(acc[0], 0)
        self.assertAlmostEqual(rhs / lhs, 1.0, delta=1e-9)
        self.assertTrue(np.all(st >= 0))

    def test_laser_ppt_and_ib_ledger(self):
        prm = _params()
        tabs = _tables()
        tab = np.log(np.full(50, 1e6))  # constant 1e6 s^-1 per neutral
        st = np.array([1e10, 0.0, 1e10 * 0.05])
        w0 = st[2]
        acc = np.zeros(9)
        k3.kinetic_cell(st, 1e-10, N_GAS, 300.0, 0.0, 1e16, prm, tabs, LN_E0, H,
                        PPT_LN_I[0], PPT_LN_I[1] - PPT_LN_I[0], tab, acc)
        lhs, rhs = _ledger_balance(acc, w0, st[2])
        self.assertAlmostEqual(rhs / lhs, 1.0, delta=1e-9)
        self.assertGreater(st[0], 1e10 + 0.5 * 1e6 * N_GAS * 1e-10)

    def test_recombination_matches_analytic_decay(self):
        beta = 1e-13
        prm = _params(P_B1_0=beta)
        tabs = _tables(iz=False, exc=False)
        n0 = 1e21
        st = np.array([n0, 0.0, n0 * 0.04])
        acc = np.zeros(9)
        t = 1e-7
        k3.kinetic_cell(st, t, N_GAS, 300.0, 0.0, 0.0, prm, tabs, LN_E0, H,
                        PPT_LN_I[0], PPT_LN_I[1] - PPT_LN_I[0], PPT_TAB, acc)
        self.assertAlmostEqual(st[0] / (n0 / (1 + beta * n0 * t)), 1.0, delta=0.02)

    def test_stiff_ion_conversion_conserves_charge_and_equilibrates(self):
        # N2-like: conversion ~2.5e12 /s, back-reaction 1e6 /s, one 15 ps step, no losses
        prm = _params(P_KCONV0=6.8e-41, P_KCONV_EXP=-1.64, P_KDISS0=1e6 / N_GAS, P_DISS_A=0.0)
        tabs = _tables(iz=False, exc=False)
        n0 = 1e15
        st = np.array([n0, 0.0, n0 * 0.04])
        acc = np.zeros(9)
        kconv = 6.8e-41 * N_GAS**2
        for _ in range(10):
            iters = k3.kinetic_cell(st, 1.5e-11, N_GAS, 300.0, 0.0, 0.0, prm, tabs, LN_E0, H,
                                    PPT_LN_I[0], PPT_LN_I[1] - PPT_LN_I[0], PPT_TAB, acc)
            self.assertLess(iters, 5)
        self.assertAlmostEqual((st[0] + st[1]) / n0, 1.0, delta=1e-12)
        self.assertAlmostEqual((st[0] / st[1]) / (1e6 / kconv), 1.0, delta=1e-3)

    def test_eedf_closure_switches_from_boltzmann_to_maxwell_with_ionization_degree(self):
        # Boltzmann rows: no ionization; Maxwell rows: ionization. Fixed eps (no heating, no losses).
        tb = _tables(iz=False, exc=False)
        tm = _tables(exc=False)
        tabs = np.vstack((tb, tm))
        prm = _params(P_EE_SCALE=1.0)
        growth = {}
        for x in (1e-12, 1e-1):
            n0 = x * N_GAS
            st = np.array([n0, 0.0, n0 * 4.0])
            k3.kinetic_cell(st, 1e-13, N_GAS, 300.0, 0.0, 0.0, prm, tabs, LN_E0, H,
                            PPT_LN_I[0], PPT_LN_I[1] - PPT_LN_I[0], PPT_TAB, np.zeros(9))
            growth[x] = (st[0] + st[1]) / n0 - 1.0
        kiz_m = 2.3e-14 * (4.0 / 1.5) ** 0.59 * math.exp(-17.44 / (4.0 / 1.5))
        full = kiz_m * N_GAS * 1e-13
        self.assertLess(growth[1e-12], 1e-3 * full)          # nu_ee << nu_eps: Boltzmann (no tail)
        self.assertAlmostEqual(growth[1e-1] / full, 1.0, delta=0.05)  # nu_ee >> nu_eps: Maxwellian

    def test_eedf_closure_weight_follows_energy_relaxation_not_a_fixed_ionization_degree(self):
        # same ne/N, eps: stronger inelastic loss (N2-like vibrational sink) keeps the EEDF Boltzmann
        elastic = np.vstack((_tables(iz=False, exc=False), _tables(exc=False)))
        lossy = elastic.copy()
        lossy[2] = 1e-13  # eV m^3/s energy-loss coefficient, ~4e12 /s relaxation at 4 eV
        ne = 1e-4 * N_GAS
        w_el = k3._wmax(ne, N_GAS, 4.0, elastic, LN_E0, H, 15.76, 1.0)
        w_ls = k3._wmax(ne, N_GAS, 4.0, lossy, LN_E0, H, 15.76, 1.0)
        nee = k3.nu_ei(ne, 4.0 / 1.5)
        nue = N_GAS * (1e-13 / 4.0 + elastic[3, 0])
        self.assertGreater(w_el, 0.99)
        self.assertAlmostEqual(w_ls, nee / (nee + nue), delta=1e-3)
        self.assertLess(w_ls, 0.5)
        self.assertEqual(k3._wmax(ne, N_GAS, 4.0, elastic, LN_E0, H, 15.76, 0.0), 0.0)

    def test_ionization_depletes_neutrals_and_never_overdraws_electron_energy(self):
        # strong RF heating with a huge ionization rate: ions cannot exceed the nuclei available
        prm = _params(P_KCONV0=2.25e-43)
        tabs = _tables()
        tabs[1] *= 1e3
        n_g = 1e20
        st = np.array([1e18, 0.0, 1e18 * 10.0])
        w0 = st[2]
        acc = np.zeros(9)
        k3.kinetic_cell(st, 1e-9, n_g, 300.0, 1e10 / k3.QE * 1e6, 0.0, prm, tabs, LN_E0, H,
                        PPT_LN_I[0], PPT_LN_I[1] - PPT_LN_I[0], PPT_TAB, acc)
        self.assertLessEqual(st[0] + 2 * st[1], (n_g + 1e18) * (1 + 1e-9))
        self.assertGreater(st[0] + st[1], 0.9 * n_g)        # nearly fully ionized
        self.assertGreater(st[2], 0.0)
        lhs, rhs = _ledger_balance(acc, w0, st[2])
        self.assertAlmostEqual(rhs / lhs, 1.0, delta=1e-9)


class FilamentTests(unittest.TestCase):
    def _filament(self, nr=40, te_core=2.0, ne_core=1e24, tg_core=300.0):
        rb = np.zeros((1, nr + 1))
        rb[0, 1:] = np.geomspace(1e-6, 1e-3 / math.sqrt(math.pi), nr)
        area = math.pi * np.diff(rb[0] ** 2)
        rc = 0.5 * (rb[0, 1:] + rb[0, :-1])
        prof = np.exp(-(rc / 20e-6) ** 2)
        heavy = (N_GAS * area)[None, :]
        tg = 300.0 + (tg_core - 300.0) * prof
        eg = (1.5 * k3.KB * tg * N_GAS * area)[None, :]
        cxp = (ne_core * prof * area)[None, :]
        cx2p = np.zeros_like(cxp)
        ue = cxp * 1.5 * te_core
        return dict(rb=rb, ub=np.zeros_like(rb), heavy=heavy, eg=eg, cxp=cxp, cx2p=cx2p, ue=ue,
                    ring_work=np.zeros_like(cxp))

    def _advance(self, f, prm, tabs, dt, steps, hydro):
        ledger = np.zeros((1, 12))
        p0 = N_GAS * k3.KB * 300.0
        for _ in range(steps):
            k3.filament_update(f["rb"], f["ub"], f["heavy"], f["eg"], f["cxp"], f["cx2p"], f["ue"],
                               f["ring_work"], np.ones(1), np.zeros(1, dtype=np.int64), np.zeros(1), dt,
                               N_GAS, 300.0, p0, prm, tabs, LN_E0, H, PPT_LN_I[0], PPT_LN_I[1] - PPT_LN_I[0],
                               PPT_TAB, np.zeros(1), np.full(1, 20e-6), ledger, hydro, 1e-3)
        return ledger

    def test_laser_march_never_absorbs_more_than_incident(self):
        # three optically thick edges in series: same-block depletion must bound total absorption
        n_e, nr, path = 3, 40, 1e-3
        one = self._filament(te_core=3.0, ne_core=2e26)
        f = {k: np.repeat(v, n_e, axis=0) for k, v in one.items()}
        prm = _params()
        tabs = _tables(iz=False, exc=False)
        trans = np.ones(1)
        laser_in = np.zeros(n_e)
        ledger = np.zeros((n_e, 12))
        w = 20e-6
        i_axis = np.full(n_e, 1e17)
        dt = 1e-11
        k3.filament_update(f["rb"], f["ub"], f["heavy"], f["eg"], f["cxp"], f["cx2p"], f["ue"], f["ring_work"],
                           trans, np.arange(n_e, dtype=np.int64), laser_in, dt, N_GAS, 300.0,
                           N_GAS * k3.KB * 300.0, prm, tabs, LN_E0, H, PPT_LN_I[0], PPT_LN_I[1] - PPT_LN_I[0],
                           PPT_TAB, i_axis, np.full(n_e, w), ledger, False, path)
        e_beam = 1e17 * 0.5 * math.pi * w * w * dt
        absorbed = laser_in.sum() * path
        self.assertGreater(absorbed, 0.5 * e_beam)                   # thick: most of the beam absorbed
        self.assertLessEqual(absorbed, e_beam * (1 + 1e-9))
        self.assertAlmostEqual(absorbed + trans[0] * e_beam, e_beam, delta=1e-9 * e_beam)
        self.assertGreater(laser_in[0], laser_in[2])                 # upstream edge shields downstream

    def test_diffusion_conserves_charge(self):
        f = self._filament()
        prm = _params(P_MU_XP=1.5e-4, P_MU_X2P=1.8e-4)
        tabs = _tables(iz=False, exc=False)
        q0 = f["cxp"].sum()
        self._advance(f, prm, tabs, 1e-9, 20, False)
        self.assertAlmostEqual(f["cxp"].sum() / q0, 1.0, delta=1e-9)
        prof = f["cxp"][0] / (math.pi * np.diff(f["rb"][0] ** 2))
        self.assertLess(prof[0], 1e24)  # core density spread outward

    def test_hydro_conserves_total_energy(self):
        f = self._filament(te_core=0.03, ne_core=1e10, tg_core=3000.0)
        prm = _params()
        tabs = _tables(km=0.0, iz=False, exc=False)
        e0 = f["eg"].sum() + f["ue"].sum() * k3.QE
        ledger = self._advance(f, prm, tabs, 2e-10, 50, True)
        e1 = f["eg"].sum() + f["ue"].sum() * k3.QE + ledger[0, 9] + ledger[0, 10]
        excess = e0 - (N_GAS * k3.KB * 300.0 * 1.5) * math.pi * f["rb"][0, -1] ** 2
        self.assertGreater(f["rb"][0, 1], 1e-6)  # hot core expanded
        self.assertLess(abs(e1 - e0) / excess, 1e-9)  # compatible (time-centred) work: exact


class PPTTests(unittest.TestCase):
    def test_tunnel_limit_matches_adk(self):
        for lam, inten in [(10e-6, 1e18), (10e-6, 3e18)]:
            self.assertLess(keldysh_gamma(inten, lam, 15.76), 0.1)
            self.assertAlmostEqual(ppt_rate_s(inten, lam, 15.76) / adk_rate_s(inten, 15.76), 1.0, delta=0.06)

    def test_multiphoton_slope_equals_photon_number(self):
        k = photons_required(1e15, 532e-9, 15.76)
        r1, r2 = ppt_rate_s(1e15, 532e-9, 15.76), ppt_rate_s(1.01e15, 532e-9, 15.76)
        self.assertEqual(k, 7)
        self.assertAlmostEqual(math.log(r2 / r1) / math.log(1.01), k, delta=0.05)


if __name__ == "__main__":
    unittest.main()
