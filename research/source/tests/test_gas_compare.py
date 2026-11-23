"""Focused tests for the He/Ar/N2 swarm tables and gas comparison (rf_plasma.swarm_gases, gas_compare)."""
import math
import unittest

import numpy as np

from rf_plasma import gas_compare as gc
from rf_plasma import swarm_gases as sg
from rf_plasma.breakdown_probability import Setup, _time_axis
from rf_plasma.ppt import photons_required
from rf_plasma.rod3d import RodConfig
from rf_plasma.stochastic_seed import Rates
from rf_plasma.streamer2d import ION_K0, SWARM_FILE, SWARM_FILES, Config
from rf_plasma.swarm import COLUMNS, CROSS_SECTION_DIR, read_table


class PowerLawRates:
    """nu_iz = a E0^s, independent of N and f (analytic threshold check)."""

    def __init__(self, a, s):
        self.a, self.s = a, s

    def nu_iz(self, e0, n, f):
        return self.a * np.asarray(e0, float) ** self.s


class TestTables(unittest.TestCase):
    def test_columns_grid_and_attachment(self):
        he = read_table(SWARM_FILE)
        for gas in ("Ar", "N2"):
            d = read_table(SWARM_FILES[gas])
            for c in COLUMNS + sg.EXTRA_COLUMNS:
                self.assertIn(c, d)
            np.testing.assert_allclose(d["E_N_Td"], he["E_N_Td"], rtol=1e-6)
            self.assertTrue(np.all(np.asarray(d["katt_m3_s"]) == 0.0))
            en = np.asarray(d["E_N_Td"])
            sel = en > 20
            w = np.asarray(d["muN_1_Vms"]) * en * 1e-21
            np.testing.assert_allclose(np.asarray(d["w_m_s"])[sel], w[sel], rtol=1e-5)

    def test_rotational_data_and_populations(self):
        p = sg.n2_rot_populations(300.0)
        self.assertAlmostEqual(p.sum(), 1.0, places=12)
        self.assertEqual(int(np.argmax(p)), 6)          # N2 at 300 K: J_max = 6 (ortho weight 6)
        rot = sg.parse_state_resolved(CROSS_SECTION_DIR / "LoKI_N2_rot_LXCat.txt")
        b0 = [b for b in rot if b["lower"] == 0][0]
        self.assertEqual(b0["upper"], 2)
        self.assertAlmostEqual(b0["threshold"], 6 * sg.B_N2_EV, delta=0.05 * 6 * sg.B_N2_EV)
        self.assertLess(sg.n2_vib_populations(300.0)[1], 1e-4)


class TestGasSelection(unittest.TestCase):
    def test_he_defaults_unchanged(self):
        c = Config()
        self.assertEqual(c.gas, "He")
        self.assertEqual(c.swarm_file, str(SWARM_FILE))
        self.assertEqual(c.K0_ion_m2_Vs, 16.7e-4)
        self.assertEqual(RodConfig().gas, "He")
        self.assertEqual(SWARM_FILES["He"], SWARM_FILE)

    def test_gas_switch(self):
        c = Config(gas="N2")
        self.assertEqual(c.swarm_file, str(SWARM_FILES["N2"]))
        self.assertEqual(c.K0_ion_m2_Vs, ION_K0["N2"])
        c2 = Config(gas="Ar", swarm_file="custom.csv", K0_ion_m2_Vs=1e-4)   # explicit values win
        self.assertEqual((c2.swarm_file, c2.K0_ion_m2_Vs), ("custom.csv", 1e-4))
        with self.assertRaises(ValueError):
            Config(gas="Xe")

    def test_he_rates_identical(self):
        a, b = Rates("He"), gc.rates("He")
        st = Setup()
        e = np.geomspace(1e6, 3e7, 25)
        np.testing.assert_array_equal(a.nu_iz(e, st.n_gas, st.f_hz), b.nu_iz(e, st.n_gas, st.f_hz))


class TestThreshold(unittest.TestCase):
    def test_power_law_envelope_analytic(self):
        st = Setup()
        a, s, g_req = 1e-30, 5.0, 21.7
        r = PowerLawRates(a, s)
        t, env = _time_axis(st)
        t_flat = st.t_rf_end_s - st.t_laser_s
        tail = st.tau_amp_s / s * -math.expm1(-s * st.n_tau_after)
        e_exact = (g_req / (a * (t_flat + tail))) ** (1 / s)
        e_num = gc.threshold_envelope(r, st, g_req)
        self.assertAlmostEqual(e_num / e_exact, 1.0, delta=2e-3)   # trapezoid on dt = 0.2 ns
        e_flat = gc.threshold_flat(r, st, g_req, 50e-9)
        self.assertAlmostEqual(e_flat / (g_req / (a * 50e-9)) ** (1 / s), 1.0, places=6)

    def test_gas_ordering_and_pressure_scaling(self):
        g_req = math.log(1e12 / gc.seed_n0(gc.setup(200.0)))
        self.assertAlmostEqual(g_req, 21.70, delta=0.02)
        for psig in (50.0, 200.0):
            st = gc.setup(psig)
            thr = {g: gc.threshold_envelope(gc.rates(g), st, g_req) for g in gc.GASES}
            self.assertGreater(thr["N2"], 4 * thr["He"])
            self.assertGreater(thr["Ar"], thr["He"])      # model result (opposite of the user's observation)
            # collisional regime: threshold E_eff/N nearly pressure independent
            en, _ = gc.eff_en_td(gc.rates("N2"), thr["N2"], st.n_gas, st.f_hz)
            self.assertTrue(90 < en < 125)


class TestValidation(unittest.TestCase):
    def test_against_phelps_and_norman(self):
        v = gc.validation()
        for x in v["N2"]["w_m_s"]:
            self.assertLess(abs(x["dev"]), 0.05)
        for x in v["N2"]["alphaN_m2"]:
            if x["E_N_Td"] >= 150:
                self.assertLess(abs(x["dev"]), 0.40)
        for x in v["Ar"]["w_m_s"]:
            self.assertLess(abs(x["dev"]), 0.20)
        for x in v["Ar"]["alphaN_m2"]:
            if x["E_N_Td"] >= 20:
                self.assertLess(abs(x["dev"]), 0.40)
        for g in ("N2", "Ar"):
            for x in v[g]["mean_energy_eV"]:
                self.assertLess(abs(x["dev"]), 0.05)
        self.assertLess(abs(v["Ar_DC_breakdown_10bar_1cm"]["dev"]), 0.10)

    def test_photon_orders(self):
        st = Setup()
        i = st.i_peak(0.75e-3)
        self.assertEqual([photons_required(i, st.wavelength_m, gc.IP_EV[g]) for g in gc.GASES], [11, 7, 7])


if __name__ == "__main__":
    unittest.main()
