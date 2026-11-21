import unittest

import numpy as np

from rf_plasma import plasma3d as k3
from rf_plasma import plasma3d_cr as kc
from rf_plasma.breakdown_probability import Setup, seed_per_length
from rf_plasma.cavity_dynamics import advance_cells, coarsen_plane, perpendicular_seed, weighted_edge_work
from rf_plasma.gas_dynamics import Euler1D, M_HE, chemical_energy, sod_star
from rf_plasma.gas_model import load_gas_model
from rf_plasma.radial_afterglow import chemistry_step, project_slice
from rf_plasma.radiation_diagnostics import free_free_estimate


class GasDynamicsTests(unittest.TestCase):
    def test_visual_coarsening_preserves_an_odd_index_seed(self):
        values = np.zeros((5, 7))
        values[:, 3] = 2
        coarse = coarsen_plane(values)
        weights = np.minimum(2, 5 - np.arange(0, 5, 2))[:, None] * np.minimum(2, 7 - np.arange(0, 7, 2))
        self.assertAlmostEqual(float((coarse * weights).sum()), float(values.sum()))
        self.assertGreater(coarse.max(), 0)

    def test_cylindrical_uniform_gas_stays_at_rest(self):
        solver = Euler1D(np.linspace(0, 2e-3, 81), cylindrical=True)
        u = np.zeros((3, 80))
        u[0] = 2.0
        u[2] = 1.48e6 / (solver.gamma - 1)
        _, frames = solver.integrate(u, 1e-6, 6)
        np.testing.assert_allclose(frames[:, 0], 2.0, rtol=1e-12)
        self.assertLess(np.max(np.abs(frames[:, 1])), 1e-9)
        np.testing.assert_allclose(frames[:, 2], u[2, 0], rtol=1e-12)

    def test_sod_shock_and_energy_boundary_accounting(self):
        positions = []
        for n in (200, 400):
            solver = Euler1D(np.linspace(0, 1, n + 1), gamma=1.4)
            u = np.zeros((3, n))
            u[0] = np.where(solver.centres < 0.5, 1.0, 0.125)
            u[2] = np.where(solver.centres < 0.5, 1.0, 0.1) / 0.4
            initial = (u * solver.volumes).sum(axis=1)
            _, frames = solver.integrate(u, 0.2, 2)
            end = frames[-1]
            final = (end * solver.volumes).sum(axis=1)
            np.testing.assert_allclose(final[[0, 2]] + solver.boundary_out[[0, 2]],
                                       initial[[0, 2]], rtol=2e-12)
            rho, _, p, _ = solver.primitive(end)
            right = np.flatnonzero(solver.centres > 0.75)
            j = right[np.argmax(-np.diff(p)[right[:-1]])]
            positions.append(solver.centres[j])
            self.assertGreater(rho.min(), 0)
        exact = 0.5 + sod_star()[2] * 0.2
        self.assertLess(abs(positions[-1] - exact), 0.012)
        self.assertLessEqual(abs(positions[-1] - exact), abs(positions[0] - exact) + 0.003)

    def test_continuum_is_zero_without_electrons_and_band_is_bounded(self):
        zero = np.zeros((6, 5))
        total, visible = free_free_estimate(zero)
        np.testing.assert_array_equal(total, 0)
        np.testing.assert_array_equal(visible, 0)
        zero[0] = 1e20
        zero[2] = 1e20 * 3.0
        total, visible = free_free_estimate(zero)
        self.assertTrue(np.all((visible > 0) & (visible < total)))
        zero[:] = 0
        zero[0], zero[2] = 1e-310, 3e-310
        with np.errstate(divide="raise", invalid="raise"):
            total, visible = free_free_estimate(zero)
        self.assertTrue(np.isfinite(total).all() and np.isfinite(visible).all())

    def test_perpendicular_seed_conserves_count_and_saturates(self):
        active = np.ones((61, 5, 5), dtype=bool)
        beam = Setup(n_sat_m3=1e16, n_ref_m3=1e17)
        density, count = perpendicular_seed(active, (30, 1, 2), 1e-3, beam)
        expected = seed_per_length(beam, beam.energy_j, (np.arange(61) - 30) * 1e-3).sum() * 1e-3
        self.assertAlmostEqual(count / expected, 1, places=12)
        np.testing.assert_array_equal(density[:, 1, 2], density[:, 2, 2])
        self.assertLess(np.count_nonzero(density), 3 * active.shape[0])
        unsaturated = Setup(n_ref_m3=1e17)
        self.assertLess(count, perpendicular_seed(active, (30, 1, 2), 1e-3, unsaturated)[1])

    def test_afterglow_source_transfers_heat_without_double_counting(self):
        gm = load_gas_model("He", 532e-9, allow_unverified=True, extended=True, cr=True)
        n0 = 1.48e6 / (k3.KB * 300)
        u = np.zeros((9, 3))
        u[0] = M_HE * n0
        u[3] = 1e19
        u[5] = 2e19
        u[6] = 1e20
        u[2] = 1.5 * k3.KB * n0 * 300 + u[5] * k3.QE + chemical_energy(u[3:], gm.prm)
        total0 = u[2].copy()
        gas0 = u[2] - u[5] * k3.QE - chemical_energy(u[3:], gm.prm)
        sink, heat, photons = chemistry_step(u, 1e-9, gm.prm, gm.tabs, gm.ln_e0, gm.h, gm.ln_i0, gm.hi, gm.ppt_tab)
        np.testing.assert_allclose(u[2] + sink, total0, rtol=1e-13)
        gas1 = u[2] - u[5] * k3.QE - chemical_energy(u[3:], gm.prm)
        np.testing.assert_allclose(gas1 - gas0, heat, rtol=2e-7, atol=1e-8)
        self.assertTrue(np.all((photons >= 0) & (photons <= sink)))

    def test_edge_work_reaches_only_occupied_nodes_and_the_kinetic_ledger(self):
        gm = load_gas_model("He", 532e-9, allow_unverified=True, extended=True, cr=True)
        shape = (3, 2, 2)
        weights = np.zeros(shape)
        weights[:, 0, 0] = (2, 6, 0)
        px, py, pz = np.zeros((2, 2, 2)), np.zeros((3, 1, 2)), np.zeros((3, 2, 1))
        px[:, 0, 0] = (8, -3)
        py[0, 0, 0], pz[1, 0, 0] = 7, -2
        work = np.zeros(shape)
        weighted_edge_work(px, py, pz, weights, work)
        np.testing.assert_array_equal(work[:, 0, 0], (9, 1, 0))
        self.assertEqual(float(work[weights == 0].sum()), 0)
        self.assertEqual(float(work.sum()), float(px.sum() + py.sum() + pz.sum()))
        n0 = 1.48e6 / (k3.KB * 300)
        state = np.zeros((6,) + shape)
        state[0] = weights * 1e19
        state[2] = state[0] * 0.5
        tg = np.full(shape, 300.0)
        gas = np.full(shape, 1.5 * k3.KB * n0 * 300)
        wp2, nu, ledger = np.zeros(shape), np.zeros(shape), np.zeros((3, 11))
        advance_cells(np.argwhere(weights > 0), state, tg, gas, work, 1e-12, n0, gm.prm,
                      gm.tabs, gm.ln_e0, gm.h, gm.ln_i0, gm.hi, gm.ppt_tab, wp2, nu, ledger)
        self.assertAlmostEqual(float(ledger[:, 0].sum()), float(work.sum()), places=12)
        np.testing.assert_allclose(tg[weights > 0],
                                   (gas / (1.5 * k3.KB * (n0-state[1]-state[4])))[weights > 0])

    def test_afterglow_passes_ground_neutrals_not_total_nuclei(self):
        gm = load_gas_model("He", 532e-9, allow_unverified=True, extended=True, cr=True)
        n0 = 1.48e6 / (k3.KB * 300)
        u = np.zeros((9, 1))
        u[0] = M_HE * n0
        u[3], u[4], u[6] = 0.4*n0, 0.1*n0, 0.25*n0
        u[5] = (u[3] + u[4]) * 2
        u[2] = 1.5*k3.KB*(n0-u[4])*300 + u[5]*k3.QE + chemical_energy(u[3:], gm.prm)
        solver = Euler1D([0, 1, 2], plasma_prm=gm.prm)
        for _ in range(3):
            expected = u[3:, 0].copy()
            neutral = n0-expected[0]-2*expected[1]-expected[3]-2*expected[4]-expected[5]
            gas = u[2, 0]-expected[2]*k3.QE-chemical_energy(expected, gm.prm)
            tg = gas/(1.5*k3.KB*(n0-expected[1]-expected[4]))
            kc.kinetic_cell_cr(expected, 1e-12, neutral, tg, 0.0, 0.0, gm.prm, gm.tabs,
                               gm.ln_e0, gm.h, gm.ln_i0, gm.hi, gm.ppt_tab, np.zeros(11))
            chemistry_step(u, 1e-12, gm.prm, gm.tabs, gm.ln_e0, gm.h, gm.ln_i0, gm.hi, gm.ppt_tab)
            np.testing.assert_allclose(u[3:, 0], expected, rtol=2e-13)
            solver.primitive(np.repeat(u, 2, axis=1))
        u[6] = n0
        with self.assertRaisesRegex(ArithmeticError, "nuclei"):
            solver.primitive(np.repeat(u, 2, axis=1))

    def test_handoff_preserves_composition_consistent_gas_temperature_and_energy(self):
        gm = load_gas_model("He", 532e-9, allow_unverified=True, extended=True, cr=True)
        n0 = 1.48e6 / (k3.KB * 300)
        state = np.zeros((6, 3, 3))
        state[1, :2, 1] = 1e22
        gas = 1.5 * k3.KB * (n0-state[1]-state[4]) * 300
        solver, u = project_slice(dict(dx=1e-3, n0=n0, focus=(0, 0, 1),
                                        state_yz=state, gas_energy_yz=gas), 80, 2e-3, gm.prm)
        _, _, _, recovered = solver.primitive(u)
        np.testing.assert_allclose(recovered/(1.5*k3.KB*(n0-u[4]-u[7])), 300, rtol=2e-13)

    def test_electron_compression_work_repartitions_conserved_total_energy(self):
        gm = load_gas_model("He", 532e-9, allow_unverified=True, extended=True, cr=True)
        solver = Euler1D(np.linspace(0, 0.1, 101), plasma_prm=gm.prm)
        n0 = 1e23
        u = np.zeros((9, 100))
        u[0] = M_HE*n0
        u[1] = -100*solver.centres*u[0]
        u[3], u[5] = 1e18, 2e18
        u[2] = 1.5*k3.KB*n0*300 + u[5]*k3.QE + chemical_energy(u[3:], gm.prm) + 0.5*u[1]**2/u[0]
        initial = (u[:3]*solver.volumes).sum(axis=1)
        result = solver.step(u, 1e-8)
        compression = result[0, 50]/u[0, 50]
        self.assertGreater(compression, 1)
        np.testing.assert_allclose(result[5, 50]/u[5, 50], compression**solver.gamma, rtol=1e-10)
        np.testing.assert_allclose((result[:3]*solver.volumes).sum(axis=1)+solver.boundary_out,
                                   initial, rtol=2e-12)

    def test_kinetic_substep_exhaustion_is_explicit(self):
        gm = load_gas_model("He", 532e-9, allow_unverified=True, extended=True, cr=True)
        with self.assertRaisesRegex(ArithmeticError, "exhausted"):
            kc.kinetic_cell_cr(np.zeros(6), 1e-9, 1e26, 300.0, 0.0, 0.0, gm.prm,
                               gm.tabs, gm.ln_e0, gm.h, gm.ln_i0, gm.hi, gm.ppt_tab,
                               np.zeros(11), max_substeps=0)

    def test_unassigned_drude_memory_work_is_not_silently_dropped(self):
        with self.assertRaisesRegex(ArithmeticError, "Drude-memory"):
            weighted_edge_work(np.ones((1, 1, 1)), np.zeros((2, 0, 1)), np.zeros((2, 1, 0)),
                               np.zeros((2, 1, 1)), np.zeros((2, 1, 1)))


if __name__ == "__main__":
    unittest.main()
