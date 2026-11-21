"""Focused tests for the axisymmetric drift-diffusion-Poisson channel model (rf_plasma.streamer2d).

Run from simulations/rf_plasma with PYTHONPATH=src:
    python -m unittest tests.test_streamer2d
"""

import math
import unittest
from dataclasses import replace

import numpy as np

from rf_plasma import streamer2d as s2


def small_cfg(**kw):
    """Small grid, same physics/numerics as production."""
    base = s2.Config(name="unit", pressure_Pa=0.2e6, E0_V_m=2.23e6, r_max_m=0.4e-3, dr_min_m=2e-6,
                     r_uniform_m=20e-6, dr_growth=1.1, dr_mid_m=10e-6, r_mid_m=0.2e-3, dr_growth_out=1.2,
                     z_half_m=0.6e-3, dz_m=10e-6, L_m=0.2e-3, a_m=15e-6)
    return replace(base, **kw)


class TestPoisson(unittest.TestCase):
    """Axisymmetric FV Laplacian vs an analytic solution with the model's boundary conditions."""

    @staticmethod
    def solve_error(nr, nz, R=1e-3, Z=1e-3, stretch=1.04):
        # nonuniform radial grid (geometric stretching, as in production) and uniform z
        dr0 = R * (stretch - 1) / (stretch ** nr - 1)
        rf = np.concatenate([[0.0], np.cumsum(dr0 * stretch ** np.arange(nr))])
        rf[-1] = R
        zf = np.linspace(-Z, Z, nz + 1)
        g = s2.Grid(rf, zf)
        kr, kz = math.pi / (2 * R), math.pi / (2 * Z)
        r, z = g.rc[:, None], g.zc[None, :]
        phi_exact = np.cos(kr * r) * np.cos(kz * z)          # phi = 0 at r = R and z = +-Z, d/dr = 0 on axis
        rho = s2.EPS0 * ((kr ** 2 + kz ** 2) * np.cos(kr * r) + kr * np.sin(kr * r) / r) * np.cos(kz * z)
        A = s2.assemble(g, s2.EPS0 * g.cr, s2.EPS0 * g.cz)
        phi = s2.factor(A).solve((rho * g.V).ravel()).reshape(g.shape)
        return np.max(np.abs(phi - phi_exact))

    def test_manufactured_solution_second_order(self):
        e1 = self.solve_error(40, 40, stretch=1.04)
        e2 = self.solve_error(80, 80, stretch=1.02)
        self.assertLess(e1, 2e-3)
        self.assertLess(e2, 6e-4)
        self.assertGreater(e1 / e2, 3.0)   # ~second order (4x) on refinement

    def test_uniform_field_of_charged_cylinder(self):
        """Long uniformly charged cylinder: E_r(r) = rho a^2 / (2 eps0 r) outside, rho r/(2 eps0) inside."""
        R, Z, a, rho0 = 2e-3, 20e-3, 0.2e-3, 1e-3
        rf = np.linspace(0, R, 201)
        g = s2.Grid(rf, np.linspace(-Z, Z, 401))
        rho = np.where(g.rc[:, None] < a, rho0, 0.0) * np.ones((1, g.nz))
        A = s2.assemble(g, s2.EPS0 * g.cr, s2.EPS0 * g.cz)
        phi = s2.factor(A).solve((rho * g.V).ravel()).reshape(g.shape)
        Er, _ = s2.face_fields(g, phi, 0.0)
        j = g.nz // 2
        r = g.rf[1:-1]
        ex = np.where(r < a, rho0 * r / (2 * s2.EPS0), rho0 * a ** 2 / (2 * s2.EPS0 * r))
        rel = np.max(np.abs(Er[1:-1, j] - ex)) / ex.max()
        self.assertLess(rel, 2e-3)


class TestChargeConservation(unittest.TestCase):
    def test_total_charge_plus_outflow_constant(self):
        cfg = small_cfg(n0_m3=1e22, n_bg_m3=1e15)
        m = s2.Streamer2D(cfg)
        Qe = s2.E_CH * float(np.sum(m.ne * m.g.V))
        resid = []
        while m.t < 0.3e-9:
            m.step()
            resid.append(abs(m.charge_residual()))
        # space charge must actually develop and particles must actually leave (non-trivial test)
        self.assertGreater(abs(m.total_charge()), 1e-6 * Qe)
        self.assertGreater(abs(m.out_e) + abs(m.out_i), 0.0)
        self.assertLess(max(resid), 1e-11 * Qe)
        self.assertGreater(m.ne.min(), 0.0)


class TestFreeDiffusion(unittest.TestCase):
    def test_gaussian_spreads_at_analytic_rate(self):
        """E0 = 0, no ionization/recombination, tenuous Gaussian (Debye length >> width): free diffusion.

        n(r,z,t) = n0 a^2/(a^2+4Dt) * L/sqrt(L^2+4Dt) * exp(-r^2/(a^2+4Dt) - z^2/(L^2+4Dt)).
        """
        cfg = small_cfg(E0_V_m=0.0, ionization=False, alpha_rec_m3_s=0.0, n0_m3=1e10, n_bg_m3=0.0,
                        a_m=20e-6, L_m=40e-6, r_max_m=0.5e-3, z_half_m=0.4e-3, dz_m=4e-6,
                        dr_min_m=2e-6, r_uniform_m=60e-6, steps_per_rf_period=1)
        m = s2.Streamer2D(cfg)
        D = float(m.D_e.max())
        self.assertAlmostEqual(D, float(m.D_e.min()), delta=1e-6 * D)   # field-free: constant LFA D_e
        t_end = 2 * cfg.a_m ** 2 / (4 * D)                               # width^2 grows threefold in r
        while m.t < t_end * (1 - 1e-12):
            m.step(dt_max=t_end - m.t)
        g = m.g
        wr2 = cfg.a_m ** 2 + 4 * D * m.t
        wz2 = cfg.L_m ** 2 + 4 * D * m.t
        Ntot = float(np.sum(m.ne * g.V))
        r2 = float(np.sum(m.ne * g.V * g.rc[:, None] ** 2)) / Ntot
        z2 = float(np.sum(m.ne * g.V * g.zc[None, :] ** 2)) / Ntot
        self.assertAlmostEqual(r2 / wr2, 1.0, delta=0.02)          # <r^2> = a^2 + 4Dt
        self.assertAlmostEqual(z2 / (wz2 / 2), 1.0, delta=0.02)    # <z^2> = (L^2 + 4Dt)/2
        peak = cfg.n0_m3 * cfg.a_m ** 2 / wr2 * cfg.L_m / math.sqrt(wz2)
        self.assertAlmostEqual(float(m.ne.max()) / peak, 1.0, delta=0.03)
        n0_count = cfg.n0_m3 * math.pi ** 1.5 * cfg.a_m ** 2 * cfg.L_m
        self.assertAlmostEqual(Ntot / n0_count, 1.0, delta=0.01)


class TestSwarmClosure(unittest.TestCase):
    def test_table_interpolation_and_einstein(self):
        N = 0.2e6 / (s2.KB * 300.0)
        sw = s2.Swarm(s2.SWARM_FILE, N)
        d = s2.read_swarm(s2.SWARM_FILE)
        k = 50
        mu, D, nu = sw.at(np.array([d["E_N_Td"][k]]))
        self.assertAlmostEqual(mu[0] * N / d["muN_1_Vms"][k], 1.0, delta=1e-3)
        self.assertAlmostEqual(D[0] / (mu[0] * 2 / 3 * d["mean_energy_eV"][k]), 1.0, delta=2e-3)
        self.assertAlmostEqual(nu[0] / (N * d["kiz_m3_s"][k]), 1.0, delta=2e-3)
        self.assertEqual(float(sw.at(np.array([0.5]))[2][0]), 0.0)   # noise branch removed


if __name__ == "__main__":
    unittest.main()
