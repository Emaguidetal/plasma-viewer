"""Restricted atomic-He demonstrations; NOT dense-helium sheet chemistry.

States: neutral ground, 2^1S, 2^1P, atomic He+, atomic He2+ (bare nucleus).
Includes excitation/superelastic pairs, sequential ionization and a nonlocal,
quasi-static plane-parallel resonance redistribution operator. No fitted
high-pressure radiation profile or molecular-ion afterglow is claimed.
"""

import argparse
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np
from bolos import parser as collision_parser
from numba import njit
from scipy.constants import c, elementary_charge as QE, electron_mass as ME, k as KB, atomic_mass
from scipy.integrate import solve_ivp
from scipy.special import voigt_profile

from .swarm import CROSS_SECTION_DIR

E_S = 166277.4403 / 8065.543937
E_P = 171134.8970 / 8065.543937
I_0 = 24.587389011
I_1 = 54.4177655282
A_PG = 1.7989e9
A_PS = 1.97e6
LAMBDA = 58.433436e-9
G_G, G_S, G_P = 1, 1, 3
MASS = 4.002602 * atomic_mass
STATE_ENERGIES = np.array([0.0, E_S, E_P, I_0, I_0 + I_1])
SOURCES = {
    "levels_A_wavelength": "https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable3_a.htm",
    "levels": "https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable5.htm",
    "ionization_energies": "https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=He+I%3BHe+II&units=1&format=0",
    "excitation": "cases/cross_sections/LoKI_He_LXCat.txt (resolved He(2S1), He(2P1), IST-Lisbon/LXCat)",
    "superelastic": "Klein-Rosseland detailed balance; LoKI Collision.m superElasticCrossSection",
    "ionization_fit": "Voronov 1997, doi:10.1006/adnd.1997.0732; https://www.pa.uky.edu/~verner/dima/col/cfit.dat",
    "line_absorption": "Einstein relations, https://physics.nist.gov/Pubs/AtSpec/latex.html, equations 16 and 19",
    "radiation_reduction": "Derived Beer-Lambert single-flight probabilities, complete redistribution, plane-parallel angular quadrature",
}


def ionization_rates(te_ev):
    """Voronov ground-target fits, SI; do not extrapolate below fitted ranges."""
    te = np.asarray(te_ev, dtype=float)
    if np.any(~np.isfinite(te)) or np.any(te < 3) or np.any(te > 2e4):
        raise ValueError("both helium ionization fits require 3 <= Te <= 20000 eV")
    rates = []
    for threshold, amplitude, pp, xx, kk in (
        (24.6, 1.75e-14, 0.0, 0.1800, 0.35),
        (54.4, 2.05e-15, 1.0, 0.2650, 0.25),
    ):
        u = threshold / te
        rates.append(amplitude * (1 + pp * np.sqrt(u)) / (xx + u) * u**kk * np.exp(-u))
    return tuple(rates)


@lru_cache(maxsize=1)
def excitation_tables():
    with (CROSS_SECTION_DIR / "LoKI_He_LXCat.txt").open(encoding="utf-8") as fh:
        processes = collision_parser.parse(fh)
    te = np.geomspace(0.03, 100.0, 600)
    xx, ww = np.polynomial.laguerre.laggauss(96)
    down = []
    shifts = {}
    for name, energy, weight in (("He(2S1)", E_S, G_S), ("He(2P1)", E_P, G_P)):
        matches = [p for p in processes if p["kind"] == "EXCITATION"
                   and p["target"] == "He" and p["product"] == name]
        if len(matches) != 1:
            raise ValueError(f"expected exactly one resolved {name} excitation channel")
        process = matches[0]
        data = np.asarray(process["data"])
        shift = energy - process["threshold"]
        shifts[name] = shift
        incident = energy + te[:, None] * xx
        sigma = np.interp(incident, data[:, 0] + shift, data[:, 1], left=0, right=0)
        # Directly integrate the inverse rate without exp(+Delta/Te) overflow.
        kd = math.sqrt(8 * QE / (math.pi * ME)) / np.sqrt(te) / weight * np.sum(
            ww * incident * sigma, axis=1)
        down.append(kd)
    return te, np.asarray(down), shifts


MANIFOLD = (("He(2S3)", 19.819614525, 3), ("He(2S1)", 20.615774823, 1),
            ("He(2P3)", 20.96408689, 9), ("He(2P1)", 21.2180227112, 3))
G_MANIFOLD = sum(g for _, _, g in MANIFOLD)
E_MANIFOLD = sum(e * g for _, e, g in MANIFOLD) / G_MANIFOLD


def manifold_tables(eps_grid):
    """Maxwellian (Te = 2 eps / 3) rate tables for the lumped n=2 He* manifold.

    Populations inside the manifold are taken as statistical (g_i / G). Returns
      k_sup:  superelastic e + He* -> e + He per He* atom (m^3/s), sum_i g_i kdown_i / G, with kdown_i
              integrated directly from the LoKI/LXCat excitation cross sections (Klein-Rosseland)
      k_step: e + He* -> He+ + 2e per He* atom (m^3/s), Thomson cross section
              sigma = pi (e^2 / 4 pi eps0)^2 (1/E)(1/I - 1/E) (Lieberman & Lichtenberg 2005, sec. 3.4),
              I_i = Ip - E_i; an order-of-magnitude classical estimate for excited-state ionization.
    """
    with (CROSS_SECTION_DIR / "LoKI_He_LXCat.txt").open(encoding="utf-8") as fh:
        processes = collision_parser.parse(fh)
    te = np.asarray(eps_grid, dtype=float) / 1.5
    xx, ww = np.polynomial.laguerre.laggauss(96)
    pref = math.sqrt(8 * QE / (math.pi * ME))
    k_sup = np.zeros_like(te)
    k_step = np.zeros_like(te)
    e2 = QE**2 / (4 * math.pi * 8.8541878128e-12)  # J m
    for name, energy, weight in MANIFOLD:
        matches = [p for p in processes if p["kind"] == "EXCITATION"
                   and p["target"] == "He" and p["product"] == name]
        if len(matches) != 1:
            raise ValueError(f"expected exactly one resolved {name} excitation channel")
        data = np.asarray(matches[0]["data"])
        shift = energy - matches[0]["threshold"]
        incident = energy + te[:, None] * xx
        sigma = np.interp(incident, data[:, 0] + shift, data[:, 1], left=0, right=0)
        kd = pref / np.sqrt(te) / weight * np.sum(ww * incident * sigma, axis=1)
        k_sup += weight * kd / G_MANIFOLD
        ion = I_0 - energy
        ee = ion + te[:, None] * xx  # eV, Maxwellian above threshold: k = pref Te^-1/2 e^{-I/Te} int x sigma
        sig_t = math.pi * (e2 / QE)**2 * (1 / ee) * (1 / ion - 1 / ee)  # m^2 with energies in eV
        kst = pref / np.sqrt(te) * np.exp(-ion / te) * np.sum(ww * ee * sig_t, axis=1)
        k_step += weight * kst / G_MANIFOLD
    return k_sup, k_step


def excitation_rates(te_ev):
    te = np.asarray(te_ev, dtype=float)
    grid, down_table, _ = excitation_tables()
    if np.any(~np.isfinite(te)) or np.any(te < grid[0]) or np.any(te > grid[-1]):
        raise ValueError("resolved excitation table requires 0.03 <= Te <= 100 eV")
    down = np.array([np.exp(np.interp(np.log(te), np.log(grid), np.log(row))) for row in down_table])
    weights = np.array([G_S, G_P]).reshape((2,) + (1,) * te.ndim)
    energy = np.array([E_S, E_P]).reshape((2,) + (1,) * te.ndim)
    up = down * weights * np.exp(-energy / te)
    return up, down


@njit(cache=True)
def _redistribution(ncells, tau, weights):
    off = np.zeros((ncells, ncells))
    escape = np.zeros(ncells)
    for q in range(len(tau)):
        tt = tau[q]
        survival = -math.expm1(-tt) / tt if tt > 1e-10 else 1 - 0.5 * tt
        absorbed = -math.expm1(-tt)
        for i in range(ncells):
            escape[i] += 0.5 * weights[q] * survival * (
                math.exp(-tt * i) + math.exp(-tt * (ncells - 1 - i)))
            for j in range(ncells):
                if i != j:
                    off[i, j] += 0.5 * weights[q] * survival * absorbed * math.exp(
                        -tt * (abs(j - i) - 1))
    return off, escape


def resonance_transport(ncells, dx_m, ground_m3, tg_k=300.0):
    """Nonlocal single-flight operator with natural+Doppler Voigt redistribution.

    Fixed ground-state opacity; infinite transverse extent; black optical ends.
    Local reabsorption cancels emission analytically. Photon flight times,
    collisional broadening, stimulated emission and lateral escape are omitted.
    """
    if ncells < 1 or dx_m <= 0 or ground_m3 < 0 or tg_k <= 0:
        raise ValueError("invalid resonance transport domain")
    doppler = c / LAMBDA / c * math.sqrt(2 * KB * tg_k / MASS)
    gamma = (A_PG + A_PS) / (4 * math.pi) / doppler
    bounds = np.r_[0.0, np.geomspace(1e-4, 1e8, 200)]
    xg, wg = np.polynomial.legendre.leggauss(3)
    x = (0.5 * (bounds[1:] + bounds[:-1]))[:, None] + np.diff(bounds)[:, None] * xg / 2
    phi = voigt_profile(x, 1 / math.sqrt(2), gamma)
    fw = (np.diff(bounds)[:, None] * wg * phi).ravel()
    integral = float(fw.sum())
    if abs(integral - 1) > 1e-4:
        raise ArithmeticError(f"line-profile quadrature failed normalization: {integral}")
    fw /= integral
    cross = (G_P / G_G) * LAMBDA**2 * A_PG / (8 * math.pi) * phi.ravel() / doppler
    mu, aw = np.polynomial.legendre.leggauss(4)
    mu, aw = (mu + 1) / 2, aw / 2
    tau = (ground_m3 * dx_m * cross[:, None] / mu).ravel()
    weights = (fw[:, None] * aw).ravel()
    off, escape = _redistribution(ncells, tau, weights)
    outgoing = off.sum(axis=1) + escape
    return off, escape, outgoing


def radiation_rhs(excited, off, escape, outgoing):
    return A_PG * (off.T @ excited - outgoing * excited), E_P * A_PG * escape * excited


def local_collisions(pop, electron_energy, density):
    """Fractions and electron energy in eV per nucleus; returns conservative derivatives."""
    ne_fraction = pop[3] + 2 * pop[4]
    te = electron_energy / (1.5 * ne_fraction)
    up, down = excitation_rates(te)
    k0, k1 = ionization_rates(te)
    ne = density * ne_fraction
    sx = ne * (up[0] * pop[0] - down[0] * pop[1])
    px = ne * (up[1] * pop[0] - down[1] * pop[2])
    iz0, iz1 = ne * k0 * pop[0], ne * k1 * pop[3]
    dp = np.array([-sx - px - iz0, sx, px, iz0 - iz1, iz1])
    du = -E_S * sx - E_P * px - I_0 * iz0 - I_1 * iz1
    return dp, du


def demonstration(output: Path):
    """Restricted low-density ionizing-phase example, separate from the 1 MPa RF test."""
    output.mkdir(parents=True, exist_ok=True)
    cells, length, density = 15, 0.02, 1e20
    z = (np.arange(cells) + 0.5) * length / cells
    off, escape, outgoing = resonance_transport(cells, length / cells, density)
    state = np.zeros((7, cells))
    hump = np.exp(-((z - length / 2) / 0.003)**2)
    state[1] = 0.01 * hump
    state[2] = 0.02 * hump
    state[3] = 0.10
    state[0] = 1 - state[:5].sum(axis=0)
    state[5] = 1.5 * state[3] * 80.0
    initial = state.copy()

    def rhs(t, flat):
        s = flat.reshape(7, cells)
        ne = s[3] + 2 * s[4]
        te = s[5] / (1.5 * ne)
        # BDF may evaluate just beyond an event. This extension is not returned.
        evaluation_energy = 1.5 * ne * np.clip(te, 3.0, 100.0)
        dp, du = local_collisions(s[:5], evaluation_energy, density)
        rad, escaped = radiation_rhs(s[2], off, escape, outgoing)
        branch = A_PS * s[2]
        dp[0] -= rad
        dp[1] += branch
        dp[2] += rad - branch
        return np.vstack((dp, du, escaped + (E_P - E_S) * branch)).ravel()

    def low_temperature(t, flat):
        s = flat.reshape(7, cells)
        return float(np.min(s[5] / (1.5 * (s[3] + 2 * s[4]))) - 3.0)
    low_temperature.terminal = True
    low_temperature.direction = -1
    sol = solve_ivp(rhs, (0, 2e-5), state.ravel(), method="BDF", rtol=2e-7, atol=1e-12,
                    events=low_temperature, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    t = np.linspace(0, sol.t[-1], 350)
    hist = sol.sol(t).reshape(7, cells, -1)
    if hist[:5].min() < -1e-8:
        raise ArithmeticError("negative collisional-radiative population")
    total = np.einsum("i,ijt->jt", STATE_ENERGIES, hist[:5]) + hist[5] + hist[6]
    baseline = STATE_ENERGIES @ initial[:5] + initial[5]
    residual = float(np.max(np.abs(total.mean(axis=0) - baseline.mean())) / baseline.mean())
    if residual > 2e-6:
        raise ArithmeticError(f"atomic/electron/radiation energy residual {residual}")
    np.savez_compressed(output / "restricted_cr.npz", t=t, z=z, states=hist, density_m3=density,
                        transport_off=off, escape=escape, energy_residual=residual)
    summary = dict(
        status="restricted_low_density_atomic_demonstration_NOT_1MPa_sheet_chemistry",
        density_m3=density, equivalent_initial_pressure_pa=density * KB * 300,
        end_time_s=float(t[-1]), stopped_at_temperature_limit=bool(len(sol.t_events[0])),
        energy_residual_rel=residual, max_bare_He2plus_fraction=float(hist[4].max()),
        max_excited_singlet_fraction=float(hist[1:3].max()), sources=SOURCES,
        cross_section_threshold_shifts_eV=excitation_tables()[2],
        limitations=[
            "Not coupled to the primary prescribed-sheet RF animation.",
            "Atomic ground, two singlets, He+ and bare He2+ only; not a complete helium mechanism.",
            "Missing excited-state ionization, molecular/excimer chemistry, recombination and triplets.",
            "Natural+Doppler profile is NOT a verified high-pressure-broadened profile.",
            "Fixed-opacity, quasi-static plane-parallel transport; no lateral escape or photon flight-time solution.",
            "No elastic/gas heating in this isolated electron/atomic/radiation energy-exchange example.",
            "Maxwellian electrons; integration terminates at the lower 3 eV fit validity boundary.",
        ])
    (output / "restricted_cr_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    print(json.dumps(demonstration(args.output), indent=2))
