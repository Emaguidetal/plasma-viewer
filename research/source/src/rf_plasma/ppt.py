"""Perelomov-Popov-Terent'ev (PPT) photoionization rate with ADK (Ammosov) coefficients.

Linear polarization, cycle-averaged rate, atomic units internally. The formula follows
the form summarized by Couairon & Mysyrowicz, Phys. Rep. 441, 47 (2007), Sec. 2.4 and
Perelomov, Popov & Terent'ev, Sov. Phys. JETP 23, 924 (1966). Molecules are treated as
effective atoms with charge Z_eff (Talebpour et al., Opt. Commun. 163, 29 (1999)).
"""

import math

import numpy as np
from scipy.constants import c, epsilon_0, hbar, h, elementary_charge as QE
from scipy.special import dawsn, gamma, factorial

AU_FIELD = 5.14220674763e11       # V/m
AU_TIME = 2.4188843265857e-17     # s
AU_ENERGY_EV = 27.211386245988


def keldysh_gamma(intensity_w_m2: float, wavelength_m: float, ip_ev: float) -> float:
    field = math.sqrt(2 * intensity_w_m2 / (c * epsilon_0)) / AU_FIELD
    omega = (h * c / wavelength_m / QE) / AU_ENERGY_EV
    return omega * math.sqrt(2 * ip_ev / AU_ENERGY_EV) / field


def _cnl2(n_star: float, l_star: float) -> float:
    return 2 ** (2 * n_star) / (n_star * gamma(n_star + l_star + 1) * gamma(n_star - l_star))


def _flm(l: int, m: int) -> float:
    m = abs(m)
    return (2 * l + 1) * factorial(l + m) / (2**m * factorial(m) * factorial(l - m))


def ppt_rate_s(intensity_w_m2: float, wavelength_m: float, ip_ev: float, z_eff: float = 1.0,
               l: int = 0, m: int = 0, max_terms: int = 4000) -> float:
    """Cycle-averaged PPT ionization rate per neutral (s^-1), linear polarization."""
    if intensity_w_m2 <= 0:
        return 0.0
    if m != 0:
        raise NotImplementedError("only m = 0 (Dawson-function) channel implemented")
    ip = ip_ev / AU_ENERGY_EV
    f = math.sqrt(2 * intensity_w_m2 / (c * epsilon_0)) / AU_FIELD
    omega = (h * c / wavelength_m / QE) / AU_ENERGY_EV
    n_star = z_eff / math.sqrt(2 * ip)
    l_star = n_star - 1
    f0 = (2 * ip) ** 1.5
    gam = omega * math.sqrt(2 * ip) / f
    sq = math.sqrt(1 + gam**2)
    g = 1.5 / gam * ((1 + 1 / (2 * gam**2)) * math.asinh(gam) - sq / (2 * gam))
    alpha = 2 * (math.asinh(gam) - gam / sq)
    beta = 2 * gam / sq
    up = f**2 / (4 * omega**2)
    nu = (ip + up) / omega
    kappa0 = math.ceil(nu)
    total = 0.0
    max_terms = int(min(2e6, max(max_terms, 40.0 / max(alpha, 1e-12))))
    for kk in range(max_terms):
        x = kappa0 + kk - nu
        term = math.exp(-alpha * x) * dawsn(math.sqrt(beta * x))
        total += term
        if kk > 5 and term < 1e-14 * total:
            break
    a_m = 4 / math.sqrt(3 * math.pi) * gam**2 / (1 + gam**2) * total
    log_rate = (math.log(ip * _cnl2(n_star, l_star) * _flm(l, m) * math.sqrt(6 / math.pi) * a_m)
                + (2 * n_star - 1.5) * math.log(2 * f0 / (f * sq)) - 2 * f0 * g / (3 * f))
    return math.exp(log_rate) / AU_TIME


def adk_rate_s(intensity_w_m2: float, ip_ev: float, z_eff: float = 1.0, l: int = 0, m: int = 0) -> float:
    """Cycle-averaged ADK tunnelling rate (quasi-static limit of PPT for gamma -> 0)."""
    ip = ip_ev / AU_ENERGY_EV
    f = math.sqrt(2 * intensity_w_m2 / (c * epsilon_0)) / AU_FIELD
    n_star = z_eff / math.sqrt(2 * ip)
    f0 = (2 * ip) ** 1.5
    rate = (ip * _cnl2(n_star, n_star - 1) * _flm(l, m) * math.sqrt(6 / math.pi)
            * (2 * f0 / f) ** (2 * n_star - abs(m) - 1.5) * math.exp(-2 * f0 / (3 * f)))
    return rate / AU_TIME


def photons_required(intensity_w_m2: float, wavelength_m: float, ip_ev: float) -> int:
    photon_ev = h * c / wavelength_m / QE
    omega = 2 * math.pi * c / wavelength_m
    up_ev = QE * (2 * intensity_w_m2 / (c * epsilon_0)) / (4 * 9.1093837015e-31 * omega**2)
    return int(math.ceil((ip_ev + up_ev) / photon_ev))


def ppt_table(wavelength_m: float, ip_ev: float, z_eff: float, l: int = 0,
              i_min_w_m2: float = 1e12, i_max_w_m2: float = 1e19, points: int = 281):
    """Log-uniform table of ln(rate) versus ln(intensity) for the kinetics kernels."""
    ln_i = np.linspace(math.log(i_min_w_m2), math.log(i_max_w_m2), points)
    rates = np.array([ppt_rate_s(math.exp(x), wavelength_m, ip_ev, z_eff, l) for x in ln_i])
    return ln_i, np.log(np.maximum(rates, 1e-300))
