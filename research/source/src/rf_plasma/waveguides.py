"""Analytic perfect-conductor guide and cavity benchmarks."""

import cmath
import math

from scipy.constants import c, epsilon_0, mu_0
from scipy.special import jn_zeros, jnp_zeros

from .config import Copper, RectangularGeometry, _finite, _positive


def rectangular_cutoff_hz(width_m: float, height_m: float, m: int, n: int) -> float:
    _positive("width_m", width_m)
    _positive("height_m", height_m)
    if not isinstance(m, int) or not isinstance(n, int) or m < 0 or n < 0 or (m == n == 0):
        raise ValueError("TE mode indices must be nonnegative integers and not both zero")
    return c / 2 * math.hypot(m / width_m, n / height_m)


def circular_cutoff_hz(radius_m: float, mode: str, m: int, n: int) -> float:
    _positive("radius_m", radius_m)
    if mode not in {"TE", "TM"}:
        raise ValueError("mode must be TE or TM")
    if not isinstance(m, int) or not isinstance(n, int) or m < 0 or n < 1:
        raise ValueError("circular mode indices require m >= 0 and n >= 1")
    roots = jnp_zeros(m, n) if mode == "TE" else jn_zeros(m, n)
    return c * float(roots[-1]) / (2 * math.pi * radius_m)


def te10_propagation_constant_per_m(frequency_hz: float, width_m: float) -> complex:
    _positive("frequency_hz", frequency_hz)
    _positive("width_m", width_m)
    k = 2 * math.pi * frequency_hz / c
    return cmath.sqrt(k**2 - (math.pi / width_m) ** 2)


def te10_power_for_peak_electric_field_w(
    frequency_hz: float, width_m: float, height_m: float, peak_field_v_m: float
) -> float:
    """Time-averaged power in an ideal, vacuum-filled guide's TE10 mode."""
    _positive("height_m", height_m)
    _finite("peak_field_v_m", peak_field_v_m)
    if peak_field_v_m < 0:
        raise ValueError("peak_field_v_m must be nonnegative")
    beta = te10_propagation_constant_per_m(frequency_hz, width_m)
    if beta.real <= 0 or beta.imag != 0:
        raise ValueError("frequency must exceed the TE10 cutoff")
    impedance = 2 * math.pi * frequency_hz * mu_0 / beta.real
    return width_m * height_m * peak_field_v_m**2 / (4 * impedance)


def te101_frequency_hz(geometry: RectangularGeometry) -> float:
    return c / 2 * math.hypot(1 / geometry.width_m, 1 / geometry.closed_cavity_length_m)


def te101_length_for_frequency_m(width_m: float, frequency_hz: float) -> float:
    _positive("width_m", width_m)
    _positive("frequency_hz", frequency_hz)
    inverse_length_squared = (2 * frequency_hz / c) ** 2 - (1 / width_m) ** 2
    if inverse_length_squared <= 0:
        raise ValueError("frequency must exceed the TE10 cutoff for this cavity width")
    return 1 / math.sqrt(inverse_length_squared)


def skin_depth_m(frequency_hz: float, copper: Copper) -> float:
    _positive("frequency_hz", frequency_hz)
    omega = 2 * math.pi * frequency_hz
    return math.sqrt(2 / (omega * mu_0 * copper.relative_permeability * copper.conductivity_s_m))


def surface_resistance_ohm(frequency_hz: float, copper: Copper) -> float:
    return 1 / (copper.conductivity_s_m * skin_depth_m(frequency_hz, copper))


def te10_copper_power_attenuation_per_m(
    frequency_hz: float, geometry: RectangularGeometry, copper: Copper
) -> float:
    """First-order TE10 power exponent alpha in P(z) = P(0) exp(-alpha*z)."""
    beta = te10_propagation_constant_per_m(frequency_hz, geometry.width_m)
    if beta.real <= 0 or beta.imag != 0:
        raise ValueError("frequency must exceed the TE10 cutoff")
    a, b = geometry.width_m, geometry.height_m
    omega = 2 * math.pi * frequency_hz
    k = omega / c
    kx = math.pi / a
    return (
        2 * surface_resistance_ohm(frequency_hz, copper) * (2 * b * kx**2 + a * k**2)
        / (a * b * beta.real * omega * mu_0)
    )


def isolated_te101_copper_q(geometry: RectangularGeometry, copper: Copper) -> float:
    """First-order good-conductor Q for an ideal closed, vacuum-filled cavity."""
    a = geometry.width_m
    b = geometry.height_m
    d = geometry.closed_cavity_length_m
    omega = 2 * math.pi * te101_frequency_hz(geometry)
    kx = math.pi / a
    kz = math.pi / d
    h_scale = 1 / (omega * mu_0)
    surface_h2 = h_scale**2 * (
        b * d * kx**2 + a * d * (kx**2 + kz**2) / 2 + a * b * kz**2
    )
    stored_energy = epsilon_0 * a * b * d / 8
    wall_power = surface_resistance_ohm(omega / (2 * math.pi), copper) * surface_h2 / 2
    return omega * stored_energy / wall_power
