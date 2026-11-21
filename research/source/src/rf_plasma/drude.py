"""Unmagnetized collisional electron response with exp(-i omega t)."""

import cmath
import math

from scipy.constants import c, elementary_charge, electron_mass, epsilon_0

from .config import _finite, _positive
from .waveguides import te10_propagation_constant_per_m


def conductivity_s_m(density_m3: float, collision_rate_s: float, frequency_hz: float) -> complex:
    _finite("density_m3", density_m3)
    _finite("collision_rate_s", collision_rate_s)
    _positive("frequency_hz", frequency_hz)
    if density_m3 < 0 or collision_rate_s < 0:
        raise ValueError("density and collision rate must be nonnegative")
    omega = 2 * math.pi * frequency_hz
    return density_m3 * elementary_charge**2 / (
        electron_mass * complex(collision_rate_s, -omega)
    )


def relative_permittivity(
    density_m3: float,
    collision_rate_s: float,
    frequency_hz: float,
    background_permittivity: float = 1.0,
) -> complex:
    _positive("background_permittivity", background_permittivity)
    omega = 2 * math.pi * frequency_hz
    sigma = conductivity_s_m(density_m3, collision_rate_s, frequency_hz)
    return background_permittivity + 1j * sigma / (epsilon_0 * omega)


def peak_phasor_heating_w_m3(
    density_m3: float, collision_rate_s: float, frequency_hz: float, field_peak_v_m: float
) -> float:
    _finite("field_peak_v_m", field_peak_v_m)
    sigma = conductivity_s_m(density_m3, collision_rate_s, frequency_hz)
    return 0.5 * sigma.real * field_peak_v_m**2


def te10_slab_scattering(
    frequency_hz: float, width_m: float, thickness_m: float,
    density_m3: float, collision_rate_s: float,
) -> tuple[complex, complex]:
    """Analytic TE10 reflection and transmission for a full-cross-section Drude slab."""
    _positive("thickness_m", thickness_m)
    beta_air = te10_propagation_constant_per_m(frequency_hz, width_m)
    if beta_air.real <= 0 or beta_air.imag != 0:
        raise ValueError("the surrounding guide must propagate TE10")
    eps = relative_permittivity(density_m3, collision_rate_s, frequency_hz)
    k = 2 * math.pi * frequency_hz / c
    beta_slab = cmath.sqrt(k**2 * eps - (math.pi / width_m) ** 2)
    sine_over_beta = thickness_m if abs(beta_slab) < 1e-14 else cmath.sin(beta_slab * thickness_m) / beta_slab
    denominator = (
        2 * cmath.cos(beta_slab * thickness_m)
        - 1j * (beta_air + beta_slab**2 / beta_air) * sine_over_beta
    )
    reflection = -1j * (beta_air - beta_slab**2 / beta_air) * sine_over_beta / denominator
    transmission = 2 / denominator
    return reflection, transmission
