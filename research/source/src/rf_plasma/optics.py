"""Energy-normalized Gaussian-pulse benchmark; no plasma propagation."""

import math

import numpy as np
from scipy.constants import c

from .config import Laser, _finite, _positive


def rayleigh_length_m(laser: Laser) -> float:
    return (
        math.pi * laser.phase_refractive_index * laser.waist_radius_m**2
        / (laser.beam_quality_m2 * laser.wavelength_m)
    )


def beam_radius_m(laser: Laser, axial_distance_m: float | np.ndarray) -> np.ndarray:
    z = np.asarray(axial_distance_m, dtype=float)
    if not np.all(np.isfinite(z)):
        raise ValueError("axial_distance_m must be finite")
    return laser.waist_radius_m * np.sqrt(1 + (z / rayleigh_length_m(laser)) ** 2)


def gaussian_intensity_w_m2(
    laser: Laser,
    radius_m: float | np.ndarray,
    axial_distance_m: float | np.ndarray,
    time_from_focus_s: float | np.ndarray,
) -> np.ndarray:
    radius = np.asarray(radius_m, dtype=float)
    time = np.asarray(time_from_focus_s, dtype=float)
    z = np.asarray(axial_distance_m, dtype=float)
    if not np.all(np.isfinite(radius)) or np.any(radius < 0):
        raise ValueError("radius_m must be finite and nonnegative")
    if not np.all(np.isfinite(time)):
        raise ValueError("time_from_focus_s must be finite")
    width = beam_radius_m(laser, z)
    retarded_time = time - laser.group_refractive_index * z / c
    temporal_normalization = math.sqrt(4 * math.log(2) / math.pi) / laser.intensity_fwhm_s
    temporal = np.exp(-4 * math.log(2) * (retarded_time / laser.intensity_fwhm_s) ** 2)
    radial = np.exp(-2 * (radius / width) ** 2)
    return 2 * laser.energy_at_focus_j * temporal_normalization * temporal * radial / (
        math.pi * width**2
    )


def ideal_focal_length_m(laser: Laser, incident_beam_radius_m: float) -> float:
    _positive("incident_beam_radius_m", incident_beam_radius_m)
    return (
        math.pi * laser.phase_refractive_index * laser.waist_radius_m
        * incident_beam_radius_m / (laser.beam_quality_m2 * laser.wavelength_m)
    )


def energy_before_window_j(laser: Laser, transmission: float) -> float:
    _finite("transmission", transmission)
    if not 0 < transmission <= 1:
        raise ValueError("transmission must be greater than zero and at most one")
    return laser.energy_at_focus_j / transmission
