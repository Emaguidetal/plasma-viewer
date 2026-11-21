"""Validated reference-case inputs, expressed in SI units."""

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any

from scipy.constants import Boltzmann

PSI_TO_PA = 6894.757293168


def _finite(name: str, value: float) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")


def _positive(name: str, value: float) -> None:
    _finite(name, value)
    if value <= 0:
        raise ValueError(f"{name} must be positive")


def gauge_psi_to_absolute_pa(gauge_psi: float, ambient_pa: float) -> float:
    _finite("gauge_psi", gauge_psi)
    _positive("ambient_pa", ambient_pa)
    absolute_pa = gauge_psi * PSI_TO_PA + ambient_pa
    _positive("absolute pressure", absolute_pa)
    return absolute_pa


@dataclass(frozen=True)
class Gas:
    species: str
    gauge_pressure_pa: float
    ambient_pressure_pa: float
    temperature_k: float

    def __post_init__(self) -> None:
        if self.species not in {"He", "N2", "Ar"}:
            raise ValueError("species must be He, N2, or Ar")
        _finite("gauge_pressure_pa", self.gauge_pressure_pa)
        _positive("ambient_pressure_pa", self.ambient_pressure_pa)
        _positive("temperature_k", self.temperature_k)
        _positive("absolute pressure", self.absolute_pressure_pa)

    @property
    def absolute_pressure_pa(self) -> float:
        return self.gauge_pressure_pa + self.ambient_pressure_pa

    @property
    def ideal_neutral_density_m3(self) -> float:
        return self.absolute_pressure_pa / (Boltzmann * self.temperature_k)


@dataclass(frozen=True)
class RFDrive:
    frequency_hz: float
    empty_focus_field_peak_v_m: float

    def __post_init__(self) -> None:
        _positive("frequency_hz", self.frequency_hz)
        _finite("empty_focus_field_peak_v_m", self.empty_focus_field_peak_v_m)
        if self.empty_focus_field_peak_v_m < 0:
            raise ValueError("empty_focus_field_peak_v_m must be nonnegative")


@dataclass(frozen=True)
class Laser:
    wavelength_m: float
    energy_at_focus_j: float
    intensity_fwhm_s: float
    waist_radius_m: float
    beam_quality_m2: float
    phase_refractive_index: float
    group_refractive_index: float

    def __post_init__(self) -> None:
        for name, value in vars(self).items():
            _positive(name, value)
        if self.beam_quality_m2 < 1:
            raise ValueError("beam_quality_m2 must be at least one")


@dataclass(frozen=True)
class Copper:
    conductivity_s_m: float
    relative_permeability: float

    def __post_init__(self) -> None:
        _positive("conductivity_s_m", self.conductivity_s_m)
        _positive("relative_permeability", self.relative_permeability)


@dataclass(frozen=True)
class RectangularGeometry:
    width_m: float
    height_m: float
    closed_cavity_length_m: float

    def __post_init__(self) -> None:
        for name, value in vars(self).items():
            _positive(name, value)
        if self.width_m <= self.height_m:
            raise ValueError("rectangular width_m must exceed height_m")


@dataclass(frozen=True)
class CircularGeometry:
    diameter_m: float
    taper_length_m: float

    def __post_init__(self) -> None:
        _positive("diameter_m", self.diameter_m)
        _positive("taper_length_m", self.taper_length_m)


@dataclass(frozen=True)
class ReferenceCase:
    status: str
    description: str
    gas: Gas
    rf: RFDrive
    laser: Laser
    copper: Copper
    rectangular: RectangularGeometry
    circular: CircularGeometry

    def __post_init__(self) -> None:
        if self.status != "illustrative":
            raise ValueError("reference case must be explicitly illustrative")
        if not isinstance(self.description, str) or not self.description.strip():
            raise ValueError("description must not be empty")


def _object(raw: dict[str, Any], key: str) -> dict[str, Any]:
    value = raw.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"{key} must be a JSON object")
    return value


def load_reference_case(path: Path) -> ReferenceCase:
    with path.open(encoding="utf-8") as handle:
        raw = json.load(handle)
    if not isinstance(raw, dict):
        raise ValueError("reference case must be a JSON object")
    expected = {
        "status", "description", "gas", "rf", "laser",
        "copper", "rectangular", "circular",
    }
    if set(raw) != expected:
        raise ValueError(f"case keys must be {sorted(expected)}; got {sorted(raw)}")
    return ReferenceCase(
        status=raw["status"],
        description=raw["description"],
        gas=Gas(**_object(raw, "gas")),
        rf=RFDrive(**_object(raw, "rf")),
        laser=Laser(**_object(raw, "laser")),
        copper=Copper(**_object(raw, "copper")),
        rectangular=RectangularGeometry(**_object(raw, "rectangular")),
        circular=CircularGeometry(**_object(raw, "circular")),
    )
