"""Conservative single-ion reaction benchmark with externally supplied rates."""

from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import solve_ivp

from .config import _finite, _positive


@dataclass(frozen=True)
class ReactionRates:
    photoionization_s: float
    impact_ionization_m3_s: float
    recombination_m3_s: float

    def __post_init__(self) -> None:
        for name, value in vars(self).items():
            _finite(name, value)
            if value < 0:
                raise ValueError(f"{name} must be nonnegative")


def single_ion_benchmark(
    total_density_m3: float,
    initial_electron_density_m3: float,
    rates: ReactionRates,
    time_s: np.ndarray,
) -> dict[str, np.ndarray]:
    """Solve an abstract closed neutral/singly charged ion/electron system.

    Rates must be supplied externally. No gas, electron-energy, transport, or
    laser-deposition model is implied.
    """
    _positive("total_density_m3", total_density_m3)
    _finite("initial_electron_density_m3", initial_electron_density_m3)
    if not 0 <= initial_electron_density_m3 <= total_density_m3:
        raise ValueError("initial electron density must lie between zero and total density")
    time = np.asarray(time_s, dtype=float)
    if time.ndim != 1 or len(time) < 2 or not np.all(np.isfinite(time)) or np.any(np.diff(time) <= 0):
        raise ValueError("time_s must contain at least two finite, strictly increasing values")
    ionization = rates.impact_ionization_m3_s * total_density_m3
    recombination = rates.recombination_m3_s * total_density_m3
    if not math.isfinite(ionization) or not math.isfinite(recombination):
        raise ValueError("rate-density products must be finite")

    def derivative(_time: float, state: np.ndarray) -> list[float]:
        fraction = state[0]
        return [
            (rates.photoionization_s + ionization * fraction) * (1 - fraction)
            - recombination * fraction**2
        ]

    solution = solve_ivp(
        derivative, (float(time[0]), float(time[-1])), [initial_electron_density_m3 / total_density_m3],
        method="Radau", t_eval=time, rtol=1e-8, atol=1e-12,
    )
    if not solution.success or solution.y.shape != (1, len(time)) or not np.all(np.isfinite(solution.y)):
        raise RuntimeError(f"Single-ion benchmark failed: {solution.message}")
    fraction = solution.y[0]
    if np.min(fraction) < -1e-10 or np.max(fraction) > 1 + 1e-10:
        raise RuntimeError("Single-ion benchmark produced a nonphysical species fraction.")
    charged = total_density_m3 * fraction
    return {
        "time_s": time,
        "electron_density_m3": charged,
        "ion_density_m3": charged.copy(),
        "neutral_density_m3": total_density_m3 * (1 - fraction),
    }
