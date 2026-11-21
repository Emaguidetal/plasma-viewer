"""Command-line entry point for explicitly labeled analytic benchmarks."""

import argparse
import json
import math
from pathlib import Path

from .config import load_reference_case
from .optics import gaussian_intensity_w_m2, ideal_focal_length_m, rayleigh_length_m
from .waveguides import (
    circular_cutoff_hz,
    isolated_te101_copper_q,
    rectangular_cutoff_hz,
    skin_depth_m,
    te10_copper_power_attenuation_per_m,
    te101_frequency_hz,
)


def describe_case(path: Path) -> dict[str, object]:
    case = load_reference_case(path)
    a = case.rectangular.width_m
    b = case.rectangular.height_m
    radius = case.circular.diameter_m / 2
    return {
        "status": "illustrative analytical checks; not a coupled plasma simulation",
        "gas": case.gas.species,
        "pressure_absolute_pa": case.gas.absolute_pressure_pa,
        "ideal_neutral_density_m3": case.gas.ideal_neutral_density_m3,
        "rf_frequency_hz": case.rf.frequency_hz,
        "empty_focus_field_peak_v_m": case.rf.empty_focus_field_peak_v_m,
        "peak_reduced_field_td": (
            case.rf.empty_focus_field_peak_v_m / case.gas.ideal_neutral_density_m3 / 1e-21
        ),
        "wr284_te10_cutoff_hz": rectangular_cutoff_hz(a, b, 1, 0),
        "wr284_te20_cutoff_hz": rectangular_cutoff_hz(a, b, 2, 0),
        "closed_cavity_te101_frequency_hz": te101_frequency_hz(case.rectangular),
        "circular_te11_cutoff_hz": circular_cutoff_hz(radius, "TE", 1, 1),
        "circular_tm01_cutoff_hz": circular_cutoff_hz(radius, "TM", 0, 1),
        "copper_skin_depth_m": skin_depth_m(case.rf.frequency_hz, case.copper),
        "te10_copper_power_attenuation_db_per_m_first_order": (
            10 / math.log(10)
            * te10_copper_power_attenuation_per_m(case.rf.frequency_hz, case.rectangular, case.copper)
        ),
        "ideal_closed_cavity_copper_q": isolated_te101_copper_q(case.rectangular, case.copper),
        "laser_rayleigh_length_m": rayleigh_length_m(case.laser),
        "laser_peak_intensity_w_m2": float(gaussian_intensity_w_m2(case.laser, 0, 0, 0)),
        "ideal_focal_length_for_1mm_input_radius_m": ideal_focal_length_m(case.laser, 1e-3),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path, help="Path to an illustrative JSON case")
    args = parser.parse_args()
    print(json.dumps(describe_case(args.case), indent=2))


if __name__ == "__main__":
    main()
