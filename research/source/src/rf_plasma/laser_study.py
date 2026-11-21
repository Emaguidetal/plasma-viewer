"""Illustrative Gaussian-focus sensitivity; no ionization or propagation solver."""

import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from .config import ReferenceCase, load_reference_case
from .optics import beam_radius_m, gaussian_intensity_w_m2, ideal_focal_length_m, rayleigh_length_m

PROJECT = Path(__file__).resolve().parents[2]
WAISTS_M = (10e-6, 20e-6, 50e-6)
PULSES_S = (300e-12, 3e-9)


def calculate(case: ReferenceCase) -> list[dict[str, float]]:
    records = []
    for waist in WAISTS_M:
        for duration in PULSES_S:
            laser = replace(case.laser, waist_radius_m=waist, intensity_fwhm_s=duration)
            records.append({
                "waist_radius_m": waist,
                "intensity_fwhm_s": duration,
                "energy_at_focus_j": laser.energy_at_focus_j,
                "rayleigh_length_m": rayleigh_length_m(laser),
                "peak_intensity_w_m2": float(gaussian_intensity_w_m2(laser, 0, 0, 0)),
                "ideal_focal_length_for_1mm_input_radius_m": ideal_focal_length_m(laser, 1e-3),
            })
    return records


def save_result(case: ReferenceCase, output: Path) -> list[dict[str, float]]:
    output = output.expanduser().resolve()
    records = calculate(case)
    output.mkdir(parents=True, exist_ok=False)
    metadata = {
        "status": "illustrative Gaussian intensity only; not a breakdown prediction",
        "wavelength_m": case.laser.wavelength_m,
        "energy_at_focus_j": case.laser.energy_at_focus_j,
        "assumptions": [
            "Gaussian spatial and temporal intensity, specified energy reaches focus.",
            "No absorption, gas nonlinearity, plasma defocusing, or window properties.",
            "Beam quality and refractive indices are illustrative reference-case inputs.",
        ],
        "cases": records,
    }
    (output / "optics-sensitivity.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    time_ps = np.linspace(-3000, 3000, 1001)
    for duration in PULSES_S:
        laser = replace(case.laser, waist_radius_m=20e-6, intensity_fwhm_s=duration)
        intensity = gaussian_intensity_w_m2(laser, 0, 0, time_ps * 1e-12) / 1e4
        axes[0].plot(time_ps, intensity, label=f"{duration * 1e12:g} ps")
    axes[0].set(xlabel="Time at focus (ps)", ylabel=r"Intensity (W/cm$^2$)",
                title="20 µm waist; 1 mJ at focus")
    axes[0].legend()
    z_mm = np.linspace(-20, 20, 1001)
    for waist in WAISTS_M:
        laser = replace(case.laser, waist_radius_m=waist)
        axes[1].plot(z_mm, beam_radius_m(laser, z_mm * 1e-3) * 1e6, label=f"{waist * 1e6:g} µm")
    axes[1].set(xlabel="Distance from focus (mm)", ylabel="1/e² beam radius (µm)",
                title="Ideal Gaussian beam")
    axes[1].legend(title="Focus waist")
    for axis in axes:
        axis.grid(alpha=0.3)
    fig.suptitle("Optical sensitivity: illustrative, no gas or plasma propagation")
    fig.savefig(output / "optical-sensitivity.png", dpi=160)
    plt.close(fig)
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=PROJECT / "cases" / "reference.json")
    parser.add_argument("--output", type=Path, required=True, help="New output directory; never overwritten")
    args = parser.parse_args()
    case = load_reference_case(args.case)
    print(json.dumps(save_result(case, args.output), indent=2))


if __name__ == "__main__":
    main()
