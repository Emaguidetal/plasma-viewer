"""Verify saved partial-benchmark metrics; do not imply full multiphysics validation."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parents[2]


def _read(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _reciprocity_error(forward: Path, reverse: Path, frequency_hz: float) -> float:
    with np.load(forward / "modal-data.npz") as left, np.load(reverse / "modal-data.npz") as right:
        frequencies = left["frequency_hz"]
        if not np.array_equal(frequencies, right["frequency_hz"]):
            raise ValueError("Forward and reverse taper frequencies do not match.")
        index = int(np.argmin(abs(frequencies - frequency_hz)))
        if not math.isclose(float(frequencies[index]), frequency_hz, rel_tol=1e-12):
            raise ValueError("Taper spectrum does not contain its reported diagnostic frequency.")
        error = abs(left["RECT-TE10_raw_power_wave"][index] - right["TE11-90_raw_power_wave"][index])
    return float(error)


def evaluate(results: Path) -> tuple[dict[str, object], dict[str, dict]]:
    cold = results / "cold"
    names = {
        "guide3": "guide-3mm-cal", "guide2": "guide-2mm-cal",
        "closed3": "closed-cavity-3mm-r2", "closed2": "closed-cavity-2mm",
        "cavity3": "cavity-3mm", "cavity2": "cavity-2mm",
        "slab3": "drude-slab-3mm-prototype", "slab2": "drude-slab-2mm",
        "slab15": "drude-slab-1p5mm", "slab1": "drude-slab-1mm",
        "slab0": "drude-slab-zero-2mm",
        "circle0": "cartesian-circle-te11-0-3mm",
        "circle90": "cartesian-circle-te11-90-3mm",
        "circletm": "cartesian-circle-tm01-3mm",
        "taper3": "transition-te11-90-3mm-20slice",
        "taper2": "transition-te11-90-2mm-20slice",
        "taper2refined": "transition-te11-90-2mm-40slice",
    }
    data = {key: _read(cold / name / "metrics.json") for key, name in names.items()}
    data["transient300"] = _read(results / "transient" / "prescribed-300ps" / "metrics.json")
    data["transient3000"] = _read(results / "transient" / "prescribed-3000ps" / "metrics.json")
    diagnostic_names = {
        "2mm-20slice": ("transition-flux-forward-2mm", "transition-flux-reverse-2mm"),
        "1p5mm-40slice": ("transition-flux-forward-1p5mm-40slice", "transition-flux-reverse-1p5mm-40slice"),
        "1mm-40slice": ("transition-flux-forward-1mm-40slice", "transition-flux-reverse-1mm-40slice"),
    }
    for label, (forward, reverse) in diagnostic_names.items():
        data[f"{label}-forward"] = _read(cold / forward / "metrics.json")
        data[f"{label}-reverse"] = _read(cold / reverse / "metrics.json")
    official = list((results / "verification").glob("official-rect-*/metrics.json"))
    if len(official) != 1:
        raise ValueError(f"Expected exactly one official rectangular tutorial metric file; found {len(official)}")
    data["official"] = _read(official[0])

    checks = []

    def check(name: str, value: float, limit: float) -> None:
        if not math.isfinite(value):
            raise ValueError(f"{name} produced a non-finite value")
        checks.append({"name": name, "value": value, "upper_limit": limit, "passed": value < limit})

    check("official tutorial power-balance error", abs(data["official"]["modal_power_sum"] - 1), 0.01)
    check("official tutorial impedance relative error", data["official"]["impedance_relative_error"], 0.01)
    check("WR-284 guide reflection amplitude", data["guide2"]["center_abs_s11"], 0.01)
    check("WR-284 guide same-mode power-balance error", abs(data["guide2"]["center_modal_power_sum"] - 1), 0.01)
    check("WR-284 guide transmission phase error in radians", abs(data["guide2"]["guide_phase_error_rad"]), 0.01)
    check("isolated cavity frequency relative error", data["closed2"]["relative_frequency_error"], 0.01)
    check("isolated cavity frequency relative change from 3 to 2 mm",
          abs(data["closed3"]["simulated_te101_frequency_hz"]
              / data["closed2"]["simulated_te101_frequency_hz"] - 1), 0.01)
    check(
        "iris-cavity reflection change from 3 to 2 mm",
        abs(data["cavity3"]["center_abs_s11"] - data["cavity2"]["center_abs_s11"]), 0.02,
    )
    check(
        "iris-cavity peak-frequency relative change from 3 to 2 mm",
        abs(data["cavity3"]["cavity_resonance"]["frequency_hz"]
            / data["cavity2"]["cavity_resonance"]["frequency_hz"] - 1), 0.01,
    )
    check("zero-density slab reflection amplitude", data["slab0"]["center_abs_s11"], 0.01)
    check("zero-density slab transmission error", abs(data["slab0"]["center_abs_s21"] - 1), 0.01)
    check("1 mm slab reflection error against analytic solution", data["slab1"]["analytic_slab"]["abs_s11_error"], 0.02)
    check("1 mm slab transmission error against analytic solution", data["slab1"]["analytic_slab"]["abs_s21_error"], 0.02)
    check("slab reflection change from 1.5 to 1 mm",
          abs(data["slab15"]["center_abs_s11"] - data["slab1"]["center_abs_s11"]), 0.02)
    check("slab transmission change from 1.5 to 1 mm",
          abs(data["slab15"]["center_abs_s21"] - data["slab1"]["center_abs_s21"]), 0.02)
    for mode in ("circle0", "circle90", "circletm"):
        check(f"{mode} straight-guide same-mode power-balance error",
              abs(data[mode]["straight_same_mode_power_sum"] - 1), 0.01)
    check("taper input reflection change from 3 to 2 mm",
          abs(data["taper3"]["input_same_mode_abs_s11"] - data["taper2"]["input_same_mode_abs_s11"]), 0.02)
    check("taper input reflection change from 20 to 40 slices",
          abs(data["taper2"]["input_same_mode_abs_s11"] - data["taper2refined"]["input_same_mode_abs_s11"]), 0.02)
    for label in ("transient300", "transient3000"):
        check(f"{label} prescribed-resonator energy residual",
              data[label]["max_energy_residual_relative"], 0.02)
    taper_diagnostics = []

    def diagnostic(resolution: str, name: str, value: float, limit: float) -> None:
        if not math.isfinite(value):
            raise ValueError(f"{resolution} {name} produced a non-finite value")
        taper_diagnostics.append({
            "resolution": resolution, "name": name, "value": value,
            "upper_limit": limit, "passed": value < limit,
        })

    for resolution, (forward, reverse) in diagnostic_names.items():
        for direction in ("forward", "reverse"):
            result = data[f"{resolution}-{direction}"]
            diagnostic(resolution, f"{direction} raw modal power-sum error",
                       abs(result["raw_taper_modal_power_sum_not_independently_calibrated"] - 1), 0.01)
            diagnostic(resolution, f"{direction} independent field-flux imbalance",
                       result["field_flux_validation"]["relative_flux_imbalance"], 0.01)
            diagnostic(resolution, f"{direction} rectangular flux/modal-power discrepancy",
                       abs(result["field_flux_validation"]["field_to_modal_output_net_ratio"] - 1), 0.01)
        frequency = data[f"{resolution}-forward"]["frequency_hz"]
        if not math.isclose(frequency, data[f"{resolution}-reverse"]["frequency_hz"], rel_tol=1e-12):
            raise ValueError(f"{resolution} forward and reverse metric frequencies do not match.")
        diagnostic(resolution, "forward/reverse complex power-wave difference",
                   _reciprocity_error(cold / forward, cold / reverse, frequency), 0.02)
    summary = {
        "status": "partial numerical benchmarks only; the requested full simulation is not complete",
        "selected_baseline_checks_passed": all(item["passed"] for item in checks),
        "baseline_checks": checks,
        "finest_taper_diagnostic_gates_passed": all(
            item["passed"] for item in taper_diagnostics if item["resolution"] == "1mm-40slice"
        ),
        "taper_diagnostics": taper_diagnostics,
        "not_validated": [
            "Full 3D laser-plasma-RF coupling and spatial screening.",
            "Gas-specific He, N2, or Ar kinetics and 532 nm laser ionization.",
            "Copper losses in the 3D irised or tapered devices.",
            "Complete power-normalized circular multimode scattering and reciprocity.",
            "Optical window, measured cavity/taper geometry, and experimental comparison.",
            "The requested two-geometry by three-gas by two-duration physical case matrix.",
        ],
    }
    return summary, data


def save(results: Path, output: Path) -> dict[str, object]:
    results, output = results.expanduser().resolve(), output.expanduser().resolve()
    summary, data = evaluate(results)
    output.mkdir(parents=True, exist_ok=False)
    (output / "benchmark-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    mesh = np.array([3, 2, 1.5, 1])
    slab = [data[key] for key in ("slab3", "slab2", "slab15", "slab1")]
    axes[0, 0].plot(mesh, [item["analytic_slab"]["abs_s11_error"] for item in slab], "o-", label="|S11|")
    axes[0, 0].plot(mesh, [item["analytic_slab"]["abs_s21_error"] for item in slab], "o-", label="|S21|")
    axes[0, 0].axhline(0.02, color="tab:red", linestyle="--", label="0.02 tolerance")
    axes[0, 0].set(xlabel="Maximum mesh spacing (mm)", ylabel="Absolute analytic error", title="Drude slab")
    axes[0, 0].legend()
    axes[0, 1].plot([3, 2], [data[k]["simulated_te101_frequency_hz"] / 1e9 for k in ("closed3", "closed2")], "o-")
    axes[0, 1].axhline(data["closed2"]["analytic_te101_frequency_hz"] / 1e9, color="tab:red",
                       linestyle="--", label="Analytic")
    axes[0, 1].set(xlabel="Maximum mesh spacing (mm)", ylabel="Frequency (GHz)", title="Closed PEC cavity")
    axes[0, 1].legend()
    axes[1, 0].plot([3, 2], [data[k]["center_abs_s11"] for k in ("cavity3", "cavity2")], "o-")
    axes[1, 0].set(xlabel="Maximum mesh spacing (mm)", ylabel="|S11| at 3 GHz", title="Iris-coupled cavity")
    labels = ["3 mm / 20", "2 mm / 20", "2 mm / 40"]
    axes[1, 1].plot(labels, [data[k]["input_same_mode_abs_s11"]
                            for k in ("taper3", "taper2", "taper2refined")], "o-")
    axes[1, 1].set(xlabel="Mesh / taper slices", ylabel="Input |S11|", title="Exploratory taper")
    for axis in axes.flat:
        axis.grid(alpha=0.3)
    fig.suptitle("Partial numerical convergence; not full-system validation")
    fig.savefig(output / "numerical-convergence.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    labels = ["2 mm / 20", "1.5 mm / 40", "1 mm / 40"]
    resolutions = ["2mm-20slice", "1p5mm-40slice", "1mm-40slice"]
    quantities = (
        ("raw_taper_modal_power_sum_not_independently_calibrated", "Raw modal power-sum error"),
        ("relative_flux_imbalance", "Independent E/H flux imbalance"),
        ("field_to_modal_output_net_ratio", "Rectangular flux/modal discrepancy"),
    )
    for axis, (key, title) in zip(axes, quantities):
        for direction in ("forward", "reverse"):
            values = []
            for resolution in resolutions:
                record = data[f"{resolution}-{direction}"]
                if key in record:
                    value = abs(record[key] - 1)
                elif key == "field_to_modal_output_net_ratio":
                    value = abs(record["field_flux_validation"][key] - 1)
                else:
                    value = record["field_flux_validation"][key]
                values.append(100 * value)
            axis.plot(labels, values, "o-", label=direction)
        axis.axhline(1, color="tab:red", linestyle="--", label="1% target")
        axis.set(xlabel="Mesh / taper slices", ylabel="Discrepancy (%)", title=title)
        axis.tick_params(axis="x", labelrotation=15)
        axis.grid(alpha=0.3)
        axis.legend()
    fig.suptitle("Taper diagnostic gates; do not infer physical plasma loss")
    fig.savefig(output / "taper-power-diagnostics.png", dpi=160)
    plt.close(fig)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=PROJECT / "results")
    parser.add_argument("--output", type=Path, required=True, help="New output directory; never overwritten")
    args = parser.parse_args()
    summary = save(args.results, args.output)
    print(json.dumps(summary, indent=2))
    if not summary["selected_baseline_checks_passed"] or not summary["finest_taper_diagnostic_gates_passed"]:
        raise SystemExit("One or more recorded numerical validation gates remain unmet.")


if __name__ == "__main__":
    main()
