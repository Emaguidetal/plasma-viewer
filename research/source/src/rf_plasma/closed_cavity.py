"""Verify the PEC WR-284 TE101 cavity eigenfrequency with openEMS."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from .cold_rf import GIB, _load_openems, available_memory_bytes, estimated_memory_bytes, memory_budget_bytes
from .config import load_reference_case
from .waveguides import isolated_te101_copper_q, te101_frequency_hz

PROJECT = Path(__file__).resolve().parents[2]


def run(case_file: Path, output: Path, mesh_m: float, threads: int) -> dict[str, object]:
    if not math.isfinite(mesh_m) or not 0.0005 <= mesh_m <= 0.003:
        raise ValueError("mesh spacing must be between 0.5 and 3 mm")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 8:
        raise ValueError("threads must be an integer between 1 and 8")
    case = load_reference_case(case_file)
    output = output.expanduser().resolve()
    a, b, d = case.rectangular.width_m, case.rectangular.height_m, case.rectangular.closed_cavity_length_m
    expected = te101_frequency_hz(case.rectangular)
    ContinuousStructure, openEMS, _ = _load_openems()
    solver = openEMS(NrTS=20000, EndCriteria=0)
    solver.SetGaussExcite(expected, 0.5e9)
    solver.SetBoundaryCond([0, 0, 0, 0, 0, 0])
    csx = ContinuousStructure()
    solver.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1)
    grid.SetLines("x", [0, a / 4, a / 2, a])
    grid.SetLines("y", [0, b / 2, b])
    grid.SetLines("z", [0, d / 4, d / 2, d])
    grid.SmoothMeshLines("all", mesh_m, ratio=1.4)
    lines = [np.asarray(grid.GetLines(axis)) for axis in "xyz"]
    shape = [int(len(axis)) for axis in lines]
    points = math.prod(shape)
    memory = estimated_memory_bytes(points)
    budget = memory_budget_bytes(available_memory_bytes())
    if memory > budget:
        raise RuntimeError(f"Estimated memory {memory / GIB:.2f} GiB exceeds the {budget / GIB:.2f} GiB budget")
    source_start, source_stop = [a / 4, 0, d / 4], [a / 4, b, d / 4]
    focus = [a / 2, b / 2, d / 2]
    csx.AddExcitation("electric_source", 0, [0, 1, 0]).AddBox(source_start, source_stop)
    csx.AddProbe("focus_electric_field", 2).AddBox(focus, focus)

    output.mkdir(parents=True, exist_ok=False)
    print(f"Closed-cavity grid {shape}, estimated memory {memory / GIB:.3f} GiB; budget {budget / GIB:.3f} GiB")
    solver.Run(str(output), cleanup=False, numThreads=threads, dump_statistics=True)
    data = np.loadtxt(output / "focus_electric_field", comments="%")
    if data.ndim != 2 or data.shape[1] < 4 or len(data) < 100:
        raise RuntimeError("Unexpected closed-cavity E-field probe output.")
    time, ey = data[:, 0], data[:, 2]
    keep = time > 8e-9
    if np.count_nonzero(keep) < 100:
        raise RuntimeError("The simulation did not record enough post-excitation cycles.")
    signal = ey[keep] * np.hanning(np.count_nonzero(keep))
    sample_dt = float(np.mean(np.diff(time)))
    nfft = 2 ** int(math.ceil(math.log2(len(signal) * 16)))
    freq = np.fft.rfftfreq(nfft, sample_dt)
    amplitude = np.abs(np.fft.rfft(signal, nfft))
    band = (freq > expected - 0.3e9) & (freq < expected + 0.3e9)
    if not np.any(band) or np.max(amplitude[band]) <= 0:
        raise RuntimeError("The expected eigenmode is absent from the computed spectrum.")
    measured = float(freq[band][np.argmax(amplitude[band])])
    result = {
        "status": "illustrative 3D isolated PEC cavity eigenfrequency benchmark; no plasma or coupling apertures",
        "mesh_max_m": mesh_m,
        "grid_lines": shape,
        "estimated_memory_bytes_not_hard_limit": memory,
        "initial_memory_budget_bytes": budget,
        "analytic_te101_frequency_hz": expected,
        "simulated_te101_frequency_hz": measured,
        "relative_frequency_error": abs(measured / expected - 1),
        "ideal_closed_cavity_copper_q_analytic_not_simulated": isolated_te101_copper_q(case.rectangular, case.copper),
        "limitations": [
            "The 3D field solve uses perfect conductors; the copper Q is a separate first-order analytic estimate.",
            "No irises, window, plasma, gas permittivity correction, or experimental validation.",
        ],
    }
    np.savez_compressed(output / "spectrum.npz", frequency_hz=freq[band], amplitude=amplitude[band])
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axis = plt.subplots(figsize=(7, 4), constrained_layout=True)
    axis.plot(freq[band] / 1e9, amplitude[band] / np.max(amplitude[band]), label="3D PEC spectrum")
    axis.axvline(expected / 1e9, color="tab:red", linestyle="--", label="Analytic TE101")
    axis.set(xlabel="Frequency (GHz)", ylabel="Normalized spectrum")
    axis.grid(alpha=0.3)
    axis.legend()
    fig.savefig(output / "eigenfrequency.png", dpi=160)
    plt.close(fig)
    (output / "metrics.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=PROJECT / "cases" / "reference.json")
    parser.add_argument("--output", type=Path, required=True, help="New output directory; never overwritten")
    parser.add_argument("--mesh-mm", type=float, default=3.0)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    print(json.dumps(run(args.case, args.output, args.mesh_mm * 1e-3, args.threads), indent=2))


if __name__ == "__main__":
    main()
