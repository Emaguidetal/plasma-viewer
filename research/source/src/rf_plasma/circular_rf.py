"""Verify TE11 and custom TM01 ports in a straight, PEC circular guide."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.constants import c, epsilon_0
from scipy.special import jn_zeros

from .cold_rf import GIB, _load_openems, available_memory_bytes, estimated_memory_bytes, memory_budget_bytes
from .config import load_reference_case
from .waveguides import circular_cutoff_hz

PROJECT = Path(__file__).resolve().parents[2]


def tm01_impedance_ohm(frequency_hz: np.ndarray, radius_m: float) -> np.ndarray:
    frequency = np.asarray(frequency_hz, dtype=float)
    if not math.isfinite(radius_m) or radius_m <= 0:
        raise ValueError("radius must be positive and finite")
    cutoff = circular_cutoff_hz(radius_m, "TM", 0, 1)
    if np.any(~np.isfinite(frequency)) or np.any(frequency <= cutoff):
        raise ValueError("all frequencies must exceed the TM01 cutoff")
    omega = 2 * math.pi * frequency
    beta = np.sqrt((omega / c) ** 2 - (float(jn_zeros(0, 1)[0]) / radius_m) ** 2)
    return beta / (omega * epsilon_0)


def run(
    case_file: Path, output: Path, mode: str, mesh_m: float, angular_cells: int, threads: int
) -> dict[str, object]:
    if mode not in {"TE11-0", "TE11-90", "TM01"}:
        raise ValueError("mode must be TE11-0, TE11-90, or TM01")
    if not math.isfinite(mesh_m) or not 0.0005 <= mesh_m <= 0.003:
        raise ValueError("mesh spacing must be between 0.5 and 3 mm")
    if isinstance(angular_cells, bool) or not isinstance(angular_cells, int) or not 24 <= angular_cells <= 160:
        raise ValueError("angular_cells must be an integer between 24 and 160")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 8:
        raise ValueError("threads must be an integer between 1 and 8")
    case = load_reference_case(case_file)
    output = output.expanduser().resolve()
    radius = case.circular.diameter_m / 2
    f0 = case.rf.frequency_hz
    low, high = f0 - 0.35e9, f0 + 0.20e9
    if low <= circular_cutoff_hz(radius, "TM", 0, 1):
        raise ValueError("analysis band must be above the TM01 cutoff")
    if high >= circular_cutoff_hz(radius, "TE", 2, 1):
        raise ValueError("analysis band must be below the TE21 cutoff")

    ContinuousStructure, openEMS, _ = _load_openems()
    from CSXCAD.CSRectGrid import CoordinateSystem
    from openEMS.ports import Port

    solver = openEMS(NrTS=80000, EndCriteria=1e-6, CoordSystem=1)
    solver.SetGaussExcite(f0, 0.35e9)
    solver.SetBoundaryCond([0, 0, 0, 0, 3, 3])
    csx = ContinuousStructure(CoordSystem=CoordinateSystem.CYLINDRICAL)
    solver.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1)
    length = 0.25
    port_in, port_out = (0.03, 0.045), (length - 0.03, length - 0.045)
    grid.SetLines("r", [0, radius])
    grid.SetLines("a", [0, 2 * math.pi])
    grid.SetLines("z", [0, *port_in, *port_out, length])
    for axis, resolution in zip("raz", (mesh_m, 2 * math.pi / angular_cells, mesh_m)):
        grid.SmoothMeshLines(axis, resolution, ratio=1.4)
    lines = [np.asarray(grid.GetLines(axis)) for axis in "raz"]
    shape = [int(len(axis)) for axis in lines]
    points = math.prod(shape)
    memory = estimated_memory_bytes(points)
    budget = memory_budget_bytes(available_memory_bytes())
    if memory > budget:
        raise RuntimeError(f"Estimated memory {memory / GIB:.2f} GiB exceeds the {budget / GIB:.2f} GiB budget")

    mode_ports = {}
    for index, (name, angle) in enumerate((("TE11-0", 0.0), ("TE11-90", math.pi / 2))):
        mode_ports[name] = [
            solver.AddCircWaveGuidePort(
                2 * index, [0, 0, port_in[0]], [radius, 2 * math.pi, port_in[1]],
                "z", radius, "TE11", angle, 1 if mode == name else 0
            ),
            solver.AddCircWaveGuidePort(
                2 * index + 1, [0, 0, port_out[0]], [radius, 2 * math.pi, port_out[1]],
                "z", radius, "TE11", angle
            ),
        ]
    kc = float(jn_zeros(0, 1)[0]) / radius
    radial = f"j1({kc:.17g}*rho)"
    e_profile = [radial, "0", "0"]
    h_profile = ["0", radial, "0"]
    mode_ports["TM01"] = [
        solver.AddWaveGuidePort(
            4, [0, 0, port_in[0]], [radius, 2 * math.pi, port_in[1]], "z",
            E_func=e_profile, H_func=h_profile, kc=kc, excite=1 if mode == "TM01" else 0
        ),
        solver.AddWaveGuidePort(
            5, [0, 0, port_out[0]], [radius, 2 * math.pi, port_out[1]], "z",
            E_func=e_profile, H_func=h_profile, kc=kc
        ),
    ]

    output.mkdir(parents=True, exist_ok=False)
    print(f"Cylindrical grid {shape}, estimated memory {memory / GIB:.3f} GiB; budget {budget / GIB:.3f} GiB")
    solver.Run(str(output), cleanup=False, numThreads=threads, dump_statistics=True)
    freq = np.linspace(low, high, 201)
    for name, ports in mode_ports.items():
        for port in ports:
            if name == "TM01":
                Port.CalcPort(port, str(output), freq, ref_impedance=tm01_impedance_ohm(freq, radius))
            else:
                port.CalcPort(str(output), freq)
    incident = mode_ports[mode][0].uf_inc
    if np.any(~np.isfinite(incident)) or np.any(np.abs(incident) == 0):
        raise RuntimeError("The excited port has an invalid or zero incident amplitude.")
    s11 = mode_ports[mode][0].uf_ref / incident
    s21 = mode_ports[mode][1].uf_ref / incident
    if not np.all(np.isfinite(s11)) or not np.all(np.isfinite(s21)):
        raise RuntimeError("The solver produced non-finite modal results.")
    center = int(np.argmin(np.abs(freq - f0)))
    cross = {
        name: {
            "input_abs_voltage_ratio": float(abs(ports[0].uf_ref[center] / incident[center])),
            "output_abs_voltage_ratio": float(abs(ports[1].uf_ref[center] / incident[center])),
        }
        for name, ports in mode_ports.items() if name != mode
    }
    result = {
        "status": "illustrative straight PEC circular-guide port benchmark; no plasma or taper",
        "excited_mode": mode,
        "frequency_hz": float(freq[center]),
        "radius_m": radius,
        "metal_model": "perfect electric conductor",
        "mesh_max_m": mesh_m,
        "angular_cells_requested": angular_cells,
        "grid_lines_r_azimuth_z": shape,
        "estimated_memory_bytes_not_hard_limit": memory,
        "initial_memory_budget_bytes": budget,
        "center_abs_s11_same_mode": float(abs(s11[center])),
        "center_abs_s21_same_mode": float(abs(s21[center])),
        "same_mode_power_sum": float(abs(s11[center]) ** 2 + abs(s21[center]) ** 2),
        "other_mode_voltage_ratios_not_power_normalized": cross,
        "limitations": [
            "Custom TM01 port requires straight-guide excitation and reciprocity verification before use in a taper.",
            "Cross-mode voltage ratios are not a complete, power-normalized multimode budget.",
            "No copper loss, optical window, plasma, taper, or experimental validation.",
        ],
    }
    np.savez_compressed(output / "modal-data.npz", frequency_hz=freq, s11=s11, s21=s21)
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axis = plt.subplots(figsize=(7, 4), constrained_layout=True)
    axis.plot(freq / 1e9, np.abs(s11), label="Reflected " + mode)
    axis.plot(freq / 1e9, np.abs(s21), label="Transmitted " + mode)
    axis.set(xlabel="Frequency (GHz)", ylabel="Same-mode amplitude", ylim=(0, 1.15))
    axis.grid(alpha=0.3)
    axis.legend()
    fig.savefig(output / "modal-amplitudes.png", dpi=160)
    plt.close(fig)
    (output / "metrics.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=PROJECT / "cases" / "reference.json")
    parser.add_argument("--output", type=Path, required=True, help="New output directory; never overwritten")
    parser.add_argument("--mode", choices=("TE11-0", "TE11-90", "TM01"), default="TE11-0")
    parser.add_argument("--mesh-mm", type=float, default=3.0)
    parser.add_argument("--angular-cells", type=int, default=50)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    result = run(args.case, args.output, args.mode, args.mesh_mm * 1e-3, args.angular_cells, args.threads)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
