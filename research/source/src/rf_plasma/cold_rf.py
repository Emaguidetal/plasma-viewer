"""Run bounded, window-free PEC WR-284 and prescribed Drude benchmarks with openEMS."""

import argparse
import ctypes
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
from scipy.constants import c, e, electron_mass, epsilon_0
from scipy.signal import find_peaks

from .config import load_reference_case
from .drude import te10_slab_scattering
from .waveguides import (
    rectangular_cutoff_hz,
    te10_power_for_peak_electric_field_w,
    te10_propagation_constant_per_m,
)

GIB = 1024**3
PROJECT = Path(__file__).resolve().parents[2]


def memory_budget_bytes(available_bytes: int) -> int:
    if isinstance(available_bytes, bool) or not isinstance(available_bytes, int) or available_bytes <= 0:
        raise ValueError("available memory must be a positive integer number of bytes")
    return min(8 * GIB, available_bytes // 4)


def available_memory_bytes() -> int:
    if sys.platform == "win32":
        from ctypes import wintypes

        class MemoryStatus(ctypes.Structure):
            _fields_ = [
                ("length", wintypes.DWORD),
                ("load", wintypes.DWORD),
                ("total_physical", ctypes.c_ulonglong),
                ("available_physical", ctypes.c_ulonglong),
                ("total_page_file", ctypes.c_ulonglong),
                ("available_page_file", ctypes.c_ulonglong),
                ("total_virtual", ctypes.c_ulonglong),
                ("available_virtual", ctypes.c_ulonglong),
                ("available_extended_virtual", ctypes.c_ulonglong),
            ]

        state = MemoryStatus()
        state.length = ctypes.sizeof(state)
        if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(state)):
            raise OSError("GlobalMemoryStatusEx failed")
        return int(state.available_physical)
    return int(os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE"))


def estimated_memory_bytes(grid_points: int) -> int:
    if isinstance(grid_points, bool) or not isinstance(grid_points, int) or grid_points <= 0:
        raise ValueError("grid_points must be a positive integer")
    return 256 * 1024**2 + 512 * grid_points


def _load_openems():
    if sys.platform == "win32":
        tool_dir = PROJECT / ".tools" / "openEMS"
        if not tool_dir.is_dir() and "CSXCAD_INSTALL_PATH" not in os.environ:
            raise RuntimeError(f"openEMS DLL directory not found: {tool_dir}")
        os.environ.setdefault("CSXCAD_INSTALL_PATH", str(tool_dir))
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS
    from openEMS.utilities import DFT_time2freq

    return ContinuousStructure, openEMS, DFT_time2freq


def _plot_results(output: Path, freq: np.ndarray, s11: np.ndarray, s21: np.ndarray) -> None:
    import h5py
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axis = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
    axis.plot(freq / 1e9, np.abs(s11), label="Reflected TE10")
    axis.plot(freq / 1e9, np.abs(s21), label="Transmitted TE10")
    axis.set(xlabel="Frequency (GHz)", ylabel="Modal voltage amplitude", ylim=(0, 1.15))
    axis.grid(alpha=0.3)
    axis.legend()
    fig.savefig(output / "modal-amplitudes.png", dpi=160)
    plt.close(fig)

    with h5py.File(output / "electric_midplane.h5", "r") as field_file:
        x = field_file["Mesh/x"][:] * 1e3
        z = field_file["Mesh/z"][:] * 1e3
        ey = np.abs(field_file["FieldData/FD/f0"][1, :, 0, :])
    if not np.all(np.isfinite(ey)) or np.max(ey) <= 0:
        raise RuntimeError("The frequency-domain field map is invalid.")
    fig, axis = plt.subplots(figsize=(8, 3.6), constrained_layout=True)
    image = axis.pcolormesh(z, x, ey / np.max(ey), shading="auto", cmap="inferno")
    axis.set(xlabel="Axial position z (mm)", ylabel="Width position x (mm)",
             title="Relative |Ey| on the midplane at the reference frequency")
    fig.colorbar(image, ax=axis, label="Normalized field magnitude")
    fig.savefig(output / "electric-midplane.png", dpi=160)
    plt.close(fig)


def _resonance(freq: np.ndarray, amplitude: np.ndarray) -> dict[str, float | None] | None:
    peaks, _ = find_peaks(amplitude)
    if len(peaks) == 0:
        return None
    peak = int(peaks[np.argmax(amplitude[peaks])])
    half_power = amplitude[peak] / math.sqrt(2)
    left = np.flatnonzero(amplitude[:peak] <= half_power)
    right = np.flatnonzero(amplitude[peak + 1 :] <= half_power) + peak + 1
    result: dict[str, float | None] = {
        "frequency_hz": float(freq[peak]),
        "loaded_q_half_power_approximate": None,
    }
    if len(left) == 0 or len(right) == 0:
        return result
    left_idx = int(left[-1])
    right_idx = int(right[0])
    lower = np.interp(half_power, amplitude[left_idx : left_idx + 2], freq[left_idx : left_idx + 2])
    upper = np.interp(
        half_power, amplitude[right_idx - 1 : right_idx + 1][::-1], freq[right_idx - 1 : right_idx + 1][::-1]
    )
    result["loaded_q_half_power_approximate"] = float(freq[peak] / (upper - lower))
    return result


def run(
    kind: str, case_file: Path, output: Path, mesh_m: float, threads: int,
    calibration_file: Path | None = None, density_m3: float = 1e17, collision_rate_s: float = 2e10,
) -> dict[str, object]:
    if kind not in {"guide", "cavity", "slab"}:
        raise ValueError("kind must be guide, cavity, or slab")
    if not math.isfinite(mesh_m) or not 0.0005 <= mesh_m <= 0.003:
        raise ValueError("mesh spacing must be between 0.5 and 3 mm")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 8:
        raise ValueError("threads must be an integer between 1 and 8")
    case = load_reference_case(case_file)
    output = output.expanduser().resolve()
    if case.rf.empty_focus_field_peak_v_m <= 0:
        raise ValueError("the cold-device normalization requires a positive peak target field")
    if kind != "cavity" and calibration_file is not None:
        raise ValueError("guide calibration is only used for the cavity benchmark")
    if kind == "slab" and (
        not math.isfinite(density_m3) or density_m3 < 0
        or not math.isfinite(collision_rate_s) or collision_rate_s <= 0
    ):
        raise ValueError("slab density must be nonnegative and collision rate positive and finite")
    if calibration_file is not None:
        calibration_file = calibration_file.expanduser().resolve()
    center_frequency = case.rf.frequency_hz
    a, b = case.rectangular.width_m, case.rectangular.height_m
    low, high = center_frequency - 0.45e9, center_frequency + 0.45e9
    next_cutoff = min(rectangular_cutoff_hz(a, b, 2, 0), rectangular_cutoff_hz(a, b, 0, 1))
    if low <= rectangular_cutoff_hz(a, b, 1, 0) or high >= next_cutoff:
        raise ValueError("the entire analysis band must be in the rectangular TE10-only passband")
    cavity_length = case.rectangular.closed_cavity_length_m
    iris_z = (0.09, 0.09 + cavity_length)
    length = iris_z[1] + 0.09
    port_in = (0.03, 0.045)
    port_out = (length - 0.03, length - 0.045)
    focus_z = 0.5 * (iris_z[0] + iris_z[1])
    iris_opening = 0.045
    iris_thickness = 0.002
    slab_thickness = 0.012
    if iris_opening >= a:
        raise ValueError("the illustrative iris opening must be narrower than the guide")

    ContinuousStructure, openEMS, dft = _load_openems()
    solver = openEMS(NrTS=40000, EndCriteria=1e-6)
    solver.SetGaussExcite(center_frequency, 0.45e9)
    solver.SetBoundaryCond([0, 0, 0, 0, 3, 3])
    csx = ContinuousStructure()
    solver.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1)
    left, right = (a - iris_opening) / 2, (a + iris_opening) / 2
    grid.AddLine("x", [0, left, a / 2, right, a])
    grid.AddLine("y", [0, b / 2, b])
    grid.AddLine("z", [0, *port_in, *port_out, focus_z, length])

    if kind == "cavity":
        iris = csx.AddMetal("illustrative_pec_irises")
        for z in iris_z:
            grid.AddLine("z", [z - iris_thickness / 2, z + iris_thickness / 2])
            iris.AddBox([0, 0, z - iris_thickness / 2], [left, b, z + iris_thickness / 2])
            iris.AddBox([right, 0, z - iris_thickness / 2], [a, b, z + iris_thickness / 2])
    if kind == "slab":
        slab_start, slab_stop = focus_z - slab_thickness / 2, focus_z + slab_thickness / 2
        grid.AddLine("z", [slab_start, slab_stop])
        if density_m3 > 0:
            from CSXCAD.CSProperties import CSPropLorentzMaterial

            material = CSPropLorentzMaterial(csx.GetParameterSet(), order=1, epsilon=1)
            material.SetName("prescribed_drude_plasma")
            csx.AddProperty(material)
            plasma_frequency = math.sqrt(density_m3 * e**2 / (epsilon_0 * electron_mass)) / (2 * math.pi)
            material.SetDispersiveMaterialProperty(
                0, eps_plasma=plasma_frequency, eps_pole_freq=0, eps_relax=1 / collision_rate_s
            )
            material.AddBox([0, 0, slab_start], [a, b, slab_stop], priority=10)

    ports = [
        solver.AddRectWaveGuidePort(
            0, [0, 0, port_in[0]], [a, b, port_in[1]], "z", a, b, "TE10", 1
        ),
        solver.AddRectWaveGuidePort(
            1, [0, 0, port_out[0]], [a, b, port_out[1]], "z", a, b, "TE10"
        ),
    ]
    focus = [a / 2, b / 2, focus_z]
    csx.AddProbe("focus_electric_field", 2).AddBox(focus, focus)
    csx.AddDump("electric_midplane", dump_type=10, dump_mode=2, file_type=1,
                frequency=[center_frequency]).AddBox([0, b / 2, 0], [a, b / 2, length])
    grid.SmoothMeshLines("all", mesh_m, ratio=1.4)
    lines = [np.asarray(grid.GetLines(axis)) for axis in "xyz"]
    grid_shape = [int(len(axis)) for axis in lines]
    grid_points = math.prod(grid_shape)
    memory_estimate = estimated_memory_bytes(grid_points)
    available = available_memory_bytes()
    budget = memory_budget_bytes(available)
    if memory_estimate > budget:
        raise RuntimeError(
            f"Estimated memory {memory_estimate / GIB:.2f} GiB exceeds the "
            f"current {budget / GIB:.2f} GiB budget; use a coarser mesh."
        )
    calibration = None
    if calibration_file is not None:
        with calibration_file.open(encoding="utf-8") as handle:
            calibration = json.load(handle)
        if (
            calibration["kind"] != "guide"
            or calibration["frequency_hz"] != center_frequency
            or calibration["target_focus_peak_field_v_m"] != case.rf.empty_focus_field_peak_v_m
            or calibration["center_abs_s11"] >= 0.01
            or abs(calibration["center_modal_power_sum"] - 1) >= 0.01
            or not np.array_equal(calibration["mesh_lines_x_m"], lines[0])
            or not np.array_equal(calibration["mesh_lines_y_m"], lines[1])
        ):
            raise ValueError("guide calibration must match frequency, target, transverse mesh and validation checks")

    output.mkdir(parents=True, exist_ok=False)
    print(f"Grid {grid_shape}, estimated memory {memory_estimate / GIB:.3f} GiB; budget {budget / GIB:.3f} GiB")
    solver.Run(str(output), cleanup=False, numThreads=threads, dump_statistics=True)
    freq = np.linspace(center_frequency - 0.45e9, center_frequency + 0.45e9, 301)
    for port in ports:
        port.CalcPort(str(output), freq)
    s11 = ports[0].uf_ref / ports[0].uf_inc
    s21 = ports[1].uf_ref / ports[0].uf_inc
    if not (np.all(np.isfinite(s11)) and np.all(np.isfinite(s21))):
        raise RuntimeError("The solver produced non-finite modal results.")
    midpoint = len(freq) // 2
    field_data = np.loadtxt(output / "focus_electric_field", comments="%")
    if field_data.ndim != 2 or field_data.shape[1] < 4:
        raise RuntimeError("Unexpected focus E-field probe output.")
    field_spectrum = dft(field_data[:, 0], field_data[:, 2], [center_frequency])[0]
    if abs(field_spectrum) == 0 or ports[0].P_inc[midpoint] <= 0:
        raise RuntimeError("Focus field or incident modal power is zero.")
    target_power = float(
        ports[0].P_inc[midpoint] * (case.rf.empty_focus_field_peak_v_m / abs(field_spectrum)) ** 2
    )
    power_sum = float(abs(s11[midpoint]) ** 2 + abs(s21[midpoint]) ** 2)
    beta = te10_propagation_constant_per_m(center_frequency, a).real
    plane_separation = port_out[1] - port_in[1]
    phase_error = float(np.angle(s21[midpoint] * np.exp(1j * beta * plane_separation)))
    analytic_power = te10_power_for_peak_electric_field_w(
        center_frequency, a, b, case.rf.empty_focus_field_peak_v_m
    )
    calibrated_power = None
    if calibration is not None and kind != "slab":
        calibrated_power = target_power * (
            calibration["straight_guide_analytic_target_power_w"]
            / calibration["target_incident_power_w_unvalidated_until_guide_normalization_check"]
        )
    analytic_slab = None
    if kind == "slab":
        r_analytic, t_analytic = te10_slab_scattering(
            center_frequency, a, slab_thickness, density_m3, collision_rate_s
        )
        analytic_slab = {
            "abs_s11": abs(r_analytic),
            "abs_s21": abs(t_analytic),
            "absorbed_power_fraction": 1 - abs(r_analytic)**2 - abs(t_analytic)**2,
            "abs_s11_error": abs(abs(s11[midpoint]) - abs(r_analytic)),
            "abs_s21_error": abs(abs(s21[midpoint]) - abs(t_analytic)),
        }
    result = {
        "status": (
            "illustrative 3D prescribed-Drude slab benchmark; not plasma dynamics"
            if kind == "slab" else "illustrative 3D cold-device PEC benchmark; no plasma, copper loss, or window"
        ),
        "kind": kind,
        "frequency_hz": center_frequency,
        "metal_model": "perfect electric conductor",
        "mesh_max_m": mesh_m,
        "grid_lines": grid_shape,
        "mesh_lines_x_m": lines[0].tolist(),
        "mesh_lines_y_m": lines[1].tolist(),
        "grid_points": grid_points,
        "minimum_cell_m": [float(np.min(np.diff(axis))) for axis in lines],
        "estimated_memory_bytes_not_hard_limit": memory_estimate,
        "initial_memory_budget_bytes": budget,
        "focus_position_m": focus,
        "port_measurement_planes_z_m": [port_in[1], port_out[1]],
        "target_focus_peak_field_v_m": case.rf.empty_focus_field_peak_v_m,
        "target_incident_power_w_unvalidated_until_guide_normalization_check": target_power if kind != "slab" else None,
        "calibrated_target_incident_power_w": calibrated_power,
        "guide_calibration_file": str(calibration_file) if calibration_file is not None else None,
        "straight_guide_analytic_target_power_w": analytic_power,
        "center_abs_s11": float(abs(s11[midpoint])),
        "center_abs_s21": float(abs(s21[midpoint])),
        "center_modal_power_sum": power_sum,
        "guide_phase_error_rad": phase_error if kind == "guide" else None,
        "cavity_resonance": _resonance(freq, np.abs(s21)) if kind == "cavity" else None,
        "iris_opening_m": iris_opening if kind == "cavity" else None,
        "iris_thickness_m": iris_thickness if kind == "cavity" else None,
        "slab_thickness_m": slab_thickness if kind == "slab" else None,
        "slab_density_m3_prescribed": density_m3 if kind == "slab" else None,
        "slab_collision_rate_s_prescribed": collision_rate_s if kind == "slab" else None,
        "analytic_slab": analytic_slab,
        "limitations": [
            "Perfect-conductor walls, vacuum background and no optical window.",
            "Iris dimensions, prescribed slab parameters and the 3 GHz frequency are illustrative.",
            "A memory estimate is not a hard process-memory limit.",
            "Field normalization is approximate and mesh convergence requires independent verification.",
            "A prescribed Drude slab does not predict laser ionization or evolving plasma.",
        ],
    }
    np.savez_compressed(output / "modal-data.npz", frequency_hz=freq, s11=s11, s21=s21)
    _plot_results(output, freq, s11, s21)
    (output / "metrics.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("guide", "cavity", "slab"))
    parser.add_argument("--case", type=Path, default=PROJECT / "cases" / "reference.json")
    parser.add_argument("--output", type=Path, required=True, help="New output directory; never overwritten")
    parser.add_argument("--mesh-mm", type=float, default=3.0)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--calibration", type=Path, help="Matching verified guide metrics for field-power calibration")
    parser.add_argument("--density-m3", type=float, default=1e17, help="Prescribed slab electron density")
    parser.add_argument("--collision-rate-s", type=float, default=2e10, help="Prescribed slab collision rate")
    args = parser.parse_args()
    result = run(args.kind, args.case, args.output, args.mesh_mm * 1e-3, args.threads,
                 args.calibration, args.density_m3, args.collision_rate_s)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
