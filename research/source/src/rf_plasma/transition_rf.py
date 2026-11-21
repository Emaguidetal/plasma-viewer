"""Exploratory Cartesian 3D circular-guide and circular-to-WR-284 taper models."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import jn_zeros

from .circular_rf import tm01_impedance_ohm
from .cold_rf import GIB, _load_openems, available_memory_bytes, estimated_memory_bytes, memory_budget_bytes
from .config import load_reference_case
from .waveguides import circular_cutoff_hz

PROJECT = Path(__file__).resolve().parents[2]


def transition_polygon(
    fraction: float, radius_m: float, width_m: float, height_m: float, angular_samples: int = 64
) -> np.ndarray:
    if not math.isfinite(fraction) or not 0 <= fraction <= 1:
        raise ValueError("fraction must be between zero and one")
    if any(not math.isfinite(v) or v <= 0 for v in (radius_m, width_m, height_m)):
        raise ValueError("cross-section dimensions must be positive and finite")
    if isinstance(angular_samples, bool) or not isinstance(angular_samples, int) or angular_samples < 16:
        raise ValueError("at least 16 angular samples are required")
    corner = math.atan2(height_m, width_m)
    angles = np.unique(np.concatenate((
        np.linspace(0, 2 * math.pi, angular_samples, endpoint=False),
        [corner, math.pi - corner, math.pi + corner, 2 * math.pi - corner],
    )))
    cos, sin = np.cos(angles), np.sin(angles)
    x_limit = np.divide(width_m / 2, np.abs(cos), out=np.full_like(cos, np.inf), where=np.abs(cos) > 1e-15)
    y_limit = np.divide(height_m / 2, np.abs(sin), out=np.full_like(sin, np.inf), where=np.abs(sin) > 1e-15)
    rectangle_radius = np.minimum(x_limit, y_limit)
    smooth = fraction**2 * (3 - 2 * fraction)
    radii = (1 - smooth) * radius_m + smooth * rectangle_radius
    return np.vstack((radii * cos, radii * sin))


def _tm01_cartesian_profiles(radius_m: float) -> tuple[list[str], list[str], float]:
    kc = float(jn_zeros(0, 1)[0]) / radius_m
    rho = "sqrt(x*x+y*y+1e-30)"
    shape = f"j1({kc:.17g}*({rho}))/({rho})*(({rho})<{radius_m:.17g})"
    ex, ey = f"({shape})*x", f"({shape})*y"
    return [ex, ey, "0"], [f"-({ey})", ex, "0"], kc


def _power_wave_ratio(
    outgoing: np.ndarray, incident: np.ndarray, incident_impedance: np.ndarray, outgoing_impedance: np.ndarray
) -> np.ndarray:
    if (
        np.any(~np.isfinite(outgoing)) or np.any(~np.isfinite(incident))
        or np.any(~np.isfinite(incident_impedance)) or np.any(~np.isfinite(outgoing_impedance))
        or np.any(np.abs(incident) == 0) or np.any(np.real(incident_impedance) <= 0)
        or np.any(np.real(outgoing_impedance) <= 0)
        or np.any(np.imag(incident_impedance) != 0) or np.any(np.imag(outgoing_impedance) != 0)
    ):
        raise ValueError("Power-wave normalization requires finite waves and positive real impedances.")
    result = outgoing / incident * np.sqrt(np.real(incident_impedance) / np.real(outgoing_impedance))
    if not np.all(np.isfinite(result)):
        raise RuntimeError("The power-wave ratio is non-finite.")
    return result


def _field_flux(output: Path, label: str, frequency_hz: float) -> float:
    import h5py

    if not math.isfinite(frequency_hz) or frequency_hz <= 0:
        raise ValueError("The field-flux frequency must be positive and finite.")
    fields = []
    meshes = []
    frequencies = []
    for component in ("electric", "magnetic"):
        with h5py.File(output / f"{label}_{component}_flux.h5", "r") as handle:
            meshes.append(tuple(handle[f"Mesh/{axis}"][:] for axis in "xyz"))
            dataset = handle["FieldData/FD/f0"]
            frequencies.append(float(dataset.attrs["frequency"]))
            if dataset.ndim != 4 or dataset.shape[-1] != 1:
                raise RuntimeError(f"The {label} field dump is not a single plane.")
            fields.append(dataset[:, :, :, 0])
    if not all(np.array_equal(a, b) for a, b in zip(*meshes)):
        raise RuntimeError(f"The {label} electric and magnetic field meshes do not match.")
    if any(not math.isclose(value, frequency_hz, rel_tol=1e-12) for value in frequencies):
        raise RuntimeError(f"The {label} field dumps do not match the modal frequency.")
    electric, magnetic = fields
    x, y, z = meshes[0]
    if (
        len(z) != 1 or electric.shape != (3, len(x), len(y)) or magnetic.shape != electric.shape
        or len(x) < 2 or len(y) < 2 or np.any(np.diff(x) <= 0) or np.any(np.diff(y) <= 0)
        or any(not np.all(np.isfinite(axis)) for axis in (x, y, z))
        or not np.all(np.isfinite(electric)) or not np.all(np.isfinite(magnetic))
    ):
        raise RuntimeError(f"The {label} field-flux data are invalid.")
    flux_density = 0.5 * np.real(electric[0] * np.conj(magnetic[1]) - electric[1] * np.conj(magnetic[0]))
    flux = float(np.trapezoid(np.trapezoid(flux_density, y, axis=1), x))
    if not math.isfinite(flux):
        raise RuntimeError(f"The {label} integrated field flux is non-finite.")
    return flux


def run(
    kind: str, case_file: Path, output: Path, mesh_m: float, slices: int, threads: int,
    excited_mode: str = "TE11-90", flux_validation: bool = False,
) -> dict[str, object]:
    if kind not in {"straight", "taper"}:
        raise ValueError("kind must be straight or taper")
    if not math.isfinite(mesh_m) or not 0.001 <= mesh_m <= 0.003:
        raise ValueError("mesh spacing must be between 1 and 3 mm")
    if isinstance(slices, bool) or not isinstance(slices, int) or not 10 <= slices <= 100:
        raise ValueError("slices must be an integer between 10 and 100")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 8:
        raise ValueError("threads must be an integer between 1 and 8")
    if excited_mode not in {"TE11-0", "TE11-90", "TM01", "RECT-TE10"}:
        raise ValueError("excited mode must be TE11-0, TE11-90, TM01, or RECT-TE10")
    if excited_mode == "RECT-TE10" and kind != "taper":
        raise ValueError("RECT-TE10 excitation requires a taper.")
    case = load_reference_case(case_file)
    output = output.expanduser().resolve()
    radius = case.circular.diameter_m / 2
    a, b = case.rectangular.width_m, case.rectangular.height_m
    if max(a, b) > 2 * radius:
        raise ValueError("the rectangle must fit within this model's circular bounding box")
    f0 = case.rf.frequency_hz
    low, high = f0 - 0.35e9, f0 + 0.20e9
    if low <= circular_cutoff_hz(radius, "TM", 0, 1) or high >= circular_cutoff_hz(radius, "TE", 2, 1):
        raise ValueError("analysis band must lie between circular TM01 and TE21 cutoffs")
    ContinuousStructure, openEMS, _ = _load_openems()
    from openEMS.ports import Port

    inlet, taper_length, outlet = 0.09, case.circular.taper_length_m, 0.09
    length = inlet + taper_length + outlet if kind == "taper" else 0.25
    input_planes = (0.03, 0.045)
    output_planes = (length - 0.03, length - 0.045)
    solver = openEMS(NrTS=40000, EndCriteria=1e-6)
    solver.SetGaussExcite(f0, 0.35e9)
    solver.SetBoundaryCond([0, 0, 0, 0, 3, 3])
    csx = ContinuousStructure()
    solver.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1)
    grid.SetLines("x", [-radius, -a / 2, 0, a / 2, radius])
    grid.SetLines("y", [-radius, -b / 2, 0, b / 2, radius])
    grid.SetLines("z", [0, *input_planes, *output_planes, length])
    metal = csx.AddMetal("pec_background")
    metal.AddBox([-radius, -radius, 0], [radius, radius, length], priority=0)
    vacuum = csx.AddMaterial("vacuum_interior", epsilon=1)
    if kind == "straight":
        vacuum.AddCylinder([0, 0, 0], [0, 0, length], radius=radius, priority=10)
    else:
        vacuum.AddCylinder([0, 0, 0], [0, 0, inlet], radius=radius, priority=10)
        vacuum.AddBox([-a / 2, -b / 2, inlet + taper_length], [a / 2, b / 2, length], priority=10)
        edges = np.linspace(inlet, inlet + taper_length, slices + 1)
        grid.AddLine("z", edges)
        for index in range(slices):
            fraction = (index + 0.5) / slices
            points = transition_polygon(fraction, radius, a, b)
            vacuum.AddLinPoly(points, "z", edges[index], edges[index + 1] - edges[index], priority=10)
    grid.SmoothMeshLines("all", mesh_m, ratio=1.4)
    lines = [np.asarray(grid.GetLines(axis)) for axis in "xyz"]
    shape = [int(len(axis)) for axis in lines]
    points = math.prod(shape)
    memory = estimated_memory_bytes(points)
    budget = memory_budget_bytes(available_memory_bytes())
    if memory > budget:
        raise RuntimeError(f"Estimated memory {memory / GIB:.2f} GiB exceeds the {budget / GIB:.2f} GiB budget")

    def circle_box(z: tuple[float, float]) -> tuple[list[float], list[float]]:
        return [-radius, -radius, z[0]], [radius, radius, z[1]]

    modes = {}
    for index, (name, angle) in enumerate((("TE11-0", 0.0), ("TE11-90", math.pi / 2))):
        start, stop = circle_box(input_planes)
        input_port = solver.AddCircWaveGuidePort(
            2 * index, start, stop, "z", radius, "TE11", angle, 1 if name == excited_mode else 0
        )
        ports = [input_port]
        if kind == "straight":
            start, stop = circle_box(output_planes)
            ports.append(solver.AddCircWaveGuidePort(2 * index + 1, start, stop, "z", radius, "TE11", angle))
        modes[name] = ports
    electric, magnetic, kc = _tm01_cartesian_profiles(radius)
    start, stop = circle_box(input_planes)
    tm_input = solver.AddWaveGuidePort(
        4, start, stop, "z", E_func=electric, H_func=magnetic, kc=kc,
        excite=1 if excited_mode == "TM01" else 0, local_origin=[0, 0, 0]
    )
    modes["TM01"] = [tm_input]
    if kind == "straight":
        start, stop = circle_box(output_planes)
        modes["TM01"].append(solver.AddWaveGuidePort(
            5, start, stop, "z", E_func=electric, H_func=magnetic, kc=kc, local_origin=[0, 0, 0]
        ))
        rectangular = None
    else:
        rectangular = solver.AddRectWaveGuidePort(
            6, [-a / 2, -b / 2, output_planes[0]], [a / 2, b / 2, output_planes[1]],
            "z", a, b, "TE10", 1 if excited_mode == "RECT-TE10" else 0
        )

    if flux_validation:
        analysis_frequencies = np.linspace(low, high, 201)
        dump_frequency = float(analysis_frequencies[np.argmin(abs(analysis_frequencies - f0))])
        for label, z in (("input", input_planes[1]), ("output", output_planes[1])):
            for component, dump_type in (("electric", 10), ("magnetic", 11)):
                csx.AddDump(f"{label}_{component}_flux", dump_type=dump_type, dump_mode=1, file_type=1,
                            frequency=[dump_frequency]).AddBox([-radius, -radius, z], [radius, radius, z])

    output.mkdir(parents=True, exist_ok=False)
    print(f"Cartesian grid {shape}, estimated memory {memory / GIB:.3f} GiB; budget {budget / GIB:.3f} GiB")
    solver.Run(str(output), cleanup=False, numThreads=threads, dump_statistics=True)
    freq = np.linspace(low, high, 201)
    for name, ports in modes.items():
        for port in ports:
            if name == "TM01":
                Port.CalcPort(port, str(output), freq, ref_impedance=tm01_impedance_ohm(freq, radius))
            else:
                port.CalcPort(str(output), freq)
    if rectangular is not None:
        rectangular.CalcPort(str(output), freq)
    excited = rectangular if excited_mode == "RECT-TE10" else modes[excited_mode][0]
    inc = excited.uf_inc
    if np.any(~np.isfinite(inc)) or np.any(np.abs(inc) == 0):
        raise RuntimeError("The excited port has an invalid or zero incident amplitude.")
    reflected = excited.uf_ref / inc
    center = int(np.argmin(abs(freq - f0)))
    cross = {
        name: float(abs(ports[0].uf_ref[center] / inc[center]))
        for name, ports in modes.items() if name != excited_mode
    }
    raw_waves = {
        name: _power_wave_ratio(ports[0].uf_ref, inc, excited.Z_ref, ports[0].Z_ref)
        for name, ports in modes.items()
    }
    if rectangular is not None:
        raw_waves["RECT-TE10"] = _power_wave_ratio(rectangular.uf_ref, inc, excited.Z_ref, rectangular.Z_ref)
    if kind == "straight":
        transmitted = modes[excited_mode][1].uf_ref / inc
        output_measure = float(abs(transmitted[center]))
        power_sum = float(abs(reflected[center])**2 + abs(transmitted[center])**2)
        raw_output_power = None
    else:
        transmitted = None
        raw_output_power = rectangular.P_ref / excited.P_inc if excited_mode != "RECT-TE10" else None
        output_measure = float(raw_output_power[center]) if raw_output_power is not None else None
        power_sum = None
    if (
        not np.all(np.isfinite(reflected))
        or (output_measure is not None and not math.isfinite(output_measure))
        or (raw_output_power is not None and not np.all(np.isfinite(raw_output_power)))
    ):
        raise RuntimeError("The transition solver produced non-finite modal results.")
    result = {
        "status": "exploratory PEC geometry; circular modal power and reciprocity not yet validated",
        "kind": kind,
        "excited_mode": excited_mode,
        "frequency_hz": float(freq[center]),
        "mesh_max_m": mesh_m,
        "grid_lines": shape,
        "estimated_memory_bytes_not_hard_limit": memory,
        "initial_memory_budget_bytes": budget,
        "taper_length_m": taper_length if kind == "taper" else None,
        "taper_staircase_slices": slices if kind == "taper" else None,
        "input_same_mode_abs_s11": float(abs(reflected[center])),
        "straight_same_mode_abs_s21": output_measure if kind == "straight" else None,
        "straight_same_mode_power_sum": power_sum,
        "taper_rectangular_raw_power_ratio_not_calibrated": output_measure if kind == "taper" else None,
        "other_input_mode_voltage_ratios_not_power_normalized": (
            cross if excited_mode != "RECT-TE10" else None
        ),
        "other_circular_side_mode_voltage_ratios_not_power_normalized": cross,
        "raw_power_normalized_wave_amplitudes_not_independently_calibrated": {
            name: float(abs(wave[center])) for name, wave in raw_waves.items()
        },
        "raw_taper_modal_power_sum_not_independently_calibrated": (
            float(sum(abs(wave[center])**2 for wave in raw_waves.values())) if kind == "taper" else None
        ),
        "limitations": [
            "PEC walls and vacuum background; no window, copper loss, or plasma.",
            "Circular cross-mode powers and the custom TM01 port are not fully calibrated.",
            "Staircased taper geometry and Cartesian curved walls require independent mesh refinement.",
            "Do not interpret unaccounted power as absorption or use this as a validated full device result.",
        ],
    }
    if kind == "taper" and excited_mode == "RECT-TE10":
        result["taper_rectangular_raw_power_ratio_not_calibrated"] = None
    if flux_validation:
        input_flux = _field_flux(output, "input", dump_frequency)
        output_flux = _field_flux(output, "output", dump_frequency)
        expected_sign = -1 if excited_mode == "RECT-TE10" else 1
        if expected_sign * input_flux <= 0 or expected_sign * output_flux <= 0:
            raise RuntimeError("The independently integrated flux has an unexpected direction or zero magnitude.")
        circular_net = sum(port[0].P_acc[center] for port in modes.values())
        output_net = (
            -rectangular.P_acc[center] if rectangular is not None
            else -sum(port[1].P_acc[center] for port in modes.values())
        )
        if not math.isfinite(circular_net) or not math.isfinite(output_net) or circular_net == 0 or output_net == 0:
            raise RuntimeError("The integrated modal net power is zero or non-finite.")
        result["field_flux_validation"] = {
            "frequency_hz": dump_frequency,
            "net_input_flux": input_flux,
            "net_output_flux": output_flux,
            "relative_flux_imbalance": float(abs(input_flux - output_flux) / max(abs(input_flux), abs(output_flux))),
            "field_to_modal_input_net_ratio": float(input_flux / circular_net),
            "field_to_modal_output_net_ratio": float(output_flux / output_net),
            "note": (
                "Input/output refer to the low/high z planes; E/H flux and modal powers use "
                "the solver's Fourier normalization."
            ),
        }
    spectra = {
        "frequency_hz": freq,
        "input_reflection": reflected,
        **{f"input_{name}_voltage_ratio": ports[0].uf_ref / inc
           for name, ports in modes.items() if name != excited_mode},
        **{f"{name}_raw_power_wave": wave for name, wave in raw_waves.items()},
    }
    if transmitted is not None:
        spectra["same_mode_transmission"] = transmitted
    if raw_output_power is not None:
        spectra["rectangular_raw_power_ratio_not_calibrated"] = raw_output_power
    np.savez_compressed(output / "modal-data.npz", **spectra)
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axis = plt.subplots(figsize=(6, 5), constrained_layout=True)
    if kind == "taper":
        for fraction in (0, 0.25, 0.5, 0.75, 1):
            polygon = transition_polygon(fraction, radius, a, b)
            polygon = np.column_stack((polygon, polygon[:, 0]))
            axis.plot(polygon[0] * 1e3, polygon[1] * 1e3, label=f"{fraction:.2f}")
        axis.legend(title="Taper fraction")
    else:
        polygon = transition_polygon(0, radius, a, b)
        polygon = np.column_stack((polygon, polygon[:, 0]))
        axis.plot(polygon[0] * 1e3, polygon[1] * 1e3)
    axis.set(xlabel="x (mm)", ylabel="y (mm)", title="Illustrative waveguide cross-sections")
    axis.set_aspect("equal")
    axis.grid(alpha=0.3)
    fig.savefig(output / "cross-sections.png", dpi=160)
    plt.close(fig)
    fig, axis = plt.subplots(figsize=(7, 4), constrained_layout=True)
    axis.plot(freq / 1e9, abs(reflected), label="Input reflection amplitude")
    if transmitted is not None:
        axis.plot(freq / 1e9, abs(transmitted), label="Same-mode transmission amplitude")
    if raw_output_power is not None:
        axis.plot(freq / 1e9, raw_output_power, label="Raw rectangular power ratio (uncalibrated)")
    axis.set(xlabel="Frequency (GHz)", ylabel="Modal response",
             title="Exploratory PEC result; not a multimode power balance")
    axis.grid(alpha=0.3)
    axis.legend()
    fig.savefig(output / "modal-responses.png", dpi=160)
    plt.close(fig)
    (output / "metrics.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("straight", "taper"))
    parser.add_argument("--case", type=Path, default=PROJECT / "cases" / "reference.json")
    parser.add_argument("--output", type=Path, required=True, help="New output directory; never overwritten")
    parser.add_argument("--mesh-mm", type=float, default=3.0)
    parser.add_argument("--slices", type=int, default=20)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--mode", choices=("TE11-0", "TE11-90", "TM01", "RECT-TE10"), default="TE11-90")
    parser.add_argument("--flux-validation", action="store_true", help="Save and compare independent E/H flux planes")
    args = parser.parse_args()
    print(json.dumps(run(args.kind, args.case, args.output, args.mesh_mm * 1e-3,
                         args.slices, args.threads, args.mode, args.flux_validation), indent=2))


if __name__ == "__main__":
    main()
