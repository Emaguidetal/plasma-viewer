"""Cycle-resolved passive resonator/current verification with a prescribed plasma."""

import argparse
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np
from scipy.constants import e, electron_mass, epsilon_0
from scipy.integrate import solve_ivp
from scipy.special import erf

from .config import ReferenceCase, _positive, load_reference_case
from .optics import rayleigh_length_m
from .waveguides import isolated_te101_copper_q, te101_frequency_hz

PROJECT = Path(__file__).resolve().parents[2]


def density_and_growth(
    time_s: float | np.ndarray, seed_m3: float, peak_m3: float, pulse_fwhm_s: float
) -> tuple[np.ndarray, np.ndarray]:
    _positive("seed_m3", seed_m3)
    _positive("peak_m3", peak_m3)
    _positive("pulse_fwhm_s", pulse_fwhm_s)
    if peak_m3 < seed_m3:
        raise ValueError("peak density must not be smaller than the seed density")
    time = np.asarray(time_s, dtype=float)
    factor = 2 * math.sqrt(math.log(2)) / pulse_fwhm_s
    phase = factor * time
    density = seed_m3 + (peak_m3 - seed_m3) * 0.5 * (1 + erf(phase))
    growth = (peak_m3 - seed_m3) * factor / math.sqrt(math.pi) * np.exp(-(phase**2))
    return density, growth


def simulate(
    case: ReferenceCase,
    pulse_fwhm_s: float,
    peak_density_m3: float,
    seed_density_m3: float,
    collision_rate_s: float,
    max_step_s: float,
) -> dict[str, object]:
    for name, value in (
        ("pulse_fwhm_s", pulse_fwhm_s),
        ("peak_density_m3", peak_density_m3),
        ("seed_density_m3", seed_density_m3),
        ("collision_rate_s", collision_rate_s),
        ("max_step_s", max_step_s),
    ):
        _positive(name, value)
    if peak_density_m3 < seed_density_m3:
        raise ValueError("peak density must not be smaller than the seed density")
    if peak_density_m3 > case.gas.ideal_neutral_density_m3:
        raise ValueError("prescribed electron density must not exceed the illustrative neutral density")
    if max_step_s > min(8e-12, pulse_fwhm_s / 40):
        raise ValueError("max step must be at most 8 ps and one-fortieth of the pulse FWHM")

    laser = replace(case.laser, intensity_fwhm_s=pulse_fwhm_s)
    geom = case.rectangular
    field_volume = geom.width_m * geom.height_m * geom.closed_cavity_length_m / 4
    plasma_volume = math.pi * laser.waist_radius_m**2 * (2 * rayleigh_length_m(laser))
    scale = math.sqrt(epsilon_0 * field_volume / 2)
    if plasma_volume >= field_volume:
        raise ValueError("the assumed compact plasma volume must be smaller than the resonator mode volume")
    omega0 = 2 * math.pi * te101_frequency_hz(geom)
    omega = 2 * math.pi * case.rf.frequency_hz
    q_external_1 = q_external_2 = 400.0
    q_copper = isolated_te101_copper_q(geom, case.copper)
    gamma1, gamma2, gamma_wall = omega0 / q_external_1, omega0 / q_external_2, omega0 / q_copper
    gamma = gamma1 + gamma2 + gamma_wall
    coupling = plasma_volume / (2 * scale)
    k_seed = seed_density_m3 * e**2 / electron_mass
    current_scale = peak_density_m3 * e**2 / electron_mass * case.rf.empty_focus_field_peak_v_m / collision_rate_s
    _positive("current scale", current_scale)
    t0, t1 = -3 * pulse_fwhm_s, 3 * pulse_fwhm_s + 2e-9
    samples = int(math.ceil((t1 - t0) / max_step_s)) + 1
    if samples > 100_000:
        raise ValueError("requested time interval and sampling exceed the bounded demonstration size")

    loading = plasma_volume * k_seed / (epsilon_0 * field_volume * (collision_rate_s - 1j * omega))
    response_denominator = gamma + 1j * (omega0**2 / omega - omega) + loading
    input_peak = (
        scale * case.rf.empty_focus_field_peak_v_m * abs(response_denominator) / math.sqrt(2 * gamma1)
    )
    a_e_phasor = math.sqrt(2 * gamma1) * input_peak / response_denominator
    a_h_phasor = -1j * omega0 / omega * a_e_phasor
    current_phasor = k_seed * (a_e_phasor / scale) / (collision_rate_s - 1j * omega)
    phase0 = np.exp(-1j * omega * t0)
    state0 = np.zeros(9)
    state0[0] = (a_e_phasor * phase0).real
    state0[1] = (a_h_phasor * phase0).real
    state0[2] = (current_phasor * phase0).real / current_scale

    def rhs(time: float, state: np.ndarray) -> np.ndarray:
        a_e, a_h, scaled_current = state[:3]
        current = scaled_current * current_scale
        density, growth = density_and_growth(time, seed_density_m3, peak_density_m3, pulse_fwhm_s)
        k = float(density) * e**2 / electron_mass
        e_field = a_e / scale
        incoming = input_peak * math.cos(omega * time)
        reflected = incoming - math.sqrt(2 * gamma1) * a_e
        transmitted = -math.sqrt(2 * gamma2) * a_e
        drift_energy = plasma_volume * current**2 / (2 * k)
        return np.array([
            omega0 * a_h - gamma * a_e + math.sqrt(2 * gamma1) * incoming - coupling * current,
            -omega0 * a_e,
            (k * e_field - collision_rate_s * current) / current_scale,
            incoming**2,
            reflected**2,
            transmitted**2,
            2 * gamma_wall * a_e**2,
            plasma_volume * collision_rate_s * current**2 / k,
            drift_energy * float(growth / density),
        ])

    time = np.linspace(t0, t1, samples)
    solution = solve_ivp(rhs, (t0, t1), state0, method="Radau", t_eval=time,
                         rtol=2e-8, atol=1e-10, max_step=max_step_s)
    if not solution.success or solution.y.shape[1] != len(time) or not np.all(np.isfinite(solution.y)):
        raise RuntimeError(f"Transient integration failed: {solution.message}")
    a_e, a_h, scaled_current = solution.y[:3]
    current = scaled_current * current_scale
    density, growth = density_and_growth(time, seed_density_m3, peak_density_m3, pulse_fwhm_s)
    k = density * e**2 / electron_mass
    field = a_e / scale
    cavity_energy = a_e**2 + a_h**2
    drift_energy = plasma_volume * current**2 / (2 * k)
    collisional_power = plasma_volume * collision_rate_s * current**2 / k
    mixing_power = drift_energy * growth / density
    input_energy, reflected_energy, transmitted_energy, wall_energy, collisional_energy, mixing_energy = solution.y[3:]
    residual = (
        cavity_energy[0] + drift_energy[0] + input_energy
        - reflected_energy - transmitted_energy - wall_energy - collisional_energy - mixing_energy
        - cavity_energy - drift_energy
    )
    energy_scale = cavity_energy[0] + drift_energy[0] + input_energy[-1]
    if energy_scale <= 0:
        raise RuntimeError("The energy-balance normalization is zero.")
    metrics = {
        "status": "prescribed-density 0D resonator/current verification; not a gas breakdown or 3D prediction",
        "pulse_fwhm_s": pulse_fwhm_s,
        "seed_density_m3_prescribed": seed_density_m3,
        "peak_density_m3_prescribed": peak_density_m3,
        "collision_rate_s_prescribed": collision_rate_s,
        "plasma_volume_m3_illustrative": plasma_volume,
        "electric_mode_volume_m3_ideal_closed_cavity": field_volume,
        "rf_frequency_hz": case.rf.frequency_hz,
        "oscillator_frequency_hz_ideal_closed_cavity": omega0 / (2 * math.pi),
        "q_external_1_assumed": q_external_1,
        "q_external_2_assumed": q_external_2,
        "q_copper_ideal_closed_cavity": q_copper,
        "initial_rf_peak_field_v_m": case.rf.empty_focus_field_peak_v_m,
        "input_average_power_w_illustrative": input_peak**2 / 2,
        "max_step_s": max_step_s,
        "integrator": "SciPy Radau",
        "final_incident_energy_j": float(input_energy[-1]),
        "final_reflected_energy_j": float(reflected_energy[-1]),
        "final_transmitted_energy_j": float(transmitted_energy[-1]),
        "final_wall_energy_j": float(wall_energy[-1]),
        "final_collisional_energy_j": float(collisional_energy[-1]),
        "final_birth_mixing_energy_j": float(mixing_energy[-1]),
        "max_energy_residual_relative": float(np.max(np.abs(residual)) / energy_scale),
        "limitations": [
            "Density is imposed as the integral of a Gaussian; it is not calculated from laser ionization or chemistry.",
            "Collision rate, seed density, external Qs and compact uniform plasma volume are illustrative assumptions.",
            "The lumped oscillator is not a validated reduction of either 3D device.",
            "No spatial screening, charge transport, gas-specific reactions, optical depletion or experimental validation.",
            "Birth-mixing energy records drift-energy dilution from adding electrons at rest.",
        ],
    }
    return {
        "metrics": metrics,
        "time_s": time,
        "density_m3": density,
        "field_v_m": field,
        "current_a_m2": current,
        "cavity_energy_j": cavity_energy,
        "drift_energy_j": drift_energy,
        "collisional_power_w": collisional_power,
        "birth_mixing_power_w": mixing_power,
        "energy_residual_j": residual,
    }


def save_result(output: Path, result: dict[str, object]) -> None:
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=False)
    metrics = result["metrics"]
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    np.savez_compressed(output / "transient.npz", **{k: v for k, v in result.items() if k != "metrics"})

    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    time_ps = np.asarray(result["time_s"]) * 1e12
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True, constrained_layout=True)
    axes[0].plot(time_ps, np.asarray(result["density_m3"]) / 1e22)
    axes[0].set_ylabel(r"Prescribed $n_e$ ($10^{22}$ m$^{-3}$)")
    axes[1].plot(time_ps, np.asarray(result["field_v_m"]) / 1e6)
    axes[1].set_ylabel("Field (MV/m)")
    axes[2].plot(time_ps, np.asarray(result["collisional_power_w"]) / 1e3, label="Collisional")
    axes[2].plot(time_ps, np.asarray(result["birth_mixing_power_w"]) / 1e3, label="Birth mixing")
    axes[2].set(xlabel="Time relative to prescribed pulse peak (ps)", ylabel="Power (kW)")
    axes[2].legend()
    for axis in axes:
        axis.grid(alpha=0.3)
    fig.suptitle("Prescribed-plasma verification, not a breakdown prediction")
    fig.savefig(output / "transient-traces.png", dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=PROJECT / "cases" / "reference.json")
    parser.add_argument("--output", type=Path, required=True, help="New output directory; never overwritten")
    parser.add_argument("--pulse-ps", type=float, default=300)
    parser.add_argument("--peak-density-m3", type=float, default=1e22)
    parser.add_argument("--seed-density-m3", type=float, default=1e15)
    parser.add_argument("--collision-rate-s", type=float, default=2e11)
    parser.add_argument("--max-step-ps", type=float, default=7.5)
    args = parser.parse_args()
    result = simulate(load_reference_case(args.case), args.pulse_ps * 1e-12, args.peak_density_m3,
                      args.seed_density_m3, args.collision_rate_s, args.max_step_ps * 1e-12)
    save_result(args.output, result)
    print(json.dumps(result["metrics"], indent=2))


if __name__ == "__main__":
    main()
