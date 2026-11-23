"""Evidence-bounded screening of Ar metastables in the existing gas comparison.

This is deliberately not a metastable-kinetics solver.  The workspace contains
state-specific ground-state excitation cross sections and secondary-source
Maxwellian rate fits for stepwise ionization, but no identified Ar* pooling,
quenching, or metastable diffusion closure at the pressures in question.  The
module extracts the existing, unmodified He/Ar/N2 baseline and calculates
transparent rate-scale diagnostics without feeding unsupported rates back into
the avalanche solver.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import gas_compare as gc
from .streamer2d import SWARM_FILES
from .stochastic_seed import Rates
from .swarm import read_table

ROOT = Path(__file__).resolve().parents[2]
BASELINE = ROOT / "results" / "gas_compare" / "gas_compare.json"
OUT = ROOT / "results" / "gas_compare" / "argon_metastable"

# NIST ASD Ar I levels (also used as thresholds in the local IST-Lisbon file).
AR_LEVELS_EV = {
    "Ar(3P2)": 11.548,
    "Ar(3P1)": 11.623,
    "Ar(3P0)": 11.723,
    "Ar(1P1)": 11.828,
}
AR_IP_EV = 15.7596119
KCONV_M6_S_300K = 2.25e-43
KEMANECI = {
    "k_exc_m1": (5.0e-15, 12.64),
    "k_exc_m2": (1.4e-15, 12.42),
    "k_step": (6.8e-15, 4.20),
}


def _fit_rates(te_eV: float) -> dict[str, float]:
    """Evaluate the Kemaneci 2017 Table 3 fits (Te in eV, coefficients in SI)."""
    if te_eV <= 0:
        raise ValueError("Te_eV must be positive")
    k_exc = sum(a * math.exp(-b / te_eV) for a, b in
                (KEMANECI["k_exc_m1"], KEMANECI["k_exc_m2"]))
    a, b = KEMANECI["k_step"]
    k_step = a * te_eV**0.67 * math.exp(-b / te_eV)
    return {"k_exc_two_4s_fits_m3_s": k_exc, "k_step_m3_s": k_step}


def _threshold_diagnostics(baseline: dict) -> dict:
    ar = Rates("Ar", path=SWARM_FILES["Ar"])
    swarm = read_table(SWARM_FILES["Ar"])
    en_grid = np.asarray(swarm["E_N_Td"], float)
    mean_e = np.asarray(swarm["mean_energy_eV"], float)
    rows = []
    for row in baseline["avalanche"]["gases"]["Ar"]:
        psig = float(row["psig"])
        st = gc.setup(psig)
        e0 = float(row["thr_env_kV_cm"]) * 1e5
        en_eff, nu_m = gc.eff_en_td(ar, e0, st.n_gas, st.f_hz)
        eps_bar = float(np.interp(math.log(en_eff), np.log(en_grid), mean_e))
        # Maxwellian-equivalent Te = (2/3)<epsilon>; this is only an explicit
        # mapping for evaluating source fits, not a Maxwellian EEDF assertion.
        te = 2.0 * eps_bar / 3.0
        fit = _fit_rates(te)
        kiz = float(ar.kiz(en_eff))
        nu_ground = st.n_gas * kiz
        # This fraction would make the fitted stepwise frequency equal the
        # baseline ground-state ionization frequency at the same field.
        f_equal = kiz / fit["k_step_m3_s"]
        kconv = KCONV_M6_S_300K * (st.t_k / 300.0) ** -0.4
        nu_conv = kconv * st.n_gas**2
        rows.append({
            "psig": psig,
            "p_abs_Pa": st.p_pa,
            "n_Ar_m3_300K": st.n_gas,
            "baseline_Ar_threshold_kV_cm_peak": float(row["thr_env_kV_cm"]),
            "threshold_Eeff_N_Td": en_eff,
            "nu_m_s_at_threshold": nu_m,
            "BOLOS_mean_energy_eV": eps_bar,
            "Te_equivalent_eV_for_Kemaneci_fits_assumption": te,
            **fit,
            "k_ground_ionization_Boltzmann_m3_s": kiz,
            "nu_ground_ionization_s": nu_ground,
            "metastable_fraction_for_equal_stepwise_and_ground_frequency": f_equal,
            "k_Arplus_to_Ar2plus_m6_s": kconv,
            "Arplus_to_Ar2plus_conversion_frequency_s": nu_conv,
            "Arplus_to_Ar2plus_conversion_time_ps": 1e12 / nu_conv,
            "stepwise_Te_bracket_3_to_5_eV": {
                str(t): _fit_rates(t)["k_step_m3_s"] for t in (3.0, 5.0)
            },
            "equal_frequency_metastable_fraction_Te_3_to_5_eV": {
                str(t): kiz / _fit_rates(t)["k_step_m3_s"] for t in (3.0, 5.0)
            },
            "qualification": (
                "Rate-scale screening only: source fits are Maxwellian/low-pressure; "
                "Te mapping is not a validated EEDF closure at 2.95 GHz or 50-200 psig."
            ),
        })
    return {
        "state_energy_eV": AR_LEVELS_EV,
        "Ar_ground_ionization_energy_eV": AR_IP_EV,
        "stepwise_energy_increment_from_Ar_3P2_eV": AR_IP_EV - AR_LEVELS_EV["Ar(3P2)"],
        "rows": rows,
    }


def build_assessment(baseline_path: Path = BASELINE) -> dict:
    """Create the screening report from the preceding same-pressure comparison."""
    baseline = json.loads(Path(baseline_path).read_text(encoding="utf-8"))
    avalanche = baseline["avalanche"]
    pressure_rows = []
    for i, psig in enumerate(gc.PSIG):
        ar = avalanche["gases"]["Ar"][i]
        he = avalanche["gases"]["He"][i]
        n2 = avalanche["gases"]["N2"][i]
        pairwise = {}
        for field_index, field_kv in enumerate(gc.FIELDS_KV_CM):
            by_gas = {
                gas: avalanche["gases"][gas][i]["growth"][field_index]["mc"]
                for gas in ("He", "Ar", "N2")
            }
            if (by_gas["Ar"]["delay_median_ns"] is not None
                    and by_gas["He"]["delay_median_ns"] is not None):
                pairwise[str(field_kv)] = {
                    "Ar_to_He_median_delay_ratio": (
                        by_gas["Ar"]["delay_median_ns"] / by_gas["He"]["delay_median_ns"]
                    ),
                    "Ar_to_He_jitter_p90_p10_ratio": (
                        by_gas["Ar"]["jitter_p90_p10_ns"] / by_gas["He"]["jitter_p90_p10_ns"]
                        if by_gas["He"]["jitter_p90_p10_ns"] else None
                    ),
                }
        pressure_rows.append({
            "psig": psig,
            "thresholds_kV_cm_peak": {
                gas: float(avalanche["gases"][gas][i]["thr_env_kV_cm"])
                for gas in ("He", "Ar", "N2")
            },
            "Ar_to_He_threshold_ratio": float(ar["thr_env_kV_cm"] / he["thr_env_kV_cm"]),
            "common_peak_field_delay_distributions": {
                gas: [{
                    "E0_kV_cm_peak": float(g["E0_kV_cm"]),
                    "breakdown_probability": g["mc"]["p"],
                    "median_delay_ns": g["mc"]["delay_median_ns"],
                    "p10_delay_ns": g["mc"]["delay_p10_ns"],
                    "p90_delay_ns": g["mc"]["delay_p90_ns"],
                    "jitter_p90_minus_p10_ns": g["mc"]["jitter_p90_p10_ns"],
                } for g in avalanche["gases"][gas][i]["growth"]]
                for gas in ("He", "Ar", "N2")
            },
            "same_field_Ar_to_He_ratios_when_both_break_down": pairwise,
            "near_threshold_delay_distributions": {
                gas: [{
                    "overvoltage": float(o["overvoltage"]),
                    "E0_kV_cm_peak": float(o["E0_kV_cm"]),
                    **o["mc"],
                } for o in avalanche["gases"][gas][i]["overvoltage"]]
                for gas in ("He", "Ar", "N2")
            },
        })

    return {
        "title": "Argon metastable screening for the same-pressure He/N2/Ar RF switch comparison",
        "status": "bounded_evidence_review_no_metastable_feedback_coupled",
        "conditions": {
            "pressures_psig": list(gc.PSIG),
            "temperature_K": 300.0,
            "frequency_Hz": 2.95e9,
            "baseline_peak_field_comparison_kV_cm": list(gc.FIELDS_KV_CM),
            "baseline_seed_assumption": {
                "model": "same He fitted k=6 seed law for all gases",
                "N0_electrons": avalanche["n0_common"],
                "required_log_gain": avalanche["g_req"],
            },
            "e_cal_used": 1.0,
            "He_fitted_e_cal_applied_to_other_gases": False,
            "baseline_source": "results/gas_compare/gas_compare.json (preceding stage; 300 Monte Carlo shots, RNG seed 1)",
        },
        "available_data": {
            "argon_excitation_cross_sections": {
                "file": "cases/cross_sections/LoKI_Ar_LXCat.txt",
                "source": "IST-Lisbon LXCat LoKI Ar set; excitation thresholds and terms from Khakoo et al. 2004",
                "resolved_ground_to_levels_eV": AR_LEVELS_EV,
                "source_condition": "Local LXCat file states 'All rights reserved'; this report does not redistribute cross-section arrays.",
            },
            "excitation_and_stepwise_rate_fits": {
                "source": "Kemaneci et al. 2017, J. Phys. D 50, 245203, Table 3, reactions 2, 3 and 8; secondary fits (reactions 2-3 based on Tachibana 1986).",
                "citation": "doi:10.1088/1361-6463/aa7093; arXiv:1612.07268",
                "rates": {
                    "Ar + e -> Ar(4s metastable) + e": "5.0e-15 exp(-12.64/Te_eV) m^3/s",
                    "Ar + e -> second Ar(4s metastable) + e": "1.4e-15 exp(-12.42/Te_eV) m^3/s",
                    "Ar(4s_m) + e -> Ar+ + 2e": "6.8e-15 Te_eV^0.67 exp(-4.20/Te_eV) m^3/s",
                    "applicability": "Maxwellian-style fits in a low-pressure microwave global model; validity ranges and state-term mapping are not stated in the local compilation.",
                },
            },
            "ionic_three_body_conversion": {
                "reaction": "Ar+ + 2Ar -> Ar2+ + Ar",
                "coefficient": "2.25e-43 (Tg/300 K)^-0.4 m^6/s",
                "source": "Kemaneci et al. 2017 Table 3 reaction 30; secondary of Kabouzi et al., Phys. Rev. E 75, 016402 (2007).",
                "caveat": "This is ion conversion, not metastable quenching or a source of electrons.",
            },
            "unavailable_critical_metastable_closures": [
                "Ar(4s) + Ar(4s) pooling / associative-ionization coefficient with identified product branches",
                "Ar(4s) + Ar quenching and three-body conversion coefficients at 300 K and 4.46-14.80 bar absolute",
                "Ar(4s)-Ar diffusion coefficient or state-specific transport cross section at these densities",
                "validated electron-impact stepwise-ionization cross section for each resolved Ar(3P2/3P0) level or a documented EEDF-resolved rate table",
            ],
        },
        "baseline_comparison": pressure_rows,
        "metastable_rate_scale_screen": _threshold_diagnostics(baseline),
        "interpretation": {
            "can_explain_lower_Ar_threshold": (
                "Plausible in principle: the Ar(3P2) state is at 11.548 eV versus Ar I ionization at "
                f"{AR_IP_EV:.6f} eV, leaving only {AR_IP_EV - AR_LEVELS_EV['Ar(3P2)']:.4f} eV to the ion threshold. "
                "The cited stepwise fit can therefore be significant if the metastable fraction is sufficiently large."
            ),
            "can_explain_longer_delay_or_jitter": (
                "The unmodified baseline reproduces the longer-delay/more-jitter direction at 200 psig and "
                "200 kV/cm (Ar median 39.67 ns and p90-p10 5.29 ns; He 29.55 ns and 3.32 ns), but it does "
                "not reproduce the lower Ar threshold. At equal relative overvoltage, Ar is faster and has "
                "similar or slightly lower jitter. Metastable lifetime and production fluctuations are "
                "unconstrained, so their incremental effects cannot be calculated."
            ),
            "no_claim_of_success": True,
            "high_pressure_extrapolation": (
                "All Kemaneci rate fits derive from a low-pressure microwave global model with unspecified "
                "validity ranges; 50-200 psig is 4.46-14.80 bar absolute, so any use as predictive kinetics "
                "would be an explicit, presently unvalidated extrapolation."
            ),
            "conservation_and_scope": (
                "No time-dependent metastable population or electron source was added. For the stepwise reaction "
                "Ar* + e -> Ar+ + 2e, one incoming electron becomes two (net +1), one Ar nucleus remains, and "
                "the ionization-potential rise is paid from electron kinetic energy plus the consumed Ar* "
                "excitation reservoir. A coupled implementation must track both reservoirs and all competing "
                "metastable losses; the current effective-field swarm closure has no state-resolved energy ledger."
            ),
            "narrowest_next_step": (
                "Measure time-resolved absolute Ar(3P2) and Ar(3P0) densities (absorption/TALIF calibrated for "
                "pressure broadening) during the first 0-100 ns at each pressure, alongside Te/EEDF; pair this "
                "with state-resolved Ar*-Ar quenching/pooling and diffusion data or a validated high-pressure "
                "collisional-radiative dataset. Then add a conservative metastable/electron/ion reservoir ODE."
            ),
        },
        "sources": [
            "IST-Lisbon LXCat / LoKI Ar local file (retrieved 2022-09-23; cross-section data rights reserved; see file header).",
            "Khakoo et al., J. Phys. B: At. Mol. Opt. Phys. 37 (2004) 247, as cited in LoKI_Ar_LXCat.txt.",
            "Kemaneci et al., J. Phys. D: Appl. Phys. 50 (2017) 245203, doi:10.1088/1361-6463/aa7093; local gas_kinetics_literature.json records accessed-source caveats.",
            "NIST Atomic Spectra Database: Ar I level energies and ionization energy (also transcribed in local gas-model provenance).",
        ],
    }


def figures(report: dict, out: Path = OUT) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out.mkdir(parents=True, exist_ok=True)
    ps = [x["psig"] for x in report["baseline_comparison"]]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
    for gas, color in (("He", "tab:blue"), ("Ar", "tab:orange"), ("N2", "tab:red")):
        vals = [x["thresholds_kV_cm_peak"][gas] for x in report["baseline_comparison"]]
        ax[0].plot(ps, vals, "o-", color=color, label=gas)
    ax[0].set(xlabel="pressure (psig)", ylabel="baseline threshold (kV/cm peak)",
              title="Existing model; no metastable correction")
    ax[0].legend()
    ax[0].grid(alpha=0.3)

    p200 = report["baseline_comparison"][-1]
    for gas, color in (("He", "tab:blue"), ("Ar", "tab:orange"), ("N2", "tab:red")):
        vals = p200["common_peak_field_delay_distributions"][gas]
        ax[1].plot([x["E0_kV_cm_peak"] for x in vals],
                   [np.nan if x["median_delay_ns"] is None else x["median_delay_ns"] for x in vals],
                   "o-", color=color, label=gas)
    ax[1].set(xlabel="common peak field (kV/cm)", ylabel="median formation delay (ns)",
              title="200 psig, same peak fields")
    ax[1].legend()
    ax[1].grid(alpha=0.3)

    rates = report["metastable_rate_scale_screen"]["rows"]
    ax[2].semilogy(ps, [r["metastable_fraction_for_equal_stepwise_and_ground_frequency"] for r in rates],
                   "o-", color="tab:purple", label="Te = 2/3 mean energy")
    ax[2].fill_between(
        ps,
        [min(r["equal_frequency_metastable_fraction_Te_3_to_5_eV"].values()) for r in rates],
        [max(r["equal_frequency_metastable_fraction_Te_3_to_5_eV"].values()) for r in rates],
        color="tab:purple", alpha=0.18, label="illustrative Te = 3–5 eV bracket",
    )
    ax[2].set(xlabel="pressure (psig)", ylabel="Ar* fraction for equal fitted stepwise / ground ν",
              title="Rate-scale diagnostic, not a prediction")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3, which="both")
    fig.suptitle("Argon metastables: evidence and closure gap (baseline uses common 376-electron seed)")
    fig.tight_layout()
    path = out / "argon_metastable_screen.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main(out: Path = OUT) -> dict:
    report = build_assessment()
    out.mkdir(parents=True, exist_ok=True)
    (out / "argon_metastable_assessment.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    figures(report, out)
    return report


if __name__ == "__main__":
    r = main()
    print(r["status"])
    for x in r["metastable_rate_scale_screen"]["rows"]:
        print(
            f"{x['psig']:3.0f} psig: Ar threshold {x['baseline_Ar_threshold_kV_cm_peak']:.2f} kV/cm, "
            f"f_equal={x['metastable_fraction_for_equal_stepwise_and_ground_frequency']:.3g}, "
            f"Ar+ conversion time={x['Arplus_to_Ar2plus_conversion_time_ps']:.2f} ps"
        )
