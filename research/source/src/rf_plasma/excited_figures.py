"""Figures for the excited-state / dense-regime breakdown study (results/excited_avalanche).

python -m rf_plasma.excited_figures
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "excited_avalanche"
COL = {"He": "tab:blue", "Ar": "tab:green", "N2": "tab:red"}
PSIG = (50.0, 100.0, 150.0, 200.0, 300.0)


def _thr():
    d = json.loads((RES / "thresholds.json").read_text())
    out = {}
    for x in d:
        if x.get("thr_kV_cm") is not None:
            out[(x["mode"], x["variant"], x["gas"], float(x["psig"]))] = x["thr_kV_cm"]
    return out


def fig_mechanisms(plt):
    t = _thr()
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for gas in ("He", "Ar", "N2"):
        for var, ls, lab in (("direct_norec", ":", "direct only"), ("nominal", "-", "+ excited states (nominal)"),
                             ("hm1", "--", "+ Hornbeck-Molnar b = 1")):
            y = [t.get(("laser", var, gas, p)) for p in PSIG]
            ax[0].plot(PSIG, [np.nan if v is None else v for v in y], ls, c=COL[gas], marker="o", ms=3,
                       label=f"{gas}: {lab}")
    ax[0].set(yscale="log", xlabel="pressure (psig)", ylabel="threshold peak field (kV/cm)",
              title="Laser-seeded (376 electrons), dilute avalanche")
    ax[0].legend(fontsize=6.5, ncol=1)
    ax[0].grid(alpha=0.3, which="both")
    variants = [("direct_norec", "direct"), ("nominal", "nominal"), ("nominal_low", "low rates"),
                ("nominal_high", "high rates"), ("hm1", "HM b=1"), ("eta1e-1", "photoemission η=0.1"),
                ("eta1e-1_hm1_high", "all maximal"), ("nominal_screen", "screening criterion")]
    for k, mode in enumerate(("laser", "rf")):
        a = ax[1 + k]
        vlist = variants + ([("eta3e-1", "photoemission η=0.3"), ("eta1", "photoemission η=1 (bound)")]
                            if mode == "rf" else [])
        for i, (v, lab) in enumerate(vlist):
            r = [(t.get((mode, v, "Ar", p)), t.get((mode, v, "He", p))) for p in PSIG]
            y = [a_ / b_ if (a_ and b_) else np.nan for a_, b_ in r]
            a.plot(PSIG, y, marker="os^v<>Dp*h"[i], ls="--" if v.startswith("eta") and mode == "rf" else "-",
                   label=lab)
        a.axhline(1.0, c="k", lw=1)
        a.set(ylim=(0.5, 3.0), xlabel="pressure (psig)", ylabel="Ar / He threshold ratio",
              title="Laser-seeded: every variant keeps Ar > He" if mode == "laser" else
              "RF-only: photon feedback η ≥ 0.3 crosses 1")
        a.grid(alpha=0.3)
        a.legend(fontsize=6, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(RES / "fig_mechanisms.png", dpi=150)
    plt.close(fig)


def fig_dense(plt):
    d = json.loads((RES / "dense_regime.json").read_text())
    data = {}
    for r in d["results"]:
        if r["x"] > 0 and not r["ei"]:
            continue
        en = [row["E_N_Td"] for row in r["rows"] if row["converged"] and row["kiz"] > 0]
        k = [row["kiz"] for row in r["rows"] if row["converged"] and row["kiz"] > 0]
        data[(r["gas"], r["x"])] = dict(zip(en, k))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    xs = [0.0, 1e-5, 1e-4, 1e-3, 1e-2]
    cm = plt.get_cmap("viridis")
    for i, x in enumerate(xs):
        c = cm(i / (len(xs) - 1))
        for gas, ls in (("He", "--"), ("Ar", "-")):
            dd = data.get((gas, x), {})
            e = sorted(dd)
            ax[0].plot(e, [dd[v] for v in e], ls, c=c, marker="o" if gas == "Ar" else "s", ms=3,
                       label=f"{gas}, x = {x:g}")
        he, ar = data.get(("He", x), {}), data.get(("Ar", x), {})
        common = sorted(set(he) & set(ar))
        ax[1].plot(common, [ar[v] / he[v] for v in common], "-o", c=c, ms=4, label=f"x = n_e/N = {x:g}")
    ax[0].set(xscale="log", yscale="log", ylim=(1e-28, 1e-16), xlabel="E_eff/N (Td)",
              ylabel="ionization rate coefficient k_iz (m³/s)",
              title="Two-term Boltzmann + e-e + e-ion (dashed He, solid Ar)")
    ax[0].legend(fontsize=6, ncol=2)
    ax[0].grid(alpha=0.3, which="both")
    ax[1].axhline(1, c="k", lw=1)
    ax[1].set(xscale="log", yscale="log", xlabel="E_eff/N (Td)", ylabel="k_iz(Ar) / k_iz(He)",
              title="Dilute: He ionizes faster; dense (x ≥ 1e-4, low E/N): Ar does")
    ax[1].legend(fontsize=7)
    ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(RES / "fig_dense_regime.png", dpi=150)
    plt.close(fig)


def fig_spark(plt):
    s = json.loads((RES / "spark.json").read_text())
    t = _thr()

    def get(gas, x0, kset="nominal", drp="nominal"):
        m = {float(r["psig"]): r["thr_kV_cm"] for r in s
             if r["gas"] == gas and r["x0"] == x0 and r["kset"] == kset and r["dr_product"] == drp}
        return [m.get(p, np.nan) if m.get(p) is not None else np.nan for p in PSIG]

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    a = ax[0]
    a.plot(PSIG, get("He", 0.0), "-o", c=COL["He"], label="He, multiphoton seed (no spark)")
    a.plot(PSIG, get("Ar", 0.0), ":o", c=COL["Ar"], label="Ar, multiphoton seed")
    a.plot(PSIG, get("Ar", 1e-4), "-s", c=COL["Ar"], label="Ar, dense spark x0 = 1e-4 (hypothesized)")
    lo, hi = np.array(get("Ar", 1e-4, "low")), np.array(get("Ar", 1e-4, "high"))
    a.fill_between(PSIG, np.fmin(lo, hi), np.fmax(lo, hi), color=COL["Ar"], alpha=0.2, label="Ar spark, low-high rates")
    a.plot(PSIG, get("He", 1e-4), "--s", c=COL["He"], label="He, same dense seed (not observed)")
    a.set(yscale="log", xlabel="pressure (psig)", ylabel="threshold peak field (kV/cm)",
          title="Seed to N_br = 1e13; field assumed inside the plasma")
    a.legend(fontsize=7)
    a.grid(alpha=0.3, which="both")
    a = ax[1]
    for lab, y in (("Ar spark / He multiphoton (nominal)", np.array(get("Ar", 1e-4)) / np.array(get("He", 0.0))),
                   ("  low rates", np.array(get("Ar", 1e-4, "low")) / np.array(get("He", 0.0, "low"))),
                   ("  high rates", np.array(get("Ar", 1e-4, "high")) / np.array(get("He", 0.0, "high"))),
                   ("  recombination -> Ar2*", np.array(get("Ar", 1e-4, "nominal", "Ar2t")) / np.array(get("He", 0.0))),
                   ("  recombination -> ground", np.array(get("Ar", 1e-4, "nominal", "ground")) / np.array(get("He", 0.0))),
                   ("same seed x0 = 1e-4 (Ar / He)", np.array(get("Ar", 1e-4)) / np.array(get("He", 1e-4))),
                   ("multiphoton seed, both", np.array(get("Ar", 0.0)) / np.array(get("He", 0.0)))):
        a.plot(PSIG, y, "-o", ms=4, label=lab)
    a.axhline(1, c="k", lw=1)
    a.set(xlabel="pressure (psig)", ylabel="Ar / He threshold ratio", ylim=(0.3, 3),
          title="Only a dense Ar spark brings Ar to or below He")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(RES / "fig_spark.png", dpi=150)
    plt.close(fig)


def fig_delays(plt):
    p = RES / "delays.json"
    if not p.exists():
        return
    d = [x for x in json.loads(p.read_text()) if "mc" in x]
    psig = 200.0
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.2))
    cases = ((("laser", "nominal"), "o", "-", "laser, multiphoton seed"),
             (("laser", "spark1e-4"), "s", "--", "laser, dense spark x0=1e-4"),
             (("rf", "nominal"), "^", ":", "RF-only"),
             (("rf", "eta1"), "v", "-.", "RF-only + photon feedback η=1"))
    for (mode, var), mk, ls, lab0 in cases:
        for gas in ("He", "Ar", "N2"):
            if var == "eta1" and gas == "N2":
                continue                           # N2 has no photoemitting reservoir: identical to RF-only
            rows = sorted([x for x in d if x["gas"] == gas and x["psig"] == psig and x["mode"] == mode
                           and x["variant"] == var], key=lambda x: x["E_kV_cm"])
            if not rows or all(x["mc"]["p_breakdown"] == 0 for x in rows):
                continue
            e = [x["E_kV_cm"] for x in rows]
            kw = dict(c=COL[gas], marker=mk, ls=ls, ms=5)
            ax[0].plot(e, [x["mc"]["p_breakdown"] for x in rows], label=f"{gas}: {lab0}", **kw)
            ax[1].plot(e, [np.nan if x["mc"]["delay_median_ns"] is None else x["mc"]["delay_median_ns"]
                           for x in rows], **kw)
            ax[2].plot(e, [np.nan if x["mc"]["jitter_p90_p10_ns"] is None else x["mc"]["jitter_p90_p10_ns"]
                           for x in rows], **kw)
    ax[0].set(xscale="log", xlabel="peak field (kV/cm)", ylabel="breakdown probability per shot",
              title="200 psig: probability")
    ax[1].set(xscale="log", yscale="log", xlabel="peak field (kV/cm)", ylabel="median delay (ns)",
              title="delay after laser (RF-only: after first electron)")
    ax[2].set(xscale="log", yscale="log", xlabel="peak field (kV/cm)", ylabel="jitter p90 - p10 (ns)",
              title="2000-shot Monte Carlo (RF power 2 %, laser 3 %)")
    for a in ax:
        a.grid(alpha=0.3, which="both")
        a.set_xticks([30, 50, 100, 150, 200, 300])
        a.set_xticklabels(["30", "50", "100", "150", "200", "300"])
    fig.legend(*ax[0].get_legend_handles_labels(), loc="lower center", ncol=4, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(RES / "fig_delays.png", dpi=150)
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig_mechanisms(plt)
    fig_dense(plt)
    fig_spark(plt)
    fig_delays(plt)


if __name__ == "__main__":
    main()
