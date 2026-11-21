"""Figures and tables for the full-physics (CR + transport + gas heating + laser propagation) runs.

    python -m rf_plasma.full_physics_figures      (writes results/coupled3d/full_physics_He_1MPa/figures/)
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .coupled3d import RESULTS

BASE = RESULTS / "laser_trigger_He_1MPa"
FULL = RESULTS / "full_physics_He_1MPa"
FIG = FULL / "figures"
CASES = {  # label: (base dir, full dir, nocr dir)
    "guide RF only, 8 MV/m": ("guide_He_rf-only_E8MVm", "guide_He_rf-only_E8MVm_full", "guide_He_rf-only_E8MVm_nocr"),
    "guide + 10 mJ/3 ns, 8 MV/m": ("guide_He_rf+laser_3000ps_E8MVm", "guide_He_rf+laser_3000ps_E8MVm_full",
                                   "guide_He_rf+laser_3000ps_E8MVm_nocr"),
    "guide + 10 mJ/3 ns, 6 MV/m": ("guide_He_rf+laser_3000ps_E6MVm", "guide_He_rf+laser_3000ps_E6MVm_full", None),
    "cavity + 10 mJ/3 ns, 8 MV/m": ("cavity_He_rf+laser_3000ps_E8MVm", "cavity_He_rf+laser_3000ps_E8MVm_full", None),
}


def _load(root, name):
    if name is None or not (root / name / "summary.json").exists():
        return None
    d = root / name
    return {"s": json.loads((d / "summary.json").read_text()), "ts": np.load(d / "timeseries.npz"),
            "rf": np.load(d / "rf_series.npz") if (d / "rf_series.npz").exists() else None, "dir": d}


def _gamma(run):
    r = run["rf"]
    g = np.abs(r["bwd"]) / np.maximum(np.abs(r["fwd"]), 1e-30)
    return r["t_phasor"], g


def histories(runs):
    fig, ax = plt.subplots(len(runs), 4, figsize=(18, 3.3 * len(runs)), squeeze=False)
    sty = {"base": ("k", "-", "base (original kernels)"), "nocr": ("C0", "--", "full minus CR"),
           "full": ("C3", "-", "full physics")}
    for i, (lab, rr) in enumerate(runs.items()):
        for key, run in rr.items():
            if run is None:
                continue
            c, ls, nm = sty[key]
            ts = run["ts"]
            t = ts["t"] * 1e9
            if run["rf"] is not None:
                tg, g = _gamma(run)
                ax[i, 0].plot(tg * 1e9, g, color=c, ls=ls, label=nm)
            ax[i, 1].semilogy(t, np.maximum(ts["bulk_ne_max"], 1e8), color=c, ls=ls)
            ax[i, 2].semilogy(t, np.maximum(ts["fil_ne_axis_focus"], 1e8), color=c, ls=ls)
            ax[i, 3].plot(t, ts["fil_te_axis_focus"], color=c, ls=ls)
        ax[i, 0].set(ylabel=f"{lab}\n|Gamma| at feed", ylim=(0, 1.05))
        ax[i, 1].set(ylabel="bulk n_e max (m^-3)")
        ax[i, 2].set(ylabel="filament n_e on axis @focus (m^-3)")
        ax[i, 3].set(ylabel="filament T_e on axis @focus (eV)")
        ax[i, 0].legend(fontsize=8)
    for a in ax[-1]:
        a.set_xlabel("t (ns)")
    fig.suptitle("He, 1 MPa, WR-284: base kernels vs full physics (CR + transport + gas heating + propagation)")
    fig.tight_layout()
    fig.savefig(FIG / "fp_histories.png", dpi=110)
    plt.close(fig)


def cr_and_heating(runs):
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.4))
    for k, (lab, rr) in enumerate(runs.items()):
        run = rr.get("full")
        if run is None:
            continue
        ts = run["ts"]
        t = ts["t"] * 1e9
        c = f"C{k}"
        ax[0].semilogy(t, np.maximum(ts["bulk_hestar_max"], 1e8), color=c, label=f"{lab}: He* (n=2)")
        ax[0].semilogy(t, np.maximum(ts["bulk_he2star_max"], 1e8), color=c, ls="--", label="He2*")
        ax[1].semilogy(t, np.maximum(ts["fil_hestar_axis_focus"], 1e8), color=c, label=f"{lab}: He*")
        ax[1].semilogy(t, np.maximum(ts["fil_alpha_axis_focus"], 1e8), color=c, ls=":", label="He2+ (alpha)")
        ax[2].semilogy(t, ts["fil_tg_axis_focus"], color=c, label=f"{lab}: filament axis")
        ax[2].semilogy(t, ts["bulk_tg_max"], color=c, ls="--", label="bulk max")
    ax[0].set(title="Bulk (3D grid): excited species maxima", xlabel="t (ns)", ylabel="density (m^-3)", ylim=(1e12, None))
    ax[1].set(title="Filament axis at focus: He* and He2+", xlabel="t (ns)", ylabel="density (m^-3)", ylim=(1e12, None))
    ax[2].set(title="Gas temperature", xlabel="t (ns)", ylabel="T_g (K)")
    for a in ax:
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "fp_cr_heating.png", dpi=110)
    plt.close(fig)


def beam(runs):
    run = runs["guide + 10 mJ/3 ns, 8 MV/m"].get("full")
    if run is None:
        return
    ts = run["ts"]
    s = run["s"]["case"]
    t = ts["t"]
    y = ts["fil_y"] * 1e3
    lam, w0 = s["wavelength_m"], s["waist_m"]
    zr = np.pi * w0 ** 2 / lam
    yf = 0.5 * (y.min() + y.max())
    on = np.where(ts["beam_r2nd"].max(axis=1) > 0)[0]
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.4))
    yy = np.linspace(y.min(), y.max(), 300)
    w_vac = w0 * np.sqrt(1 + ((yy - yf) * 1e-3 / zr) ** 2)
    ax[0].plot(yy, w_vac * 1e6, "k--", label="vacuum Gaussian w(y)")
    for tt in (4e-9, 5e-9, 6.5e-9, 8e-9):
        i = on[np.argmin(abs(t[on] - tt))] if on.size else None
        if i is None:
            continue
        m = ts["beam_r2nd"][i] > 0
        ax[0].plot(y[m], ts["beam_r2nd"][i][m] * 1e6, "o-", ms=3, label=f"t = {t[i] * 1e9:.1f} ns")
        ax[1].semilogy(y[m], ts["beam_axis_I"][i][m], "o-", ms=3, label=f"t = {t[i] * 1e9:.1f} ns")
    ax[0].set(xlabel="position along beam (mm, entry -> exit)", ylabel="beam radius sqrt(2<r^2>) (um)",
              title="Laser propagation through the filament", yscale="log")
    ax[1].set(xlabel="position along beam (mm)", ylabel="on-axis intensity (W/m^2)", title="On-axis intensity")
    ax[2].plot(t * 1e9, ts["beam_r2nd_focus"] * 1e6, label="radius at focus (um)")
    a2 = ax[2].twinx()
    a2.semilogy(t * 1e9, np.maximum(ts["beam_axis_I_focus"], 1e6), color="C3", label="on-axis I at focus")
    ax[2].axhline(w0 * 1e6, color="k", ls="--", lw=0.8, label="vacuum waist")
    ax[2].set(xlabel="t (ns)", ylabel="radius at focus (um)", title="Focal spot vs time")
    a2.set_ylabel("I (W/m^2)", color="C3")
    for a in ax[:2]:
        a.legend(fontsize=8)
    ax[2].legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fp_beam.png", dpi=110)
    plt.close(fig)


def maps(runs):
    fig, ax = plt.subplots(2, 2, figsize=(15, 6.5))
    for k, lab in enumerate(("guide RF only, 8 MV/m", "guide + 10 mJ/3 ns, 8 MV/m")):
        run = runs[lab].get("full")
        if run is None:
            continue
        fr = np.load(run["dir"] / "frames.npz")
        i = len(fr["t"]) - 1
        ext = (0, fr["tg_yz"].shape[2] * float(fr["dx"]) * 1e3, 0, fr["tg_yz"].shape[1] * float(fr["dx"]) * 1e3)
        im = ax[0, k].imshow(np.log10(np.maximum(fr["hestar_yz"][i], 1e10)), origin="lower", aspect="auto",
                             extent=ext, cmap="viridis", vmin=12)
        fig.colorbar(im, ax=ax[0, k], label="log10 He* (m^-3)")
        im = ax[1, k].imshow(fr["tg_yz"][i], origin="lower", aspect="auto", extent=ext, cmap="inferno")
        fig.colorbar(im, ax=ax[1, k], label="T_g (K)")
        ax[0, k].set(title=f"{lab}: He* at t = {fr['t'][i] * 1e9:.1f} ns (y-z mid-plane)", ylabel="y (mm)")
        ax[1, k].set(title="bulk gas temperature", xlabel="z (mm, along guide)", ylabel="y (mm)")
    fig.tight_layout()
    fig.savefig(FIG / "fp_maps.png", dpi=110)
    plt.close(fig)


def table(runs):
    rows = ["| case | physics | Gamma_end | bulk n_e max (end) | fil. n_e max | fil. T_e max (eV) | fil. T_g max (K) "
            "| bulk T_g max (K) | bulk He* max | laser absorbed (J) | filament closure | bulk closure |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    out = {}
    for lab, rr in runs.items():
        for key in ("base", "nocr", "full"):
            run = rr.get(key)
            if run is None:
                continue
            s = run["s"]
            p = s["peaks"]
            rec = {"gamma_end": s["rf"].get("gamma_end"), "bulk_ne_end": p["bulk_ne_max_end"],
                   "fil_ne_max": p["fil_ne_axis_max"], "fil_te_max": p["fil_te_axis_max"],
                   "fil_tg_max": p["fil_tg_axis_max"], "bulk_tg_max": p.get("bulk_tg_max_max"),
                   "bulk_hestar_max": p.get("bulk_hestar_max_max"), "laser_abs_J": s["laser"]["absorbed_J"],
                   "fil_closure": s["ledger_filament_J"].get("closure_residual_rel"),
                   "bulk_closure": s["ledger_bulk_J"].get("closure_residual_rel")}
            out[f"{lab} | {key}"] = rec
            f = lambda v, fmt: "-" if v is None else format(v, fmt)  # noqa: E731
            rows.append(f"| {lab} | {key} | {f(rec['gamma_end'], '.3f')} | {f(rec['bulk_ne_end'], '.2e')} | "
                        f"{f(rec['fil_ne_max'], '.2e')} | {f(rec['fil_te_max'], '.1f')} | {f(rec['fil_tg_max'], '.0f')} | "
                        f"{f(rec['bulk_tg_max'], '.1f')} | {f(rec['bulk_hestar_max'], '.2e')} | "
                        f"{f(rec['laser_abs_J'], '.2e')} | {f(rec['fil_closure'], '.1e')} | {f(rec['bulk_closure'], '.1e')} |")
    (FIG / "fp_table.md").write_text("\n".join(rows) + "\n")
    (FIG / "fp_table.json").write_text(json.dumps(out, indent=1))
    print("\n".join(rows))


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    runs = {lab: {"base": _load(BASE, b), "full": _load(FULL, f), "nocr": _load(FULL, n)}
            for lab, (b, f, n) in CASES.items()}
    histories(runs)
    cr_and_heating(runs)
    beam(runs)
    maps(runs)
    table(runs)


if __name__ == "__main__":
    main()
