"""Figures and comparison tables for the coupled 3D runs (reads results/coupled3d/<case>/)."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
import numpy as np  # noqa: E402

from .coupled3d import RESULTS  # noqa: E402

GAS_COLOR = {"Ar": "tab:blue", "He": "tab:orange", "N2": "tab:green"}


def _load(case_dir: Path):
    s = json.loads((case_dir / "summary.json").read_text(encoding="utf-8"))
    ts = dict(np.load(case_dir / "timeseries.npz"))
    rf = dict(np.load(case_dir / "rf_series.npz")) if (case_dir / "rf_series.npz").exists() else None
    sn = dict(np.load(case_dir / "snapshots.npz"))
    return s, ts, rf, sn


def case_figures(case_dir: Path):
    s, ts, rf, sn = _load(case_dir)
    t = ts["t"] * 1e9
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    a = ax[0, 0]
    a.semilogy(t, np.maximum(ts["fil_ne_axis_focus"], 1e6), label="filament axis (focus)")
    a.semilogy(t, np.maximum(ts["bulk_ne_max"], 1e6), label="bulk max")
    a.set_xlabel("t (ns)"); a.set_ylabel("n_e (m$^{-3}$)"); a.legend(); a.set_title("Electron density")
    a = ax[0, 1]
    a.plot(t, ts["fil_te_axis_focus"], label="T$_e$ axis (eV)")
    a2 = a.twinx(); a2.plot(t, ts["fil_tg_axis_focus"], "r", label="T$_g$ axis (K)"); a2.set_ylabel("T$_g$ (K)", color="r")
    a.set_xlabel("t (ns)"); a.set_ylabel("T$_e$ (eV)"); a.set_title("Filament axis temperatures")
    a = ax[0, 2]
    a.plot(t, ts["fil_r50_focus"] * 1e6, label="r$_{50}$ electrons")
    a.plot(t, (ts["fil_rb_out_focus"] - ts["fil_rb_out_focus"][0]) * 1e6, label="outer ring displacement")
    a.set_xlabel("t (ns)"); a.set_ylabel("µm"); a.legend(); a.set_title("Filament radius / expansion")
    a = ax[1, 0]
    if rf is not None:
        tp = rf["t_phasor"] * 1e9
        a.plot(tp, np.abs(rf["bwd"] / rf["fwd"]), label="|Γ| at port")
        if np.any(rf["p_out"]):
            pf = np.abs(rf["fwd"]) ** 2
            a.plot(tp, rf["p_out"] / s["rf"]["P_forward_W"], label="P$_{out}$/P$_{fwd}$")
        a.set_ylabel("ratio"); a.legend()
    a.set_xlabel("t (ns)"); a.set_title("RF response")
    a = ax[1, 1]
    if rf is not None:
        tt = np.arange(len(rf["ey_focus"])) * float(rf["dt"]) * 1e9
        a.plot(tt, rf["ey_focus"] / 1e6, lw=0.3)
        a.set_ylabel("E$_y$ focus edge (MV/m)")
    a.set_xlabel("t (ns)"); a.set_title("Field at the focus edge")
    a = ax[1, 2]
    a.plot(t, ts["laser_absorbed"] * 1e3, label="laser absorbed (mJ)")
    a.plot(t, (ts["rf_to_fil"] + ts["rf_to_nodes"]) * 1e3, label="RF absorbed (mJ)")
    a.set_xlabel("t (ns)"); a.set_ylabel("mJ"); a.legend(); a.set_title("Energy deposition")
    fig.suptitle(s["name"] + "  —  " + s["layout"], fontsize=9)
    fig.tight_layout()
    fig.savefig(case_dir / "timeseries.png", dpi=110)
    plt.close(fig)

    # filament space-time map and snapshot maps
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    y = ts["fil_y"] * 1e3
    m = np.maximum(ts["fil_line_e"].T, 1e-3)
    im = ax[0].pcolormesh(t, y, m, norm=LogNorm(vmin=max(m.max() * 1e-8, 1e-3), vmax=m.max()), shading="auto")
    fig.colorbar(im, ax=ax[0], label="electrons per metre")
    ax[0].set_xlabel("t (ns)"); ax[0].set_ylabel("y along laser (mm)"); ax[0].set_title("Filament line density")
    tags = sorted({k.split("_")[0] for k in sn if k.startswith("s") and k.endswith("_t")}, key=lambda q: int(q[1:]))
    if tags:
        last = tags[-1]
        ne = np.maximum(sn[f"{last}_ne_yz"], 1e6)
        im = ax[1].imshow(ne.T, origin="lower", aspect="auto", norm=LogNorm(vmin=max(ne.max() * 1e-6, 1e6), vmax=ne.max()),
                          extent=[0, ne.shape[0], 0, ne.shape[1]])
        fig.colorbar(im, ax=ax[1], label="bulk n$_e$ (m$^{-3}$)")
        ax[1].set_xlabel("y (cells)"); ax[1].set_ylabel("z (cells)")
        ax[1].set_title(f"bulk n$_e$, y-z plane through focus, t={float(sn[last + '_t']) * 1e9:.1f} ns")
    fig.tight_layout()
    fig.savefig(case_dir / "maps.png", dpi=110)
    plt.close(fig)

    # filament radial profiles at the focus edge for each snapshot
    if tags:
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        e_f = int(np.argmin(np.abs(ts["fil_y"] - ts["fil_y"].mean())))
        for tg in tags:
            rb = sn[f"{tg}_fil_rb"][e_f]
            rc = 0.5 * (rb[1:] + rb[:-1]) * 1e6
            lab = f"{float(sn[tg + '_t']) * 1e9:.2f} ns"
            ax[0].semilogy(rc, np.maximum(sn[f"{tg}_fil_ne"][e_f], 1e6), label=lab)
            ax[1].plot(rc, sn[f"{tg}_fil_te"][e_f], label=lab)
            ax[2].plot(rc, sn[f"{tg}_fil_tg"][e_f], label=lab)
        for a_, yl in zip(ax, ("n$_e$ (m$^{-3}$)", "T$_e$ (eV)", "T$_g$ (K)")):
            a_.set_xscale("log"); a_.set_xlabel("r (µm)"); a_.set_ylabel(yl); a_.legend(fontsize=7)
        fig.suptitle(f"{s['name']}: radial profiles at the focus edge")
        fig.tight_layout()
        fig.savefig(case_dir / "radial.png", dpi=110)
        plt.close(fig)


def matrix_table(root: Path = RESULTS) -> list:
    rows = []
    for d in sorted(p for p in root.iterdir() if (p / "summary.json").exists()):
        s = json.loads((d / "summary.json").read_text(encoding="utf-8"))
        c = s["case"]
        rows.append(dict(
            name=s["name"], geometry=c["geometry"], gas=c["gas"], rf=c["rf"], laser=c["laser"],
            fwhm_ps=c["fwhm_s"] * 1e12 if c["laser"] else None, energy_mJ=c["energy_j"] * 1e3 if c["laser"] else None,
            status=s["status"],
            fil_ne_max=s["peaks"]["fil_ne_axis_max"], fil_te_max=s["peaks"]["fil_te_axis_max"],
            fil_tg_max=s["peaks"]["fil_tg_axis_max"], fil_electrons_end=s["peaks"]["fil_electrons_end"],
            bulk_ne_max_end=s["peaks"]["bulk_ne_max_end"],
            laser_abs_frac=s["laser"]["fraction_absorbed"],
            rf_absorbed_mJ=(s["rf"].get("rf_absorbed_J", 0.0) * 1e3) if s["rf"] else 0.0,
            gamma_start=s["rf"].get("gamma_start") if s["rf"] else None,
            gamma_end=s["rf"].get("gamma_end") if s["rf"] else None,
            e_focus_end_MV_m=(s["rf"].get("focus_peak_end_V_m") or 0) / 1e6 if s["rf"] else None,
            e_focus_min_MV_m=(s["rf"].get("focus_peak_min_V_m") or 0) / 1e6 if s["rf"] else None,
            delta_em_mJ=s["rf"].get("delta_em_J", 0.0) * 1e3 if s["rf"] else None,
            rf_balance=s["rf"].get("balance_residual") if s["rf"] else None,
            fil_closure=s["ledger_filament_J"].get("closure_residual_rel"),
            tag=c.get("tag", ""), ee_scale=c.get("ee_scale"), block=c.get("block"), rings=c.get("rings"),
            wall_s=s["wall_s"]))
    (root / "matrix_summary.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    return rows


def _is_matrix(r):
    # sensitivity runs carry an extra tag suffix (e.g. E30mJ_ee0.1); matrix tags are '' or E<n>mJ
    t = r.get("tag") or ""
    return t == "" or (t.startswith("E") and t.endswith("mJ") and "_" not in t)


def threshold_figure(rows: list, root: Path = RESULTS):
    lo = [r for r in rows if r["laser"] and not r["rf"] and r["geometry"] == "taper" and _is_matrix(r)]
    if not lo:
        return
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for fw, ls in ((300.0, "-o"), (3000.0, "--s")):
        for gas in ("Ar", "N2", "He"):
            pts = sorted((r["energy_mJ"], r["fil_ne_max"], r["laser_abs_frac"]) for r in lo
                         if r["gas"] == gas and abs(r["fwhm_ps"] - fw) < 1)
            if not pts:
                continue
            e, n, a = zip(*pts)
            ax[0].plot(e, n, ls, color=GAS_COLOR[gas], label=f"{gas} {fw:g} ps")
            ax[1].plot(e, np.maximum(a, 1e-18), ls, color=GAS_COLOR[gas], label=f"{gas} {fw:g} ps")
    ax[0].axhline(3.9e27, color="0.5", lw=0.8); ax[0].text(1.1, 2.5e27, "532 nm critical density", fontsize=7)
    ax[0].set_ylabel("peak on-axis n$_e$ (m$^{-3}$)"); ax[1].set_ylabel("laser fraction absorbed")
    for a in ax:
        a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("laser energy at focus (mJ)"); a.grid(alpha=0.3)
    ax[0].legend(fontsize=7, ncol=2)
    ax[0].set_title("Laser-only breakdown, 790.8 kPa, w$_0$ = 20 $\\mu$m (no RF)")
    fig.tight_layout()
    fig.savefig(root / "laser_thresholds.png", dpi=120)
    plt.close(fig)


def matrix_figure(rows: list, root: Path = RESULTS):
    comb = [r for r in rows if r["rf"] and _is_matrix(r)]
    if not comb:
        return
    fig, ax = plt.subplots(1, 5, figsize=(22, 4.4))
    labels = [f"{r['geometry'][:3]} {r['gas']} " + (f"{r['fwhm_ps']:.0f}ps {r['energy_mJ']:g}mJ" if r["laser"]
                                                    else "RF only") for r in comb]
    x = np.arange(len(comb))
    col = [GAS_COLOR[r["gas"]] for r in comb]
    ax[0].bar(x, [r["fil_ne_max"] for r in comb], color=col); ax[0].set_yscale("log"); ax[0].set_title("peak filament n$_e$ (m$^{-3}$)")
    ax[1].bar(x, [r["laser_abs_frac"] for r in comb], color=col); ax[1].set_title("laser fraction absorbed")
    ax[2].bar(x, [max(r["rf_absorbed_mJ"], 1e-9) for r in comb], color=col); ax[2].set_yscale("log")
    ax[2].set_title("RF energy absorbed (mJ)")
    ax[3].bar(x - 0.2, [r["gamma_start"] for r in comb], 0.4, color="0.6", label="start")
    ax[3].bar(x + 0.2, [r["gamma_end"] for r in comb], 0.4, color=col, label="end"); ax[3].set_title("|Γ| at port"); ax[3].legend()
    ax[4].bar(x, [r["e_focus_min_MV_m"] for r in comb], color=col); ax[4].axhline(10, color="0.5", lw=0.8)
    ax[4].set_title("min RF peak |E$_y$| at focus (MV/m)")
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    fig.tight_layout()
    fig.savefig(root / "matrix_overview.png", dpi=120)
    plt.close(fig)


def sensitivity_figure(rows: list, root: Path = RESULTS):
    """Relative change of key outputs for tagged sensitivity runs against their untagged reference."""
    by_name = {r["name"]: r for r in rows}
    keys = (("fil_ne_max", "peak filament n$_e$"), ("rf_absorbed_mJ", "RF absorbed"),
            ("gamma_end", "|Γ| end"), ("e_focus_min_MV_m", "min focus |E|"), ("laser_abs_frac", "laser absorbed"))
    pairs = []
    for r in rows:
        if _is_matrix(r):
            continue
        ref = by_name.get(r["name"].rsplit("_", 1)[0])
        if ref is not None and ref is not r:
            pairs.append((r["name"].rsplit("_", 1)[1], ref, r))
    if not pairs:
        return
    fig, ax = plt.subplots(figsize=(9, 4.2))
    w = 0.8 / len(keys)
    x = np.arange(len(pairs))
    for i, (k, lab) in enumerate(keys):
        d = [100 * (r[k] - ref[k]) / ref[k] if ref[k] else 0.0 for _, ref, r in pairs]
        ax.bar(x + (i - len(keys) / 2 + 0.5) * w, d, w, label=lab)
    ax.axhspan(-5, 5, color="0.9", zorder=0, label="±5 % gate")
    ax.set_xticks(x); ax.set_xticklabels([f"{ref['gas']} {lab}" for lab, ref, _ in pairs], fontsize=8)
    ax.set_ylabel("change vs reference (%)"); ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3, axis="y")
    ax.set_title("Sensitivity of coupled cavity runs (reference: ν$_{ee}$ scale 1, block 8, 48 rings, seed 1e9 m$^{-3}$)")
    fig.tight_layout()
    fig.savefig(root / "sensitivity.png", dpi=120)
    plt.close(fig)


def main():
    for d in sorted(p for p in RESULTS.iterdir() if (p / "summary.json").exists()):
        case_figures(d)
    rows = matrix_table()
    matrix_figure(rows)
    threshold_figure(rows)
    sensitivity_figure(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
