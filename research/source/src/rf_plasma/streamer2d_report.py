"""Aggregate the streamer2d case runs into results/streamer2d/summary.json and figures.

Post-processing only; it reads results/streamer2d/<run>/{summary.json, history.npz, snap_*.npz}.

For each physical case it prefers the +/-8 mm axial-domain rerun ("_z8") when it exists and
otherwise uses the +/-4 mm run.  Where both exist, the overlap gives a domain-size sensitivity.

Usage:  python -m rf_plasma.streamer2d_report
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from rf_plasma.streamer2d import CASES, RESULTS

# ---------------------------------------------------------------- measured reference values
MEASURED = {
    "air_filament_diameter_m": {
        "range": [0.2e-3, 0.6e-3],
        "gas": "air",
        "source": "Technion (Krasik/Leopold group) filament imaging; literature/experimental-data-krasik-leopold.md",
        "note": "AIR, not helium. Used here only as an order-of-magnitude reference.",
    },
    "air_growth_speed_along_E_m_s": {
        "value": 5e5,
        "gas": "air",
        "source": "~5e7 cm/s, Technion; literature/experimental-data-krasik-leopold.md",
        "note": "AIR, not helium.",
    },
    "He_ne_peak_m3": {
        "value": 6e22,
        "gas": "He 0.2 MPa",
        "source": "Stark broadening of He I 388.9 nm at ~300 kW output, Beilin 2015 Fig. 3 / Yatom 2016 Fig. 16",
        "note": "Line-of-sight average through a 5 mm hole, 0.25 ns resolution, laser direction mixed; the model value "
                "is local (on axis at z=0).",
    },
    "He_delta_delay_s": {
        "value": 30e-9, "uncertainty": 15e-9,
        "gas": "He 0.2 MPa",
        "definition": "delay(16.6 kV/cm) - delay(22.3 kV/cm), laser entry -> RF output switching",
        "source": "Beilin 2015; Yatom 2016 sec. 4.2 / Fig. 14",
        "note": "The measured delay covers laser entry to output-pulse switching (tens of ns).  The model starts from "
                "an already-formed laser channel and runs <=10 ns, so only a proxy is computed (see delay_proxies).",
    },
}

# physical case -> run directories in order of preference
PREFER = ["_z8", ""]
REPORT_TIMES = (0.0, 1e-9, 2e-9, 5e-9, 10e-9)
D_THRESH = 0.2e-3  # m: lower end of the measured (air) diameter range, used for the delay proxy


def _load_run(name):
    d = RESULTS / name
    if not (d / "summary.json").exists():
        return None
    s = json.loads((d / "summary.json").read_text())
    h = dict(np.load(d / "history.npz"))
    snaps = []
    for f in sorted(d.glob("snap_t*.npz")):
        z = np.load(f)
        rc, zc, ne = z["rc"], z["zc"], z["ne"].astype(float)
        j0 = int(np.argmin(np.abs(zc)))
        snaps.append({"t_s": float(z["t_s"]), "ne_axis_mid": float(ne[0, j0]),
                      "ne_axis_mid_avg_r100um": _area_avg(rc, z["rf"], ne[:, j0], 100e-6),
                      "file": f.name})
    snaps.sort(key=lambda x: x["t_s"])
    return {"name": name, "summary": s, "hist": h, "snaps": snaps}


def _area_avg(rc, rf, col, R):
    """Area-weighted average of a radial profile over r <= R (cell-face based)."""
    a = np.pi * (np.minimum(rf[1:], R) ** 2 - np.minimum(rf[:-1], R) ** 2)
    a = np.clip(a, 0, None)
    return float(np.sum(a * col) / np.sum(a))


def _at(h, key, t):
    tt = h["t_s"]
    if t > tt[-1] + 1e-15:
        return None
    return float(np.interp(t, tt, h[key]))


def _first_crossing(t, y, thr):
    idx = np.nonzero(np.asarray(y) >= thr)[0]
    if idx.size == 0:
        return None
    i = idx[0]
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (thr - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1]))


def _rate_fit(t, y, t_min):
    m = (t >= t_min) & np.isfinite(y)
    if m.sum() < 5:
        return None
    return float(np.polyfit(t[m], y[m], 1)[0])


def case_metrics(run):
    s, h = run["summary"], run["hist"]
    t = h["t_s"]
    r_mid, r_max, r_sig = h["r_1e20_mid_m"], h["r_1e20_max_m"], h["r_sigma_mid_m"]
    G = h["line_conductance_mid_S_m"]
    t_d_mid = _first_crossing(t, 2 * r_mid, D_THRESH)
    t_d_max = _first_crossing(t, 2 * r_max, D_THRESH)
    t_G2 = _first_crossing(t, G / G[0], 2.0)
    # linear extrapolation of the midplane 1e20 radius (fit over the last half of the run)
    vr_mid = _rate_fit(t, r_mid, 0.5 * t[-1])
    t_d_mid_extrap = None
    if t_d_mid is None and vr_mid and vr_mid > 0:
        t_d_mid_extrap = float(t[-1] + (D_THRESH / 2 - r_mid[-1]) / vr_mid)
    tv = s["results"]["tip_velocity_1e20"]
    # late-time tip speed (last 1 ns) to show acceleration
    tl = max(t[-1] - 1e-9, 1e-9)
    vlate = None
    if t[-1] > 1.5e-9:
        vp, vm = _rate_fit(t, h["z_plus_1e20_m"], tl), _rate_fit(t, -h["z_minus_1e20_m"], tl)
        if vp is not None and vm is not None:
            vlate = 0.5 * (vp + vm)

    def at_times(key):
        d = {f"{x*1e9:g}ns": _at(h, key, x) for x in REPORT_TIMES}
        d["final"] = float(h[key][-1])
        return d

    return {
        "run": run["name"],
        "pressure_MPa": s["config"]["pressure_Pa"] / 1e6,
        "E0_kV_cm": s["config"]["E0_V_m"] / 1e5,
        "E0_over_N_Td": s["derived"]["E0_over_N_Td"],
        "seed": {"n0_m3": s["config"]["n0_m3"], "a_m": s["config"]["a_m"], "L_m": s["config"]["L_m"],
                 "n_bg_m3": s["config"]["n_bg_m3"]},
        "domain_z_half_m": s["derived"]["grid"]["z_half_m"],
        "t_final_s": s["status"]["t_final_s"],
        "stop_reason": s["status"]["stop_reason"],
        "wall_s": s["status"]["wall_s"],
        "steps": s["status"]["steps"],
        "max_rel_charge_residual": s["status"]["max_rel_charge_residual"],
        "radius_1e20_mid_m": at_times("r_1e20_mid_m"),
        "radius_1e20_maxz_m": at_times("r_1e20_max_m"),
        "radius_sigma_mid_m": at_times("r_sigma_mid_m"),
        "length_1e20_m": at_times("length_1e20_m"),
        "radial_growth_speed_mid_1e20_m_s": vr_mid,
        "tip_velocity_mean_fit_m_s": tv.get("v_mean_m_s"),
        "tip_velocity_fit_window_s": tv.get("fit_window_s"),
        "tip_velocity_last_1ns_m_s": vlate,
        "ne_axis_z0_m3": [{"t_s": x["t_s"], "ne": x["ne_axis_mid"],
                           "ne_avg_r<=100um": x["ne_axis_mid_avg_r100um"]} for x in run["snaps"]],
        "peak_ne_any_m3": float(np.max(h["peak_ne_m3"])),
        "peak_ne_final_m3": float(h["peak_ne_m3"][-1]),
        "peak_E_V_m": float(np.max(h["peak_E_V_m"])),
        "delay_proxies": {
            "t_diam1e20_mid_ge_0p2mm_s": t_d_mid,
            "t_diam1e20_mid_ge_0p2mm_extrapolated_s": t_d_mid_extrap,
            "t_diam1e20_anyz_ge_0p2mm_s": t_d_max,
            "t_line_conductance_mid_x2_s": t_G2,
        },
    }


def compare(name_a, name_b, times=(1e-9, 2e-9, 3e-9)):
    """Compare run a against reference run b at common times (rel_diff = a/b - 1)."""
    a, b = _load_run(name_a), _load_run(name_b)
    if a is None or b is None:
        return None
    ha, hb = a["hist"], b["hist"]
    tmax = min(ha["t_s"][-1], hb["t_s"][-1])
    out = {"run": name_a, "reference": name_b, "common_window_s": float(tmax)}
    for t in [x for x in times if x <= tmax + 1e-15]:
        row = {}
        for k in ("r_1e20_mid_m", "r_1e20_max_m", "r_sigma_mid_m", "length_1e20_m", "peak_ne_m3"):
            va, vb = _at(ha, k, t), _at(hb, k, t)
            row[k] = {"run": va, "ref": vb, "rel_diff": (va / vb - 1) if vb else None}
        out[f"{t*1e9:g}ns"] = row
    # tip speed over a common window
    w0 = min(1e-9, 0.5 * tmax)
    va = [_rate_fit(h["t_s"][h["t_s"] <= tmax], h["z_plus_1e20_m"][h["t_s"] <= tmax], w0) for h in (ha, hb)]
    out["tip_speed_plus_fit_window_s"] = [w0, float(tmax)]
    out["tip_speed_plus_m_s"] = {"run": va[0], "ref": va[1],
                                 "rel_diff": (va[0] / va[1] - 1) if (va[0] and va[1]) else None}
    out["ne_axis_z0_last_common_snapshot"] = None
    sa = {round(x["t_s"] * 1e12): x["ne_axis_mid"] for x in a["snaps"]}
    sb = {round(x["t_s"] * 1e12): x["ne_axis_mid"] for x in b["snaps"]}
    common = sorted(set(sa) & set(sb))
    if common:
        k = common[-1]
        out["ne_axis_z0_last_common_snapshot"] = {"t_s": k * 1e-12, "run": sa[k], "ref": sb[k],
                                                  "rel_diff": sa[k] / sb[k] - 1}
    return out


def domain_check(case):
    """Compare +/-4 mm and +/-8 mm runs over their common time window."""
    a, b = _load_run(case), _load_run(case + "_z8")
    if a is None or b is None:
        return None
    ha, hb = a["hist"], b["hist"]
    tmax = min(ha["t_s"][-1], hb["t_s"][-1])
    out = {"common_window_s": float(tmax)}
    for k in ("r_1e20_mid_m", "r_1e20_max_m", "r_sigma_mid_m", "length_1e20_m", "peak_ne_m3"):
        va, vb = _at(ha, k, tmax), _at(hb, k, tmax)
        out[k] = {"z4mm": va, "z8mm": vb, "rel_diff": (va - vb) / vb if vb else None}
    return out


def pick_runs():
    runs = {}
    for c in CASES:
        for suf in PREFER:
            r = _load_run(c + suf)
            if r is not None:
                runs[c] = r
                break
    return runs


def build(write=True):
    runs = pick_runs()
    metrics = {c: case_metrics(r) for c, r in runs.items()}
    # delay-proxy difference low vs high field at 0.2 MPa
    dd = None
    lo, hi = metrics.get("He0p2MPa_E16p6"), metrics.get("He0p2MPa_E22p3")
    if lo and hi:
        dd = {}
        for k in ("t_diam1e20_mid_ge_0p2mm_s", "t_line_conductance_mid_x2_s", "t_diam1e20_anyz_ge_0p2mm_s"):
            a, b = lo["delay_proxies"][k], hi["delay_proxies"][k]
            if a is not None and b is not None:
                dd[k] = {"value_s": a - b}
            elif b is not None:
                dd[k] = {"lower_bound_s": lo["t_final_s"] - b, "note": "low-field case had not crossed by t_final"}
            else:
                dd[k] = None
        a = lo["delay_proxies"]["t_diam1e20_mid_ge_0p2mm_s"] or lo["delay_proxies"]["t_diam1e20_mid_ge_0p2mm_extrapolated_s"]
        b = hi["delay_proxies"]["t_diam1e20_mid_ge_0p2mm_s"] or hi["delay_proxies"]["t_diam1e20_mid_ge_0p2mm_extrapolated_s"]
        if a is not None and b is not None:
            dd["t_diam1e20_mid_ge_0p2mm_incl_linear_extrapolation_s"] = {"value_s": a - b}
    out = {
        "description": "streamer2d: 2D axisymmetric drift-diffusion-Poisson side-model, laser-seeded He filament in a "
                       "uniform 2.95 GHz field parallel to the channel. Simulation results (not validated against He "
                       "channel-growth data, which do not exist in the workspace).",
        "measured_reference": MEASURED,
        "delay_proxy_definition": "t at which the n_e>=1e20 m^-3 diameter at z=0 reaches 0.2 mm (lower end of the "
                                  "measured air filament diameter); secondary: time for the z=0 line conductance "
                                  "int(sigma 2 pi r dr) to double. Both are measured from the laser-formed initial state.",
        "delta_delay_proxy_16p6_minus_22p3": dd,
        "cases": metrics,
        "domain_size_check": {c: domain_check(c) for c in CASES if domain_check(c) is not None},
        "sensitivity": {
            "recombination_alpha_1e-15_vs_1e-14": compare("He0p2MPa_E22p3_alpha1e-15", "He0p2MPa_E22p3"),
            "grid_dz10um_drmin2um_vs_dz5um_drmin1um": compare("He0p2MPa_E22p3_coarse", "He0p2MPa_E22p3"),
            "grid_1MPa_dz10_vs_dz5": compare("He1MPa_46Td_coarse", "He1MPa_46Td", times=(0.5e-9, 1e-9, 1.5e-9)),
            "grid_1p825MPa_dz10_vs_dz5": compare("He1p825MPa_46Td_coarse", "He1p825MPa_46Td",
                                                 times=(0.5e-9, 0.96e-9)),
        },
        "background_preionization_note": (
            "Pure He has no photoionization. A 1e12 m^-3 background (vs 1e9) was not simulated: with the cycle-"
            "averaged nu_iz <= 8e8 s^-1 at 0.2 MPa / 46 Td it grows by at most exp(8) ~ 3e3 in 10 ns, to <~3e15 "
            "m^-3, 7 orders below the 1e20 contour, so it cannot change channel growth on this time scale."),
        "stats": {"total_wall_s_all_runs": sum(
            json.loads((d / "summary.json").read_text())["status"]["wall_s"]
            for d in RESULTS.iterdir() if (d / "summary.json").exists() and d.name != "figures")},
    }
    if write:
        (RESULTS / "summary.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out, runs


# ---------------------------------------------------------------- figures
def figures(runs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    figdir = RESULTS / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    made = []
    for c, run in runs.items():
        d = RESULTS / run["name"]
        files = [d / x["file"] for x in run["snaps"]]  # time-ordered
        times = np.array([x["t_s"] for x in run["snaps"]])
        tf = times[-1]
        want = [0.0] + [w for w in (2e-9, 5e-9, 10e-9) if w < tf - 0.2e-9] + [tf]
        if len(want) < 4:
            want = sorted(set(want) | {tf / 3, 2 * tf / 3})
        pick = []
        for w in want:
            k = int(np.argmin(np.abs(times - w)))
            if k not in pick:
                pick.append(k)
        pick = sorted(pick)[:5]
        zh = run["summary"]["derived"]["grid"]["z_half_m"]
        r_win = min(max(0.5e-3, 1.15 * float(np.nanmax(run["hist"]["r_1e20_max_m"]))), 1.5e-3)
        fig, axs = plt.subplots(len(pick), 1, figsize=(10, 1.9 * len(pick) + 0.6), sharex=True,
                                constrained_layout=True)
        axs = np.atleast_1d(axs)
        for ax, k in zip(axs, pick):
            z = np.load(files[k])
            rc, zc, ne = z["rc"], z["zc"], z["ne"].astype(float)
            m = rc <= r_win
            R = np.concatenate([-rc[m][::-1], rc[m]])
            Z = np.concatenate([ne[m][::-1], ne[m]], axis=0)
            pc = ax.pcolormesh(zc * 1e3, R * 1e3, np.maximum(Z, 1e14), norm=LogNorm(1e17, 1e24),
                               cmap="inferno", shading="auto", rasterized=True)
            ax.contour(zc * 1e3, R * 1e3, Z, levels=[1e20], colors="c", linewidths=0.7)
            ax.set_ylabel("r (mm)")
            ax.set_title(f"t = {float(z['t_s'])*1e9:.2f} ns   (E_app = {float(z['E_app'])/1e5:+.1f} kV/cm)",
                         fontsize=9)
        axs[-1].set_xlabel("z (mm)  [domain ±%.0f mm; grounded ends]" % (zh * 1e3))
        fig.colorbar(pc, ax=axs, label="n_e (m$^{-3}$), cyan = 1e20 contour", shrink=0.85)
        mtr = run["summary"]
        fig.suptitle(f"{run['name']}: He {mtr['config']['pressure_Pa']/1e6:g} MPa, E0 = "
                     f"{mtr['config']['E0_V_m']/1e5:.1f} kV/cm peak ({mtr['derived']['E0_over_N_Td']:.1f} Td). "
                     "Simulation; r axis mirrored", fontsize=9)
        p = figdir / f"ne_map_{c}.png"
        fig.savefig(p, dpi=120); plt.close(fig); made.append(p)

    # radius / length / n_e curves
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.5), constrained_layout=True)
    axs = axs.ravel()
    for i, (c, run) in enumerate(runs.items()):
        h = run["hist"]
        t = h["t_s"] * 1e9
        col = f"C{i}"
        axs[0].plot(t, h["r_1e20_mid_m"] * 1e6, color=col, lw=1.3, label=run["name"])
        axs[0].plot(t, h["r_1e20_max_m"] * 1e6, color=col, lw=0.8, ls="--")
        axs[1].plot(t, h["r_sigma_mid_m"] * 1e6, color=col, lw=1.2)
        axs[2].plot(t, h["length_1e20_m"] * 1e3, color=col, lw=1.2)
        sn = run["snaps"]
        axs[3].semilogy([x["t_s"] * 1e9 for x in sn], [x["ne_axis_mid"] for x in sn], "o-", color=col, lw=1, ms=3)
    axs[0].axhspan(100, 300, color="0.85", zorder=0, label="measured AIR radius 0.1-0.3 mm")
    axs[0].set(xlabel="t (ns)", ylabel="radius (µm)",
               title="n_e ≥ 1e20 m$^{-3}$ half-width: solid at z=0, dashed max over z")
    axs[0].legend(fontsize=7)
    axs[1].set(xlabel="t (ns)", ylabel="σ-weighted RMS radius at z=0 (µm)", title="conductivity-weighted radius")
    axs[2].set(xlabel="t (ns)", ylabel="axial extent of n_e ≥ 1e20 (mm)",
               title="channel length (curves end when the tip nears the domain end)")
    axs[3].axhline(6e22, color="k", ls="--", lw=1, label="He measured ~6e22 (Stark, line-of-sight)")
    axs[3].set(xlabel="t (ns)", ylabel="n_e(r=0, z=0) (m$^{-3}$)", title="on-axis electron density at mid-plane")
    axs[3].legend(fontsize=7)
    p = figdir / "radius_length_vs_time.png"
    fig.savefig(p, dpi=120); plt.close(fig); made.append(p)
    return made


def figures_summary(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figdir = RESULTS / "figures"
    made = []
    fig, ax = plt.subplots(figsize=(7, 4.8), constrained_layout=True)
    mk = {0.2: "o", 0.446: "s", 1.0: "D", 1.825: "^"}
    for i, (c, m) in enumerate(out["cases"].items()):
        v, vl = m["tip_velocity_mean_fit_m_s"], m["tip_velocity_last_1ns_m_s"]
        if v is None or not np.isfinite(v):
            continue
        p = round(m["pressure_MPa"], 3)
        ax.plot(m["E0_over_N_Td"], v, mk.get(p, "o"), color=f"C{i}", ms=7,
                label=f"{c}: fit t≥1 ns" + (f", last ns {vl:.1e}" if vl else ""))
        if vl:
            ax.plot([m["E0_over_N_Td"]] * 2, [v, vl], color=f"C{i}", lw=1)
            ax.plot(m["E0_over_N_Td"], vl, "_", color=f"C{i}", ms=12, mew=2)
    ax.axhline(5e5, color="k", ls="--", label="measured AIR filament growth ~5e5 m/s (not He)")
    ax.set_yscale("log")
    ax.set(xlabel="applied peak E0/N (Td)", ylabel="tip speed of n_e ≥ 1e20 contour (m/s)",
           title="Model tip speed vs E/N (bar: fit mean → last-ns speed)")
    ax.legend(fontsize=6.5)
    p = figdir / "tip_velocity_vs_EN.png"
    fig.savefig(p, dpi=120); plt.close(fig); made.append(p)

    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    names = list(out["cases"])
    for k, c in enumerate(names):
        m = out["cases"][c]
        for j, (key, colr, lab) in enumerate((("radius_1e20_mid_m", "C0", "n_e≥1e20 at z=0"),
                                              ("radius_1e20_maxz_m", "C2", "n_e≥1e20 max over z"),
                                              ("radius_sigma_mid_m", "C1", "σ-weighted at z=0"))):
            vals = [(tk, v) for tk, v in m[key].items() if v is not None and tk.endswith("ns") and tk != "0ns"]
            for tk, v in vals:
                ax.plot(k + (j - 1) * 0.22, 2 * v * 1e3, "o", color=colr, ms=3 + 1.2 * float(tk[:-2]) ** 0.5,
                        alpha=0.8, label=lab if (k == 0 and tk == vals[0][0]) else None)
    ax.axhspan(0.2, 0.6, color="0.85", zorder=0, label="measured AIR filament diameter 0.2-0.6 mm (not He)")
    ax.set_xticks(range(len(names)), names, rotation=25, ha="right", fontsize=7)
    ax.set_yscale("log")
    ax.set_ylabel("channel diameter (mm)")
    ax.set_title("Model diameter at 1, 2, 5, 10 ns (marker size grows with t) vs measured air value")
    ax.legend(fontsize=7)
    p = figdir / "diameter_vs_measured.png"
    fig.savefig(p, dpi=120); plt.close(fig); made.append(p)
    return made


def main():
    out, runs = build()
    made = figures(runs) + figures_summary(out)
    print(RESULTS / "summary.json")
    for p in made:
        print(p)
    hdr = f"{'case':28s} {'run':30s} {'t_f ns':>6s} {'r20mid 2/5/10 ns (um)':>24s} {'rσ 2/5/10':>18s} " \
          f"{'v_tip':>9s} {'v_last':>9s} {'ne_ax_final':>11s} {'t_d0.2mm':>9s}"
    print(hdr)
    f = lambda v, s=1e6: "  -  " if v is None else f"{v*s:5.0f}"
    for c, m in out["cases"].items():
        r, rs = m["radius_1e20_mid_m"], m["radius_sigma_mid_m"]
        dp = m["delay_proxies"]
        td = dp["t_diam1e20_mid_ge_0p2mm_s"]
        tde = dp["t_diam1e20_mid_ge_0p2mm_extrapolated_s"]
        print(f"{c:28s} {m['run']:30s} {m['t_final_s']*1e9:6.2f} {f(r['2ns'])}/{f(r['5ns'])}/{f(r['10ns'])}"
              f"          {f(rs['2ns'])}/{f(rs['5ns'])}/{f(rs['10ns'])} "
              f"{(m['tip_velocity_mean_fit_m_s'] or float('nan')):9.2e} {(m['tip_velocity_last_1ns_m_s'] or float('nan')):9.2e} "
              f"{m['ne_axis_z0_m3'][-1]['ne']:11.2e} "
              + (f"{td*1e9:7.2f}ns" if td else (f"~{tde*1e9:.1f}ns(ex)" if tde else "   -")))
    print("delta delay proxy:", json.dumps(out["delta_delay_proxy_16p6_minus_22p3"]))
    print("domain check:", json.dumps(out["domain_size_check"], indent=0))


if __name__ == "__main__":
    main()
