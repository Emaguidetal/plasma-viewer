"""Figures and the interactive 3D viewer for the plasma-switch study in a WR-284 cavity.

* rod growth along E (rod3d): results/rod3d/rod3d_overview.png and results/rod3d/<case>/rod3d_3d.html
* continued RF drive of a bridged rod (rod_drive): results/rod_drive/rod_drive.png
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
R3 = ROOT / "results" / "rod3d"
RD = ROOT / "results" / "rod_drive"
A_WG, B_WG, D_CAV = 72.136e-3, 34.036e-3, 69.27e-3
PLOTLY = "https://cdn.plot.ly/plotly-2.35.2.min.js"


def _load3(name):
    d = R3 / name
    if not (d / "rod3d.npz").exists():
        return None, None
    return np.load(d / "rod3d.npz"), json.loads((d / "summary.json").read_text(encoding="utf-8"))


def rod3d_overview(field_cases=("E50", "E75", "E100", "E150", "E200", "E300"),
                   var_cases=("E200", "E200_s2", "E200_s3", "E200_ph90", "E200_h025"),
                   laser_cases=("E100", "E100_0p5mJ", "E100_1mJ", "E100_f10", "E100_f20")):
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    cmap = plt.get_cmap("viridis")
    rows = []
    for k, name in enumerate(field_cases):
        h, s = _load3(name)
        if h is None:
            continue
        col = cmap(k / max(len(field_cases) - 1, 1))
        t = 1e9 * h["hist_t_s"]
        ax[0, 0].plot(t, 1e3 * h["hist_max_len_m"], color=col, label=name)
        ax[0, 1].plot(t, h["hist_e0_v_m"] / 1e5, color=col, label=name)
        rows.append(s)
    ax[0, 0].axhline(1e3 * B_WG, color="k", ls="--", lw=0.8)
    ax[0, 0].text(1, 1e3 * B_WG + 0.5, "broad-wall gap b = 34 mm", fontsize=8)
    ax[0, 0].set(xlabel="time after laser (ns)", ylabel="longest rod (mm)", title="(a) rod growth along E")
    ax[0, 1].axvline(100, color="grey", ls=":", lw=0.8)
    ax[0, 1].set(xlabel="time after laser (ns)", ylabel="cavity field E0 (kV/cm)",
                 title="(b) cavity field (drive off at 100 ns)")
    for a in ax[0, :2]:
        a.legend(fontsize=8)
    # (c) rod pattern along the laser line at end of run
    for k, s in enumerate(rows):
        for r in s["rods"]:
            ax[0, 2].plot([r["x_mm"]] * 2, [k - 0.4 * r["length_mm"] / 34.0, k + 0.4 * r["length_mm"] / 34.0],
                          color="tab:red" if r["bridged_ns"] is not None else "tab:blue", lw=3)
    ax[0, 2].set_yticks(range(len(rows)), [s["name"] for s in rows])
    ax[0, 2].set(xlabel="x along the laser line (mm, 0 = field maximum)",
                 title="(c) rods at end (bar ~ length; red = bridged)")
    # (d) shot-to-shot / numerics at 200 kV/cm, (e) laser energy and focus at 100 kV/cm
    for a, cases, ttl in ((ax[1, 0], var_cases, "(d) 200 kV/cm: RNG seed, RF phase, element h"),
                          (ax[1, 1], laser_cases, "(e) 100 kV/cm: laser energy and focus offset")):
        for name in cases:
            h, s = _load3(name)
            if h is None:
                continue
            a.plot(1e9 * h["hist_t_s"], 1e3 * h["hist_max_len_m"], label=name)
        a.axhline(1e3 * B_WG, color="k", ls="--", lw=0.8)
        a.set(xlabel="time after laser (ns)", ylabel="longest rod (mm)", title=ttl)
        a.legend(fontsize=8)
    # (f) summary table
    ax[1, 2].axis("off")
    lines = ["case        birth(ns) rods long bridge(ns) E0min*",
             "------------------------------------------------"]
    for name in dict.fromkeys(field_cases + var_cases + laser_cases):
        h, s = _load3(name)
        if s is None:
            continue
        b = s["t_first_birth_ns"]
        tb = s["t_bridge_ns"]
        on = h["hist_t_s"] < s["config"]["drive_off_s"]
        e0min = float(h["hist_e0_v_m"][on].min()) / 1e5 if on.any() else float("nan")
        lines.append(f"{name:<11} {('%.1f' % b) if b else '  -':>8} {s['n_rods']:>4} {s['n_long_rods']:>4} "
                     f"{('%.1f' % tb) if tb else '   -':>9} {e0min:6.0f}")
    lines.append("* min cavity field (kV/cm) while driven")
    ax[1, 2].text(0, 1, "\n".join(lines), family="monospace", fontsize=8, va="top")
    ax[1, 2].set_title("(f) summary (long = length > b/4)")
    fig.suptitle("3D rod growth along E in a WR-284 cavity, 1.48 MPa He, laser perpendicular to E "
                 "(tree model calibrated on resolved 2D runs)")
    fig.tight_layout()
    out = R3 / "rod3d_overview.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def rod3d_html(name, max_frames=160):
    """Self-contained Plotly animation of the channel tree (nodes coloured by log10 n_e) in the cavity."""
    h, s = _load3(name)
    if h is None:
        return None
    pos, rad, mf, lam, t = h["pos"], h["rad"], h["m_frame"], h["lam"], h["t_s"]
    idx = np.unique(np.linspace(0, len(t) - 1, min(max_frames, len(t))).round().astype(int))
    frames = []
    for k in idx:
        m = int(mf[k])
        ne = lam[k, :m] / (math.pi * np.maximum(rad[:m], 1e-9) ** 2)
        keep = ne > 1e18
        p = pos[:m][keep]
        frames.append(dict(t=round(1e9 * float(t[k]), 3), e0=round(float(h["hist_e0_v_m"][min(k, len(h["hist_e0_v_m"]) - 1)]) / 1e5, 1),
                           x=np.round(1e3 * (p[:, 0] - A_WG / 2), 3).tolist(), y=np.round(1e3 * p[:, 1], 3).tolist(),
                           z=np.round(1e3 * (p[:, 2] - D_CAV / 2), 3).tolist(),
                           c=np.round(np.log10(ne[keep]), 2).tolist()))
    a, b = 12.0, 1e3 * B_WG                         # zoom window: |x|, |z| < 12 mm around the laser crossing
    box = dict(x=[-a, a, a, -a, -a, None, -a, a, a, -a, -a], y=[0, 0, 0, 0, 0, None, b, b, b, b, b],
               z=[-a, -a, a, a, -a, None, -a, -a, a, a, -a])
    on = h["hist_t_s"] < s["config"]["drive_off_s"]
    summ = {k: s[k] for k in ("t_first_birth_ns", "n_rods", "n_long_rods", "t_bridge_ns")}
    summ["e0_min_kv_cm"] = float(h["hist_e0_v_m"][on].min()) / 1e5
    data = dict(name=name, frames=frames, box=box, b=b, a=a, summary=summ, e0=s["config"]["e0_v_m"] / 1e5)
    html = f"""<!doctype html><html><head><meta charset="utf-8"><title>3D rod growth {name}</title>
<script src="{PLOTLY}"></script>
<style>body{{background:#111;color:#ddd;font-family:Segoe UI,Arial;margin:12px}} #p{{width:100%;height:78vh}}
a{{color:#7cc4ff}} .s{{font-size:13px;color:#aaa}}</style></head><body>
<h3>3D rod growth along E: WR-284 cavity, 1.48 MPa He, E0 = {data['e0']:.0f} kV/cm, laser along x (perpendicular to E)</h3>
<div class="s">Points are channel elements coloured by log10 n<sub>e</sub> (m<sup>-3</sup>); the frame is the cavity
(zoomed to |x|, |z| &lt; 12 mm of the 72 x 69 mm cavity; y is along E between the broad walls, b = 34 mm). The laser line is along x at y = b/2.
Tree model calibrated on resolved axisymmetric runs; see the report for assumptions.
<span id="info"></span></div>
<input type="range" id="sl" min="0" max="{len(frames) - 1}" value="0" style="width:60%">
<button id="pl">play</button> <span id="tt"></span>
<div id="p"></div>
<script>
const D={json.dumps(data, separators=(',', ':'))};
const box={{type:'scatter3d',mode:'lines',x:D.box.x,y:D.box.z,z:D.box.y,line:{{color:'#888',width:2}},name:'broad walls (y = 0, b)',hoverinfo:'skip'}};
const laser={{type:'scatter3d',mode:'lines',x:[-D.a,D.a],y:[0,0],z:[D.b/2,D.b/2],line:{{color:'#2f2',width:3,dash:'dash'}},name:'laser line'}};
function tr(f){{return {{type:'scatter3d',mode:'markers',x:f.x,y:f.z,z:f.y,name:'channel',
 marker:{{size:3,color:f.c,cmin:18,cmax:24,colorscale:'Hot',colorbar:{{title:'log10 ne',len:0.6}}}}}};}}
const lay={{paper_bgcolor:'#111',font:{{color:'#ddd'}},margin:{{l:0,r:0,t:0,b:0}},
 scene:{{xaxis:{{title:'x (mm)',range:[-D.a,D.a]}},yaxis:{{title:'z (mm)',range:[-D.a,D.a]}},zaxis:{{title:'y along E (mm)',range:[0,D.b]}},
 aspectmode:'manual',aspectratio:{{x:0.8,y:0.8,z:1.13}},camera:{{eye:{{x:1.4,y:1.4,z:0.7}}}}}}}};
Plotly.newPlot('p',[box,laser,tr(D.frames[0])],lay);
const s=D.summary; document.getElementById('info').innerHTML=' First rod '+(s.t_first_birth_ns?s.t_first_birth_ns.toFixed(1)+' ns':'none')+
 '; rods '+s.n_rods+' (long '+s.n_long_rods+'); bridge '+(s.t_bridge_ns?s.t_bridge_ns.toFixed(1)+' ns':'not reached')+'; min E0 while driven '+s.e0_min_kv_cm.toFixed(0)+' kV/cm.';
function show(k){{const f=D.frames[k];Plotly.restyle('p',{{x:[f.x],y:[f.z],z:[f.y],'marker.color':[f.c]}},[2]);
 document.getElementById('tt').textContent='t = '+f.t.toFixed(2)+' ns, E0 = '+f.e0+' kV/cm';}}
const sl=document.getElementById('sl'); sl.oninput=()=>show(+sl.value); show(0);
let timer=null; document.getElementById('pl').onclick=()=>{{if(timer){{clearInterval(timer);timer=null;return;}}
 timer=setInterval(()=>{{sl.value=(+sl.value+1)%D.frames.length;show(+sl.value);}},120);}};
</script></body></html>"""
    out = R3 / name / "rod3d_3d.html"
    out.write_text(html, encoding="utf-8")
    return out


def rod_drive_figure(cases=("chain_E200", "chain_E300", "chain_E200_cw", "chain_E300_cw", "chain_E200_cw_r100",
                            "chain_E200_cw_lowP", "chain_E200_cw_dt2"), ref="chain_E200_cw"):
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    for name in cases:
        f = RD / name / "rod_drive.npz"
        if not f.exists():
            continue
        d = np.load(f)
        sm = json.loads((RD / name / "summary.json").read_text(encoding="utf-8"))
        t = 1e9 * d["t_s"]
        ls = "--" if name.endswith("dt2") else "-"
        ax[0, 0].semilogx(t, d["e_gap_v_m"] / 1e5, ls, label=name)
        c = ax[0, 0].lines[-1].get_color()
        ax[0, 1].loglog(t, d["r_ohm"], ls, color=c, label=name)
        ax[0, 2].loglog(t, np.maximum(d["p_rod_w"], 1), ls, color=c, label=name)
        ax[1, 0].loglog(t, d["ne_axis_m3"], ls, color=c, label=name)
        ax[1, 1].semilogx(t, d["tg_axis_k"], ls, color=c, label=name)
        for key, mk in (("t_contraction_s", "o"), ("t_skin_invalid_s", "x")):
            if sm.get(key):
                i = int(np.searchsorted(d["t_s"], sm[key]))
                ax[1, 1].plot(t[i], d["tg_axis_k"][i], mk, color=c, ms=7)
    ax[1, 1].plot([], [], "ko", label="thermal contraction (arc onset)")
    ax[1, 1].plot([], [], "kx", label="skin depth < 2 r_c (validity ends)")
    d = np.load(RD / ref / "rod_drive.npz")
    ax[0, 1].loglog(1e9 * d["t_s"], d["x_ohm"], "k:", label=f"post reactance X ({ref})")
    ax[0, 2].loglog(1e9 * d["t_s"], np.maximum(d["p_refl_w"], 1), "k:", label=f"reflected ({ref})")
    ax[0, 0].set(ylabel="kV/cm", title="(a) TE10 modal field |V|/b at the post")
    ax[0, 1].set(ylabel="ohm", title="(b) rod resistance R and post reactance X")
    ax[0, 2].set(ylabel="W", title="(c) power into the rod / reflected (nominal power)")
    ax[1, 0].set(ylabel="m^-3", title="(d) axis electron density")
    ax[1, 1].set(ylabel="K", title="(e) axis gas temperature")
    for a in ax.flat[:5]:
        a.set_xlabel("time after bridging (ns)")
        a.set_xlim(0.1, 1000)
        a.legend(fontsize=6)
    r = 1e3 * d["r_m"]
    for k in np.unique(np.searchsorted(d["t_s"], np.array([0, 3, 30, 150, 240, 260, 400, 1000]) * 1e-9)
                       .clip(0, len(d["t_s"]) - 1)):
        ax[1, 2].semilogy(r, np.maximum(d["prof_ne"][k], 1e15), label=f"{1e9 * d['t_s'][k]:.0f} ns")
    a2 = ax[1, 2].twinx()
    for k in np.unique(np.searchsorted(d["t_s"], np.array([0, 3, 30, 150, 240, 260, 400, 1000]) * 1e-9)
                       .clip(0, len(d["t_s"]) - 1)):
        a2.plot(r, d["prof_rho"][k], ls=":", lw=0.8)
    a2.set_ylabel("rho / rho0 (dotted)")
    ax[1, 2].set(xlabel="r (mm)", ylabel="ne (m^-3, solid)", xlim=(0, 1.5), ylim=(1e18, 1e25),
                 title=f"(f) {ref} radial profiles")
    ax[1, 2].legend(fontsize=7)
    fig.suptitle("Bridged rod under RF drive, chained from rod3d (WR-284, 1.48 MPa He). Nominal, uncalibrated power; "
                 "values after the x markers are outside model validity. Simulated, not measured.", fontsize=10)
    fig.tight_layout()
    out = RD / "rod_drive.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main():
    print(rod3d_overview())
    for p in sorted(R3.iterdir()):
        if (p / "rod3d.npz").exists():
            print(rod3d_html(p.name))
    print(rod_drive_figure())


if __name__ == "__main__":
    main()
