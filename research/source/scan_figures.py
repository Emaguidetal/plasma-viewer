"""Figures and experiment-comparison table for the p x E, laser, seed and anchor scans.

usage: python results/scan_figures.py      (after scan_analysis.py scan_pE scan_laser anchor conv)
Writes results/coupled3d/figures_scan/*.png and exp_compare.md / scan_summary.json.
"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))
sys.path.insert(0, str(HERE))
from rf_plasma.stochastic_seed import Rates  # noqa: E402
import scan_analysis  # noqa: E402

KB = 1.380649e-23
BASE = HERE / "coupled3d"
OUT = BASE / "figures_scan"
EXP = HERE.parents[2] / "literature" / "experimental_data.json"
F_RF = 2.95e9


def rows_for(root):
    out = []
    for d in sorted((BASE / root).glob("*")):
        if (d / "summary.json").exists():
            r = scan_analysis.analyse(d)
            r["dir"] = d
            out.append(r)
    return out


def surrogate_delay(rates, e_kvcm, p_pa, n_br, n0):
    n = p_pa / (KB * 300.0)
    nu = float(rates.nu_iz(e_kvcm * 1e5, n, F_RF))
    return math.log(max(n_br / n0, 1.0)) / nu * 1e9 if nu > 0 else np.inf


def fig_pE(pe, rates):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ps = sorted({round(r["p_psig"]) for r in pe})
    es = sorted({round(r["E_kVcm"]) for r in pe})
    g = np.full((len(ps), len(es)), np.nan)
    lab = {}
    for r in pe:
        i, j = ps.index(round(r["p_psig"])), es.index(round(r["E_kVcm"]))
        g[i, j] = r["gamma_end"]
        lab[i, j] = "no BD" if r["delay03_ns"] is None else f"{r['delay03_ns']:.1f} ns"
    im = ax[0].imshow(g, origin="lower", cmap="viridis", vmin=0, vmax=0.5, aspect="auto")
    for (i, j), s in lab.items():
        ax[0].text(j, i, s, ha="center", va="center", color="w" if g[i, j] < 0.3 else "k", fontsize=8)
    ax[0].set_xticks(range(len(es)), [str(e) for e in es])
    ax[0].set_yticks(range(len(ps)), [str(p) for p in ps])
    ax[0].set_xlabel("RF field at focus E0 (kV/cm, peak)")
    ax[0].set_ylabel("He pressure (psig)")
    ax[0].set_title("|Γ| at end (colour) and delay to |Γ|>0.3 after laser")
    fig.colorbar(im, ax=ax[0], label="|Γ| end")
    # delay vs E/N with the 0D surrogate (N_br from the 3D runs, N0 = laser seed electrons)
    nbr = np.median([r["n_electrons_at_br"] for r in pe if r["n_electrons_at_br"]])
    n0 = np.median([r["seed_electrons"] for r in pe if r["seed_electrons"]])
    cm = plt.get_cmap("plasma")
    for k, p in enumerate(ps):
        sel = sorted([r for r in pe if round(r["p_psig"]) == p and r["delay03_ns"]], key=lambda r: r["E_N_Td"])
        c = cm(k / max(len(ps) - 1, 1))
        if sel:
            ax[1].plot([r["E_N_Td"] for r in sel], [r["delay03_ns"] for r in sel], "o-", color=c,
                       label=f"3D, {p} psig")
        p_pa = p * 6894.757 + 101325.0
        ee = np.linspace(min(es), max(es), 60)
        td = ee * 1e5 / (p_pa / (KB * 300)) * 1e21
        ax[1].plot(td, [surrogate_delay(rates, e, p_pa, nbr, n0) for e in ee], ":", color=c, lw=1)
    ax[1].set_yscale("log")
    ax[1].set_ylim(0.5, 300)
    ax[1].set_xlabel("E0/N (Td)")
    ax[1].set_ylabel("delay after laser (ns)")
    ax[1].set_title(f"3D (solid) vs 0D avalanche ln(N_br/N0)/ν_iz (dotted)\nN_br={nbr:.1e}, N0={n0:.1e}")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "scan_pE.png", dpi=130)
    plt.close(fig)
    return dict(n_br_median=nbr, n_seed_median=n0)


def fig_laser(la, ref):
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    base = [r for r in la + ref if abs(r["E_mJ"] - 1) < 1e-6 and abs(r["laser_seed_k"] - 6) < 1e-6
            and round(r["p_psig"]) == 100 and round(r["E_kVcm"]) == 80 and r["root"] in ("scan_laser", "scan_pE")]
    for fw, m in ((0.5, "s"), (1.0, "o"), (2.0, "^")):
        sel = sorted({r["waist_um"]: r for r in base if abs(r["fwhm_ns"] - fw) < 1e-6}.values(),
                     key=lambda r: r["waist_um"])
        if not sel:
            continue
        w = [r["waist_um"] for r in sel]
        ax[0].plot(w, [r["delay03_ns"] or np.nan for r in sel], m + "-", label=f"{fw} ns")
        ax[1].plot(w, [r["r50_um_end"] for r in sel], m + "-", label=f"{fw} ns")
        ax[2].plot(w, [r["fil_tg_axis_max"] for r in sel], m + "-", label=f"{fw} ns")
    for r in la:
        if abs(r["E_mJ"] - 0.5) < 1e-6:
            ax[0].plot(r["waist_um"], r["delay03_ns"], "kx", ms=9)
            ax[0].annotate("0.5 mJ", (r["waist_um"], r["delay03_ns"]), fontsize=7, xytext=(4, -10),
                           textcoords="offset points")
        if abs(r["laser_seed_k"] - 3) < 1e-6:
            ax[0].plot(r["waist_um"], r["delay03_ns"], "rP", ms=9)
            ax[0].annotate("k=3", (r["waist_um"], r["delay03_ns"]), fontsize=7, xytext=(4, 4),
                           textcoords="offset points")
    for a, t in zip(ax, ("delay to |Γ|>0.3 (ns)", "channel r50 at end (µm)", "filament axis Tg max (K)")):
        a.set_xscale("log")
        a.set_xlabel("1/e² waist (µm)")
        a.set_ylabel(t)
        a.legend(title="FWHM", fontsize=8)
    ax[2].set_yscale("log")
    fig.suptitle("Laser focus / pulse-length scan: guide, He 100 psig, 80 kV/cm, 1 mJ (k=6 seed)")
    fig.tight_layout()
    fig.savefig(OUT / "scan_laser.png", dpi=130)
    plt.close(fig)


def fig_channel(pe, la):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    pick = [r for r in pe if round(r["p_psig"]) == 50] + \
        [r for r in la if round(r["fwhm_ns"] * 10) == 10 and abs(r["E_mJ"] - 1) < 1e-6
         and abs(r["laser_seed_k"] - 6) < 1e-6]
    pick = sorted(pick, key=lambda r: (r["root"], r["E_kVcm"], r["waist_um"]))
    for r in pick:
        ts = np.load(r["dir"] / "timeseries.npz")
        t = (ts["t"] - (5e-9 if r["laser"] else 0)) * 1e9
        lab = (f"50 psig, {r['E_kVcm']:.0f} kV/cm" if r["root"] == "scan_pE"
               else f"w={r['waist_um']:.0f} µm")
        ls = "-" if r["root"] == "scan_pE" else "--"
        ax[0].plot(t, ts["fil_r50_focus"] * 1e6, ls, label=lab)
        ax[1].plot(t, ts["fil_ne_axis_focus"], ls, label=lab)
    ax[0].axhspan(100, 300, color="0.85", label="air, measured 0.2–0.6 mm diam.")
    ax[0].set_ylabel("channel half-electron radius r50 (µm)")
    ax[1].set_yscale("log")
    ax[1].set_ylim(1e18, 1e26)
    ax[1].axhline(6e22, color="k", ls=":", label="He Stark 6e22 m⁻³ (different averaging)")
    ax[1].set_ylabel("axis n_e at focus (m⁻³)")
    for a in ax:
        a.set_xlabel("time after laser peak (ns)")
        a.set_xlim(left=-2)
        a.legend(fontsize=7)
    fig.suptitle("Channel growth in the 3D sub-cell filament (guide, 1 mJ / 1 ns)")
    fig.tight_layout()
    fig.savefig(OUT / "channel_growth.png", dpi=130)
    plt.close(fig)


def anchor_label(r):
    n = r["run"]
    for t in ("wide", "N1e5", "N1e8"):
        if n.endswith(t):
            return t
    return "nominal"


def fig_seed(conv, anc, rates, nominal):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    seeds = [r for r in conv if "seedx" in r["run"]] + ([nominal] if nominal else [])
    seeds = sorted([r for r in seeds if r["seed_electrons"]], key=lambda r: r["seed_electrons"])
    if seeds:
        ax[0].semilogx([r["seed_electrons"] for r in seeds], [r["delay03_ns"] or np.nan for r in seeds], "o-",
                       label="3D guide 100 psig, 80 kV/cm")
        r0 = seeds[0]
        n = r0["p_MPa"] * 1e6 / (KB * 300.0)
        nu = float(rates.nu_iz(r0["E_kVcm"] * 1e5, n, F_RF))
        x = np.array([seeds[0]["seed_electrons"], seeds[-1]["seed_electrons"]])
        ax[0].semilogx(x, r0["delay03_ns"] - np.log(x / x[0]) / nu * 1e9, "k:",
                       label=f"0D slope −1/ν_iz ({1e9 / nu:.2f} ns per e-fold)")
    ax[0].set_xlabel("laser seed electrons N0")
    ax[0].set_ylabel("delay to |Γ|>0.3 (ns)")
    ax[0].set_title("Seed-number sensitivity (3D)")
    ax[0].legend(fontsize=8)
    # anchors: delay difference 16.6 - 22.3 kV/cm
    groups = {}
    for r in anc:
        groups.setdefault(anchor_label(r), {})[round(r["E_kVcm"], 1)] = r
    res = {}
    xs = []
    for k, gr in groups.items():
        if 16.6 in gr and 22.3 in gr:
            lo, hi = gr[16.6], gr[22.3]
            d = {m: (lo[m] - hi[m]) if lo[m] is not None and hi[m] is not None else None
                 for m in ("t_Nhalf_ns", "t_Ecollapse_ns")}
            res[k] = dict(N0=hi["seed_electrons"], waist_um=hi["waist_um"], E_mJ=hi["E_mJ"],
                          t_Ecollapse_16p6=lo["t_Ecollapse_ns"], t_Ecollapse_22p3=hi["t_Ecollapse_ns"],
                          t_Nhalf_16p6=lo["t_Nhalf_ns"], t_Nhalf_22p3=hi["t_Nhalf_ns"],
                          dt_Ecollapse=d["t_Ecollapse_ns"], dt_Nhalf=d["t_Nhalf_ns"],
                          N_at_Ecollapse=(lo["n_electrons_at_Ecollapse"], hi["n_electrons_at_Ecollapse"]),
                          gamma_end=(lo["gamma_end"], hi["gamma_end"]))
            if d["t_Ecollapse_ns"] is not None:
                xs.append(hi["seed_electrons"])
                ax[1].plot(hi["seed_electrons"], d["t_Ecollapse_ns"], "o", ms=9,
                           label=f"3D {k} (w={hi['waist_um']:.0f} µm)")
    dd = (30.0, 15.0)  # Beilin 2015: delay(16.6) - delay(22.3 kV/cm), He 0.2 MPa
    ax[1].axhspan(dd[0] - dd[1], dd[0] + dd[1], color="0.85",
                  label="published light-onset difference: 30±15 ns")
    n0 = np.logspace(3, 12, 50)
    nbr = 2.5e13  # geometric mean of N at field collapse in the nominal cavity anchors (1.1e13, 5.7e13)
    p = 0.2e6
    sur = [surrogate_delay(rates, 16.6, p, nbr, x) - surrogate_delay(rates, 22.3, p, nbr, x) for x in n0]
    ax[1].semilogx(n0, sur, "k:", label=f"0D avalanche, N_br={nbr:.1e}")
    ax[1].set_xlabel("seed electrons N0")
    ax[1].set_ylabel("delay(16.6) − delay(22.3 kV/cm)  (ns)")
    ax[1].set_title("Field-collapse proxy vs measured light onset")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "seed_anchor.png", dpi=130)
    plt.close(fig)
    return res


def load_exp():
    try:
        return json.loads(EXP.read_text(encoding="utf-8"))
    except Exception:
        return None


KEEP = ("root", "run", "geometry", "p_psig", "p_MPa", "E_kVcm", "E_N_Td", "E_mJ", "fwhm_ns", "waist_um",
        "laser_seed_k", "seed_electrons", "gamma_end", "delay03_ns", "t_Nhalf_ns", "t_Ecollapse_ns",
        "fil_ne_axis_max", "fil_te_axis_max", "fil_tg_axis_max", "r50_um_end", "closure_fil", "dx_mm")


def _history(r, n=240):
    ts = np.load(r["dir"] / "timeseries.npz")
    rf = np.load(r["dir"] / "rf_series.npz")
    t_l = 5e-9 if r["laser"] else 0.0
    tt = (ts["t"] - t_l) * 1e9
    idx = np.unique(np.linspace(0, len(tt) - 1, n).astype(int))
    fwd = np.abs(rf["fwd"])
    g = np.abs(rf["bwd"]) / np.where(fwd > 0, fwd, np.inf)
    tg = (rf["t_phasor"] - t_l) * 1e9
    gi = np.unique(np.linspace(0, len(tg) - 1, n).astype(int))
    f = lambda a: [float(f"{x:.4g}") for x in a]  # noqa: E731
    return dict(t=f(tt[idx]), ne=f(ts["fil_ne_axis_focus"][idx]), r50=f(ts["fil_r50_focus"][idx] * 1e6),
                te=f(ts["fil_te_axis_focus"][idx]), N=f((ts["fil_electrons"] + ts["bulk_electrons"])[idx]),
                tg=f(tg[gi]), g=f(g[gi]))


def label(r):
    if r["root"] == "anchor":
        return f"cavity 0.2 MPa {r['E_kVcm']:.1f} kV/cm, {anchor_label(r)} (N0={r['seed_electrons']:.1e})"
    if r["root"] == "conv":
        return f"conv {r['run'].split('_')[-1]}"
    s = f"{r['p_psig']:.0f} psig {r['E_kVcm']:.0f} kV/cm"
    if r["root"] == "scan_laser":
        s += f", {r['E_mJ']:g} mJ {r['fwhm_ns']:g} ns w{r['waist_um']:.0f}" + (" k3" if r["laser_seed_k"] == 3 else "")
    return s


def exp_rows(s, pe, anc):
    """Model vs measurement (Beilin 2015, Yatom 2016, Shlapakovski 2016). Strings for md/html."""
    a = s.get("anchor", {})
    nom, wide = a.get("nominal", {}), a.get("wide", {})
    pk = max((r["fil_ne_axis_max"] for r in anc), default=None)
    end = [r for r in anc if anchor_label(r) == "nominal"]
    ne_end = ", ".join(f"{r['E_kVcm']:.1f}: {np.load(r['dir'] / 'timeseries.npz')['fil_ne_axis_focus'][-1]:.1e}"
                       for r in end)
    fmt = lambda v: "—" if v is None else f"{v:.1f}"  # noqa: E731
    ab = [(k, v) for k, v in a.items() if v.get("dt_Ecollapse") is not None]
    dd = "; ".join(f"{k} (N0 {v['N0']:.0e}): {fmt(v['dt_Ecollapse'])}" for k, v in ab)
    return [
        ("Delay difference, low vs high power (He 0.2 MPa)", "Laser-to-light onset: +30 ± 15 ns",
         dd or "pending", "Field-collapse proxy is not light onset; ±15 ns is the sum of reported jitters; seed is assumed"),
        ("Delay jitter", "±6 ns (22.3), ±9 ns (16.6)", "σ_P ≈ 5–10 % magnetron power (0D Monte Carlo)",
         "deterministic 3D runs give no jitter; see stochastic_seed"),
        ("Channel n_e (He)", "6×10²² m⁻³ (Stark)", f"peak axis {pk:.1e}; end ({ne_end})" if pk else "—",
         "model is axis value of a 20–40 µm core; Stark has different spatial and temporal averaging"),
        ("Channel diameter", "0.2–0.6 mm (air)", "3D r50 45–100 µm; 2D streamer 0.15–0.21 mm at 10 ns",
         "air measurement and He simulation are not matched"),
        ("Ionization-front speed", "5×10⁵ m/s (air)", "2D: 4.5–7×10⁵ m/s (22.3), ~1×10⁵ (16.6)",
         "air measurement and He simulation are not matched"),
        ("Breakdown threshold (guide, 15 ns)", "—", "E/N ≈ 20–30 Td (see p×E map)",
         "published operating fields are not measured breakdown thresholds"),
        ("Focused laser along E (Shlapakovski 2016)",
         "Air, center: peak-power ratio 0.47; back wall: output delay < 1 ns, ratio ≈ 0.35",
         "He proxy: 1–4 ns at ≥ 75 Td, 20 µm focus",
         "different gases and focus positions; peak-power ratio is not integrated energy efficiency"),
    ]


def write_md(rows):
    h = ("The He anchors match pressure and nominal field values only. The project cavity is at 2.95 GHz; "
         "the published compressor is about 2.767 GHz, with different geometry, seed and detector. "
         "Published field values are calculated from power; peak versus RMS is unclear. "
         "The measured +30 ± 15 ns is a laser-to-light-onset difference between output-power groups; "
         "±15 ns is the sum of the reported jitters. Field-collapse and light-onset delays are different observables.\n\n"
         "| Quantity | Measured | Model | Comment |\n|---|---|---|---|\n")
    (OUT / "exp_compare.md").write_text(h + "".join(f"| {a} | {b} | {c} | {d} |\n" for a, b, c, d in rows),
                                        encoding="utf-8")


PAGE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>He laser-triggered RF breakdown: parameter scans</title>
<style>body{font-family:system-ui,Segoe UI,sans-serif;background:#0c1622;color:#e3edf6;margin:0;padding:14px}
h1{font-size:20px;margin:4px 0}h2{font-size:15px;margin:4px 0 8px}.card{background:#13222f;border:1px solid #294158;border-radius:8px;padding:10px;margin:10px 0}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,560px),1fr));gap:10px}canvas{width:100%;height:auto;background:#0f1b27;border-radius:6px}
img{width:100%;border-radius:6px;background:#fff}table{border-collapse:collapse;font-size:12px;width:100%}td,th{border-bottom:1px solid #294158;padding:3px 5px;text-align:left}
.small{font-size:12px;color:#a9bccd}select,button{padding:4px 6px}a{color:#7cc4ff}.tw{overflow-x:auto}</style></head><body>
<h1>Laser-triggered RF breakdown in helium — parameter scans (coupled 3D FDTD + kinetics + sub-cell filament)</h1>
<div class="small">Guide: matched WR-284, 2.95 GHz, He, 1 mm mesh, laser 532 nm, localized k-order seed.
The cavity anchors use the published 0.2 MPa and nominal 16.6 / 22.3 kV/cm values, but differ in geometry, frequency,
seed and measured observable. The published fields are calculated from power; peak versus RMS is unclear.
The field-collapse proxy is not the measured laser-to-light-onset delay and is not independent validation.
<a href="index.html">← main 3D viewer</a></div>
<div class="grid">
<section class="card"><h2>p × E map</h2><label>Metric <select id="met">
<option value="gamma_end">|Γ| at end</option><option value="delay03_ns">delay to |Γ|>0.3 (ns)</option><option value="t_Ecollapse_ns">focus-field collapse (ns)</option>
<option value="r50_um_end">channel r50 end (µm)</option><option value="fil_ne_axis_max">axis n_e max (m⁻³)</option><option value="fil_tg_axis_max">axis T_gas max (K)</option></select></label>
<canvas id="map" width="640" height="380"></canvas><div class="small">Click a cell to plot its histories below. Grey = not run.</div></section>
<section class="card"><h2>Histories (select runs)</h2><select id="runs" multiple size="6" style="width:100%"></select>
<label>Quantity <select id="q"><option value="g">|Γ| at port</option><option value="ne">axis n_e at focus</option><option value="r50">channel r50 (µm)</option><option value="te">axis T_e (eV)</option><option value="N">total electrons</option></select></label>
<canvas id="hist" width="640" height="330"></canvas></section></div>
<section class="card"><h2>Comparison with experiment</h2><div class="tw" id="exp"></div></section>
<div class="grid" id="figs"></div>
<section class="card"><h2>All runs</h2><div class="tw" id="tab"></div></section>
<script>
const R=__ROWS__,H=__HIST__,X=__EXP__,F=__FIGS__;
const $=i=>document.getElementById(i),fmt=v=>v==null?'—':(Math.abs(v)>=1e4||Math.abs(v)<1e-2&&v!=0?v.toExponential(2):(+v).toPrecision(3));
const pe=R.filter(r=>r.root=='scan_pE'),P=[...new Set(pe.map(r=>Math.round(r.p_psig)))].sort((a,b)=>a-b),E=[...new Set(pe.map(r=>Math.round(r.E_kVcm)))].sort((a,b)=>a-b);
function col(f){const s=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]],x=Math.max(0,Math.min(.999,f))*4,i=Math.floor(x),u=x-i;return `rgb(${s[i].map((c,k)=>Math.round(c+(s[i+1][k]-c)*u))})`}
let cells=[];function map(){const c=$('map'),g=c.getContext('2d'),m=$('met').value,lg=/ne_|tg_/.test(m);g.clearRect(0,0,c.width,c.height);
const L=70,T=20,W=(c.width-L-20)/E.length,Hh=(c.height-T-50)/P.length,v=pe.map(r=>r[m]).filter(x=>x!=null).map(x=>lg?Math.log10(x):x),lo=Math.min(...v),hi=Math.max(...v);cells=[];
g.font='12px sans-serif';P.forEach((p,i)=>E.forEach((e,j)=>{const r=pe.find(r=>Math.round(r.p_psig)==p&&Math.round(r.E_kVcm)==e),x=L+j*W,y=T+(P.length-1-i)*Hh;
if(!r){g.fillStyle='#223';g.fillRect(x+1,y+1,W-2,Hh-2);return}const val=r[m],fr=val==null?0:((lg?Math.log10(val):val)-lo)/(hi-lo||1);g.fillStyle=val==null?'#555':col(fr);g.fillRect(x+1,y+1,W-2,Hh-2);
g.fillStyle=fr<.45?'#fff':'#000';g.textAlign='center';g.fillText(val==null?'no BD':fmt(val),x+W/2,y+Hh/2-2);g.fillText(r.E_N_Td.toFixed(0)+' Td',x+W/2,y+Hh/2+13);cells.push([x,y,W,Hh,r.run])}));
g.fillStyle='#e3edf6';g.textAlign='right';P.forEach((p,i)=>g.fillText(p+' psig',L-6,T+(P.length-1-i)*Hh+Hh/2+4));g.textAlign='center';E.forEach((e,j)=>g.fillText(e+' kV/cm',L+j*W+W/2,c.height-30));g.fillText('RF field at focus (peak)',c.width/2,c.height-10)}
$('map').onclick=ev=>{const c=$('map'),b=c.getBoundingClientRect(),x=(ev.clientX-b.left)*c.width/b.width,y=(ev.clientY-b.top)*c.height/b.height;
const h=cells.find(([a,b2,w,hh])=>x>a&&x<a+w&&y>b2&&y<b2+hh);if(h){[...$('runs').options].forEach(o=>{if(o.value==h[4])o.selected=!o.selected});hist()}};
const pal=['#7cc4ff','#ff8080','#9be37c','#ffd166','#c49bff','#ff9f43','#4dd0e1','#f06292'];
function hist(){const c=$('hist'),g=c.getContext('2d'),q=$('q').value,sel=[...$('runs').selectedOptions].map(o=>o.value),lg=q=='ne'||q=='N';g.clearRect(0,0,c.width,c.height);
const L=60,B=40,W=c.width-L-10,Hh=c.height-B-10;let xs=[],ys=[];sel.forEach(k=>{const h=H[k];const t=q=='g'?h.tg:h.t;h[q].forEach((v,i)=>{if(!lg||v>0){xs.push(t[i]);ys.push(lg?Math.log10(v):v)}})});
if(!xs.length){g.fillStyle='#a9bccd';g.fillText('select one or more runs',L+20,40);return}let x0=Math.min(...xs),x1=Math.max(...xs),y0=lg?Math.max(Math.min(...ys),Math.max(...ys)-10):Math.min(0,...ys),y1=Math.max(...ys);
const X_=t=>L+(t-x0)/(x1-x0||1)*W,Y_=v=>10+Hh-((v-y0)/(y1-y0||1))*Hh;g.strokeStyle='#294158';g.strokeRect(L,10,W,Hh);g.fillStyle='#e3edf6';g.font='11px sans-serif';
for(let i=0;i<=5;i++){const v=y0+(y1-y0)*i/5;g.textAlign='right';g.fillText(lg?'1e'+v.toFixed(1):fmt(v),L-4,Y_(v)+4);const t=x0+(x1-x0)*i/5;g.textAlign='center';g.fillText(t.toFixed(1),X_(t),c.height-B+16)}
g.fillText('time after laser peak (ns)',L+W/2,c.height-6);sel.forEach((k,n)=>{const h=H[k],t=q=='g'?h.tg:h.t;g.strokeStyle=pal[n%pal.length];g.lineWidth=1.6;g.beginPath();let s=0;
h[q].forEach((v,i)=>{if(lg&&v<=0)return;const y=lg?Math.log10(v):v;if(y<y0)return;s++?g.lineTo(X_(t[i]),Y_(y)):g.moveTo(X_(t[i]),Y_(y))});g.stroke();g.fillStyle=pal[n%pal.length];g.textAlign='left';g.fillText(R.find(r=>r.run==k).label,L+8,24+14*n)})}
R.forEach(r=>{const o=new Option(r.label,r.run);$('runs').add(o)});[...$('runs').options].slice(0,3).forEach(o=>o.selected=true);
$('runs').onchange=hist;$('q').onchange=hist;$('met').onchange=map;
$('exp').innerHTML='<table><tr><th>Quantity</th><th>Measured</th><th>Model</th><th>Comment</th></tr>'+X.map(r=>'<tr>'+r.map(c=>'<td>'+c+'</td>').join('')+'</tr>').join('')+'</table>';
$('figs').innerHTML=F.map(([t,s])=>`<section class="card"><h2>${t}</h2><img src="data:image/png;base64,${s}" alt="${t}"></section>`).join('');
const K=['label','gamma_end','delay03_ns','t_Ecollapse_ns','E_N_Td','fil_ne_axis_max','fil_tg_axis_max','r50_um_end','seed_electrons','closure_fil'];
$('tab').innerHTML='<table><tr>'+K.map(k=>'<th>'+k+'</th>').join('')+'</tr>'+R.map(r=>'<tr>'+K.map(k=>'<td>'+(k=='label'?r[k]:fmt(r[k]))+'</td>').join('')+'</tr>').join('')+'</table>';
map();hist();
</script></body></html>"""


def write_html(rows, xrows):
    import base64
    rr = []
    hist = {}
    for r in rows:
        d = {k: r.get(k) for k in KEEP}
        d["label"] = label(r)
        rr.append(d)
        hist[r["run"]] = _history(r)
    figs = [(t, base64.b64encode((OUT / f).read_bytes()).decode()) for t, f in (
        ("p × E: |Γ|, delay and 0D surrogate", "scan_pE.png"), ("Laser focus / pulse length / energy", "scan_laser.png"),
        ("Channel growth (3D filament)", "channel_growth.png"), ("Seed number and cavity anchors", "seed_anchor.png"),
        ("2D streamer: radius and length vs time", "../../streamer2d/figures/radius_length_vs_time.png"),
        ("2D streamer: diameter vs measured", "../../streamer2d/figures/diameter_vs_measured.png"),
        ("Jitter sensitivity (0D Monte Carlo)", "../../stochastic_seed/jitter_sensitivity.png")) if (OUT / f).exists()]
    js = lambda o: json.dumps(o, separators=(",", ":"), default=float)  # noqa: E731
    jn = lambda o: js(o).replace("NaN", "null")  # noqa: E731  (numbers only; never on base64 payloads)
    html = (PAGE.replace("__ROWS__", jn(rr)).replace("__HIST__", jn(hist))
            .replace("__EXP__", js(xrows)).replace("__FIGS__", js(figs)))
    (OUT / "scans.html").write_text(html, encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rates = Rates("He")
    pe, la, anc, conv = (rows_for(r) for r in ("scan_pE", "scan_laser", "anchor", "conv"))
    expj = load_exp()
    s = dict(pE=fig_pE(pe, rates))
    fig_laser(la, pe + conv)
    fig_channel(pe, la)
    nominal = next((r for r in pe if round(r["p_psig"]) == 100 and round(r["E_kVcm"]) == 80), None)
    s["anchor"] = fig_seed(conv, anc, rates, nominal)
    s["seed_scan"] = [dict(run=r["run"], N0=r["seed_electrons"], delay03_ns=r["delay03_ns"],
                           gamma_end=r["gamma_end"]) for r in conv + ([nominal] if nominal else [])]
    s["threshold"] = sorted(
        [dict(p_psig=round(r["p_psig"]), E_kVcm=r["E_kVcm"], E_N_Td=round(r["E_N_Td"], 1),
              breakdown=r["delay03_ns"] is not None) for r in pe], key=lambda d: (d["p_psig"], d["E_kVcm"]))
    s["closures_max"] = max(abs(r["closure_fil"] or 0) for r in pe + la + anc + conv)
    s["experimental_json_loaded"] = expj is not None
    (OUT / "scan_summary.json").write_text(json.dumps(s, indent=1, default=float))
    xr = exp_rows(s, pe, anc)
    write_md(xr)
    write_html(pe + la + conv + anc, xr)
    print(json.dumps(s, indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
