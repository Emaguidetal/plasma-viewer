"""Unified 3D visualisation of the laser-induced plasma switch in the WR-284 cavity.

One time slider from the laser pulse (-1.5 ns) to the end of the chained drive run (0.4 us with the RF switched
off 100 ns after the laser, or 1 us with the RF kept on), stitched from the executed models:
  * laser pulse and seed: Gaussian 1 ns FWHM pulse centred at t = 0. rod3d deposits the whole seed at t = 0, so
    the seed shown before t = 0 is the cumulative I^k fraction (k = 6 or 4 for the rod3d seed model; display only);
  * avalanche and rod growth along E: rod3d frames (thin-wire tree model, cavity energy feedback);
  * bridged rod under drive: rod_drive chained run started from the rod3d hand-off (results/rod_drive/chain_*).
The RF cones are the TE101 modal envelope |Ey| = E(t) sin(pi x/a) sin(pi z/d) scaled by the modal amplitude
(rod3d: E0; rod_drive: |V|/b at the post plane); the local field around the channels and the post-loaded mode
shape are not drawn. Nothing here is measured.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .rod3d import A_WG, B_WG, D_CAV

ROOT = Path(__file__).resolve().parents[2] / "results"
OUT = ROOT / "story3d"
TAU_L = 1.0e-9                      # laser FWHM (s)
NE_MIN = 1e14                       # display threshold (keeps the ~1e15 m^-3 laser seed visible)
ZOOM = 4.0                          # mm, half width of the side view around the seeded region
Z_R_MM = 2.36                       # vacuum Rayleigh range, w0 = 20 um at 532 nm


def _r_at(r, prof, thr):
    """Outermost radius where prof >= thr (prof decreasing outward); 0 if never."""
    idx = np.nonzero(prof >= thr)[0]
    return float(r[idx[-1]]) if len(idx) else 0.0


def _tg_color(tg):
    """False colour for the axis gas temperature 300 K -> 15 kK (not an emission prediction)."""
    s = min(max((math.log(max(tg, 300)) - math.log(300)) / (math.log(15000) - math.log(300)), 0.0), 1.0)
    stops = [(0, (60, 60, 160)), (0.35, (220, 40, 20)), (0.7, (255, 200, 40)), (1.0, (255, 255, 255))]
    for (s0, c0), (s1, c1) in zip(stops[:-1], stops[1:]):
        if s <= s1:
            f = (s - s0) / (s1 - s0)
            return "rgb(%d,%d,%d)" % tuple(int(a + f * (b - a)) for a, b in zip(c0, c1))
    return "rgb(255,255,255)"


def _drive_phase(on, t, ds, r_ohm, x_ohm, q_l, tg):
    if not on:
        return "RF off: afterglow, recombination of the rod"
    tc = ds.get("t_contraction_s")
    ts = ds.get("t_skin_invalid_s")
    if t < 20e-9 and q_l < 20:
        ph = f"rod loads the cavity: modal field collapses (rod-loaded Q = {q_l:.0f})"
    elif tc is None or t < tc:
        ph = (f"RF on: diffuse glow column, gas still cool (Tg = {tg / 1e3:.1f} kK)" if tg < 2000 else
              f"RF on: diffuse glow column, gas heating toward contraction (Tg = {tg / 1e3:.1f} kK)")
    else:
        ph = "RF on: thermally contracted arc channel"
        if r_ohm < x_ohm:
            ph += "; R below the post reactance X (reactance-limited post)"
    if ts is not None and t >= ts:
        ph += " [outside model validity: skin depth < 2 r_c]"
    return ph


def case_data(rod_case, chain_case=None, n_drive_frames=60):
    with np.load(ROOT / "rod3d" / rod_case / "rod3d.npz") as z:
        h = {k: z[k] for k in z.files}
    s = json.loads((ROOT / "rod3d" / rod_case / "summary.json").read_text(encoding="utf-8"))
    cfg = s["config"]
    e_nom = cfg["e0_v_m"] / 1e5
    pos, rad, kind, rod, lam, t3 = h["pos"], h["rad"], h["kind"], h["rod"], h["lam"], h["t_s"]
    xc, zc = A_WG / 2, D_CAV / 2
    line = np.nonzero(kind == 0)[0]
    seed = h["seed_electrons"]
    k_ord = 4 if cfg.get("seed_model") == "organic" else 6
    sig_k = TAU_L / math.sqrt(8 * k_ord * math.log(2))    # std of I^k for a Gaussian pulse
    ne_seed = seed / (cfg["dx_c"] * math.pi * cfg["r_s"] ** 2)
    x_focus = cfg["x_focus_m"] * 1e3
    frames, hist = [], dict(t=[], e=[], ea=[], ne=[], tg=[])

    def laser(t):
        return math.exp(-4 * math.log(2) * (t / TAU_L) ** 2)

    def nodes_frame(k, keep_rods=None):
        m = int(h["m_frame"][k])
        ne = lam[k, :m] / (math.pi * np.maximum(rad[:m], 1e-9) ** 2)
        keep = ne > NE_MIN
        if keep_rods is not None:
            keep &= np.isin(rod[:m], keep_rods) | (kind[:m] == 0)
        p = pos[:m][keep]
        return dict(x=np.round(1e3 * (p[:, 0] - xc), 3).tolist(), y=np.round(1e3 * p[:, 1], 3).tolist(),
                    z=np.round(1e3 * (p[:, 2] - zc), 3).tolist(), c=np.round(np.log10(ne[keep]), 2).tolist()), ne

    # 1) laser pulse before t = 0
    for t in np.array([-1.5, -1.2, -0.9, -0.7, -0.5, -0.35, -0.2, -0.1]) * 1e-9:
        f = 0.5 * (1 + math.erf(t / (sig_k * math.sqrt(2))))
        ok = seed * f >= 0.5
        frames.append(dict(t=t * 1e9, phase="laser pulse: seed buildup is display-only; model deposits seed at t=0",
                           e=e_nom, il=laser(t),
                           pts=dict(x=np.round(1e3 * (pos[line[ok], 0] - xc), 3).tolist(),
                                    y=[round(1e3 * B_WG / 2, 3)] * int(ok.sum()), z=[0.0] * int(ok.sum()),
                                    c=np.round(np.log10(ne_seed[ok] * f), 2).tolist(),
                                    s=np.round(2 + 2 * np.log10(seed[ok] * f), 1).tolist()),
                           rod=None))
        hist["t"].append(t * 1e9); hist["e"].append(e_nom); hist["ea"].append(None)
        hist["ne"].append(float(np.log10(max(ne_seed.max() * f, 1.0)))); hist["tg"].append(None)
    # 2) rod3d avalanche and growth (until the bridge if a chained drive run follows)
    tb = (s["t_bridge_ns"] or 1e9) * 1e-9
    stop = tb if chain_case else t3[-1]
    last_k = 0
    for k in range(len(t3)):
        if t3[k] > stop + 1e-15 or t3[k] <= 0:
            continue
        fr, ne = nodes_frame(k)
        if t3[k] < (s["t_first_birth_ns"] or 1e9) * 1e-9:
            ph = "seed deposited (rod3d instant-seed approximation); avalanche on the laser line"
        elif s["t_bridge_ns"] is None:
            ph = "rods evolve along E; no wall-to-wall bridge in this run"
        else:
            ph = "rods grow along E toward both broad walls"
        frames.append(dict(t=float(t3[k] * 1e9), phase=ph, e=float(h["hist_e0_v_m"][k] / 1e5),
                           il=laser(t3[k]), pts=fr, rod=None))
        hist["t"].append(float(t3[k] * 1e9)); hist["e"].append(float(h["hist_e0_v_m"][k] / 1e5))
        hist["ea"].append(None); hist["ne"].append(float(np.log10(max(ne.max(), 1.0)))); hist["tg"].append(None)
        last_k = k
    info = dict(case=rod_case, e_nom=e_nom, t_birth=s["t_first_birth_ns"], t_bridge=s["t_bridge_ns"],
                n_rods=s["n_rods"], x_view=x_focus, seed_model=cfg.get("seed_model", "power6"),
                n_seed=float(seed.sum()), chain=None)
    # 3) bridged rod under drive (chained run)
    if chain_case and (ROOT / "rod_drive" / chain_case / "rod_drive.npz").exists():
        with np.load(ROOT / "rod_drive" / chain_case / "rod_drive.npz") as z:
            d = {k: z[k] for k in z.files}
        ds = json.loads((ROOT / "rod_drive" / chain_case / "summary.json").read_text(encoding="utf-8"))
        r = d["r_m"]
        m = int(h["m_frame"][last_k])
        ext = {q: np.ptp(pos[:m][rod[:m] == q, 1]) for q in np.unique(rod[:m][rod[:m] >= 0])}
        rb = max(ext, key=ext.get)
        x_rod = float(np.mean(pos[:m][rod[:m] == rb, 0]) - xc) * 1e3
        z_rod = float(np.mean(pos[:m][rod[:m] == rb, 2]) - zc) * 1e3
        others, _ = nodes_frame(last_k, keep_rods=[q for q in ext if q != rb])
        none = dict(x=[], y=[], z=[], c=[])
        td = d["t_s"]
        t_off = ds["config"]["drive_off_s"]
        want = list(np.geomspace(0.1e-9, td[-1], n_drive_frames))
        for tk in (ds.get("t_contraction_s"), t_off if t_off >= 0 else None):
            if tk:                                     # dense frames around contraction and drive-off
                want += list(tk + np.linspace(-30e-9, 40e-9, 15))
        sel = np.unique(np.clip(np.searchsorted(td, [w for w in want if 0 < w <= td[-1]]), 1, len(td) - 1))
        for i in sel:
            ne_ax = float(d["ne_axis_m3"][i]); tg = float(d["tg_axis_k"][i])
            r_e = _r_at(r, d["prof_ne"][i], 1e21)
            hot = np.nonzero(d["prof_rho"][i] < 0.7)[0]
            r_hot = float(r[hot[-1]]) if len(hot) else 0.0
            on = t_off < 0 or td[i] < t_off
            q_l = float(d["q_loaded"][i]) if "q_loaded" in d else float("nan")
            ph = _drive_phase(on, td[i], ds, d["r_ohm"][i], d["x_ohm"][i], q_l, tg)
            lne = math.log10(max(ne_ax, 1.0))
            t_ns = tb * 1e9 + float(td[i]) * 1e9
            frames.append(dict(t=t_ns, phase=ph, e=float(d["e_gap_v_m"][i] / 1e5), il=0.0,
                               pts=others if td[i] < 20e-9 else none,
                               rod=dict(x=x_rod, z=z_rod, ne=round(lne, 2), tg=round(tg),
                                        te=round(float(d["te_axis_ev"][i]), 2), r_e=round(r_e * 1e3, 4),
                                        r_hot=round(r_hot * 1e3, 4), rc=round(float(d["r_cond_m"][i]) * 1e3, 4),
                                        col=_tg_color(tg), op=round(min(max((lne - 18) / 4, 0.15), 1.0), 2),
                                        p_kw=round(float(d["p_rod_w"][i]) / 1e3, 1), rohm=round(float(d["r_ohm"][i])),
                                        xohm=round(float(d["x_ohm"][i])), q=round(q_l, 1),
                                        ea=round(float(d["e_axial_v_m"][i]) / 1e5, 2))))
            hist["t"].append(t_ns); hist["e"].append(float(d["e_gap_v_m"][i] / 1e5))
            hist["ea"].append(float(d["e_axial_v_m"][i] / 1e5))
            hist["ne"].append(lne); hist["tg"].append(tg / 1e3)
        info["x_view"] = x_rod
        tc = ds.get("t_contraction_s")
        info["chain"] = dict(name=chain_case, e_start=ds["config"]["e_start_v_m"] / 1e5,
                             n0=ds["config"]["n0_m3"], r0_mm=ds["config"]["r0_m"] * 1e3,
                             p_fwd_mw=ds["p_forward_w"] / 1e6, w_cav0_mj=ds["w_cav0_j"] * 1e3,
                             t_off_ns=(tb + t_off) * 1e9 if t_off >= 0 else None,
                             t_contraction_ns=(tb + tc) * 1e9 if tc is not None else None,
                             t_r_below_x_ns=(tb + ds["t_r_below_x_s"]) * 1e9 if ds.get("t_r_below_x_s") else None,
                             t_skin_invalid_ns=(tb + ds["t_skin_invalid_s"]) * 1e9 if ds.get("t_skin_invalid_s")
                             else None,
                             e_min_20ns=ds.get("e_gap_min_first_20ns_kv_cm"),
                             tg_max=ds["tg_axis_max_k"], r_end=ds["r_end_ohm"],
                             resid=ds["circuit_energy_residual_rel"])
    for f in frames:
        f["t"] = round(f["t"], 3); f["e"] = round(f["e"], 2); f["il"] = round(f["il"], 3)
    return dict(info=info, frames=frames, hist=hist)


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Laser-induced plasma switch in 3D: WR-284 cavity, 1.48 MPa He</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
*{box-sizing:border-box}body{margin:0;background:#0b111b;color:#e5edf5;font:15px/1.5 system-ui,Segoe UI,Arial,sans-serif}
main{max-width:1500px;margin:auto;padding:16px 20px 36px}a{color:#8ed2ff}h1,h2,p{margin:0}
nav{font-size:.85rem;margin-bottom:8px}h1{font-size:clamp(1.45rem,3vw,2.1rem);line-height:1.15}
.subtitle,.hint{color:#aabbd0}.subtitle{margin:5px 0 12px}.notice{border-left:3px solid #ffc66d;background:#302818;padding:9px 12px;margin:12px 0;color:#ffdfaa}
.controls,.card,details{background:#121d2b;border:1px solid #29394c;border-radius:10px}
.controls{padding:12px;margin-bottom:14px}.controls label{display:block;font-weight:600;margin-bottom:4px}
select,button,input{font:inherit}select{width:100%;max-width:100%;min-width:0;background:#0b1420;color:#eef6ff;border:1px solid #50647b;border-radius:6px;padding:8px}
.playback{display:flex;align-items:center;gap:12px;margin-top:10px}button{background:#7bd2ff;color:#07141e;border:0;border-radius:6px;padding:8px 18px;font-weight:700;cursor:pointer}
button:focus-visible,select:focus-visible,input:focus-visible,a:focus-visible,summary:focus-visible{outline:2px solid #ffe08a;outline-offset:2px}
input[type=range]{flex:1;min-width:0;accent-color:#7bd2ff}.readout{margin-top:8px;display:flex;gap:8px 18px;flex-wrap:wrap}
#tt{color:#d9edfc}#ph{color:#ffd27a;font-weight:600}#st{color:#c6e8df;overflow-wrap:anywhere}
#error{color:#ffd5d5;background:#502a31;padding:10px;border-radius:6px;margin:10px 0}#error:empty{display:none}
.viewer-grid{display:grid;grid-template-columns:minmax(0,1.35fr) minmax(0,1fr);gap:14px}
.card{min-width:0;overflow:hidden;padding:12px}.scene-card{grid-row:span 2}.card h2{font-size:1rem}
.hint{font-size:.84rem;margin:3px 0 0}.plot{width:100%;min-width:0}#scenePlot{height:590px}#sidePlot{height:260px}#historyPlot{height:245px}
.legend{display:flex;flex-wrap:wrap;gap:4px 14px;color:#c5d3e2;font-size:.8rem;margin:3px 0}
.legend i{display:inline-block;width:15px;height:3px;vertical-align:middle;margin-right:5px;background:var(--swatch)}
details{padding:10px 12px;margin-top:12px}summary{cursor:pointer;font-weight:650}details p{margin-top:8px;color:#c1cedc}
@media(max-width:850px){main{padding:12px}.viewer-grid{grid-template-columns:minmax(0,1fr)}.scene-card{grid-row:auto}
 #scenePlot{height:390px}#sidePlot{height:270px}#historyPlot{height:250px}}
@media(max-width:360px){main{padding:10px}.card{padding:9px}#scenePlot{height:340px}}
</style></head><body><main>
<nav><a href="https://emaguidetal.github.io/plasma-viewer/">Research home</a> / <a href="https://emaguidetal.github.io/plasma-viewer/switch/">Switch study</a></nav>
<h1>3D plasma switch</h1>
<p class="subtitle">Laser-triggered WR-284 &middot; 1.48 MPa helium &middot; 2.95 GHz</p>
<p class="notice"><strong>Simulation, not measurement.</strong> Fields and seed are uncertain. <a href="#interpretation">Limitations</a>.</p>
<section class="controls" aria-label="Simulation controls">
 <label for="cs">Scenario</label><select id="cs"></select>
 <div class="playback"><button id="pl" type="button" aria-pressed="false">Play</button>
  <input type="range" id="sl" min="0" value="0" aria-label="Simulation time"></div>
 <div class="readout"><span id="tt"></span><span id="ph"></span></div>
</section>
<div id="error" role="alert"></div><noscript><p class="notice">JavaScript is required to display the simulation.</p></noscript>
<section class="viewer-grid" aria-label="Simulation views">
 <article class="card scene-card"><h2>3D cavity</h2>
  <p class="hint">Drag to rotate. Blue: RF envelope, not carrier cycles. Green: laser.</p>
  <div id="scenePlot" class="plot" aria-label="Interactive three-dimensional cavity and plasma"></div></article>
 <article class="card"><h2>Side view</h2>
  <p class="hint">The laser crosses x; plasma grows along y. Horizontal scale is stretched.</p>
  <div id="sidePlot" class="plot" aria-label="Side view of the seeded region"></div></article>
 <article class="card"><h2>Evolution</h2>
  <p class="hint">Dashed vertical line: bridge. White marker: current frame. Time is shifted by 2 ns for the log axis.</p>
  <div id="historyPlot" class="plot" aria-label="Model history over time"></div>
  <div class="legend"><span><i style="--swatch:#77ccff"></i>modal field</span><span><i style="--swatch:#a9e3ff"></i>rod axial field (dotted)</span>
   <span><i style="--swatch:#ffaa44"></i>electron density (scaled)</span><span><i style="--swatch:#ff6666"></i>gas temperature (scaled)</span></div>
  <p class="hint">Only field uses the kV/cm axis; density and temperature are scaled to show their trends.</p></article>
</section>
<details><summary>Current rod properties (model)</summary><p id="st">No bridged rod in this frame.</p></details>
<details id="interpretation"><summary>What is shown and model limitations</summary>
 <p>Blue cones represent the TE101 modal envelope, not individual RF cycles or the local field around a channel. The green line is the laser; the marker and bar show the assumed focus and &plusmn;2.4 mm Rayleigh range. Plasma points are coloured by log10 electron density (m<sup>-3</sup>).</p>
 <p>After bridging, the rod colour represents axis gas temperature, <strong>not emitted light</strong>. Its width is a display contour at n<sub>e</sub> &gt; 10<sup>21</sup> m<sup>-3</sup>; the dashed side-view outline marks the rarefied core. Partial rods are shown in grey at their bridge-time state only, not as later calculated evolution.</p>
 <p>The model prescribes a weak seed (about 10<sup>15</sup> m<sup>-3</sup> in this example); it is not a measurement of the laser seed. Pre-zero seed buildup is display-only: the rod-growth model deposits the entire seed at t=0. The laser alone produces no visible plasma in the reported experiment, but that does not prove that no electrons are present.</p>
 <p>Nominal fields and power are not independently calibrated. An illustrative 3.91 field multiplier was fitted to an assumed breakdown probability, not measured as a field calibration; a lower-power sensitivity scales nominal power by 1/3.91<sup>2</sup>. The reported approximately 10% RF-only breakdown has no shot count. Contraction and afterglow after the skin-depth warning are not quantitatively reliable.</p>
 <p>The side view uses a stretched x scale. The history scales log10 density from 14 to 25 and rod-axis temperature from 0 to 15 kK onto the field axis; density is the maximum node before bridging and the rod axis afterwards. See the <a href="https://emaguidetal.github.io/plasma-viewer/research/literature/experimental-data-krasik-leopold.md">experimental-data review</a> for the distinct published 0.2 MPa helium benchmark.</p>
</details>
<script>
const ALL=__DATA__; const A=__A__, B=__B__, D=__D__, ZM=__ZOOM__, ZR=__ZR__;
const cs=document.getElementById('cs'),sl=document.getElementById('sl'),play=document.getElementById('pl');
Object.keys(ALL).forEach(k=>{const o=document.createElement('option');o.value=k;o.textContent=ALL[k].label;cs.appendChild(o);});
let C=null,timer=null;
const gx=[],gy=[],gz=[],sh=[];
for(let i=1;i<=5;i++)for(let j=1;j<=5;j++)for(let k=0;k<3;k++){const x=-A/2+A*i/6,z=-D/2+D*j/6;gx.push(x);gz.push(z);gy.push(B*(k+1)/4);
 sh.push(Math.sin(Math.PI*(x+A/2)/A)*Math.sin(Math.PI*(z+D/2)/D));}
function walls(){const a=A/2,d=D/2;const X=[],Y=[],Z=[];
 for(const y of [0,B]){X.push(-a,a,a,-a,-a,null);Z.push(-d,-d,d,d,-d,null);Y.push(y,y,y,y,y,null);}
 for(const [x,z] of [[-a,-d],[a,-d],[a,d],[-a,d]]){X.push(x,x,null);Z.push(z,z,null);Y.push(0,B,null);}
 return {type:'scatter3d',mode:'lines',x:X,y:Z,z:Y,line:{color:'#8295a8',width:3},name:'cavity walls',hoverinfo:'skip'};}
const W=walls();
function traces(f){const en=C.info.e_nom,xf=C.info.x_view,fo=C.info.focus;
 const u=sh.map(()=>0),w=sh.map(s=>f.e*s);
 const cone={type:'cone',x:gx,y:gz,z:gy,u:u,v:u,w:w,sizemode:'raw',sizeref:6/en,anchor:'center',
  colorscale:[[0,'#164c78'],[1,'#78d2ff']],cmin:0,cmax:en,showscale:false,name:'RF modal envelope',
  hovertemplate:'Envelope %{w:.0f} kV/cm<extra></extra>',opacity:0.55};
 const las={type:'scatter3d',mode:'lines',x:[-A/2,A/2],y:[0,0],z:[B/2,B/2],name:'laser',
  line:{color:f.il>0.02?'#46ef87':'#32694d',width:1+12*f.il,dash:f.il>0.02?'solid':'dot'},hoverinfo:'skip'};
 const foc={type:'scatter3d',mode:'lines+markers',x:[fo-ZR,fo,fo+ZR],y:[0,0,0],z:[B/2,B/2,B/2],name:'focus and Rayleigh range',
  line:{color:'#a6ffbf',width:7},marker:{size:[0,5,0],color:'#a6ffbf'},hoverinfo:'skip'};
 const P=f.pts,part=!!f.rod;
 const pts={type:'scatter3d',mode:'markers',x:P.x,y:P.z,z:P.y,name:part?'partial rods at bridging':'plasma / seed',
  marker:part?{size:3,color:'#999',opacity:0.35}:{size:P.s||3,color:P.c,cmin:14,cmax:24.5,colorscale:'Hot',
   colorbar:{title:'log10 ne',len:0.45,x:0.96,y:0.75,thickness:10}},
  hovertemplate:part?'partial rod at bridging<extra></extra>':'log10 ne %{marker.color}<extra></extra>'};
 const scene=[W,cone,las,foc,pts];
 if(f.rod){const R=f.rod;scene.push({type:'scatter3d',mode:'lines',x:[R.x,R.x],y:[R.z,R.z],z:[0,B],name:'bridged rod: axis gas temperature',
  opacity:R.op,line:{color:R.col,width:Math.max(4,Math.min(26,R.r_e*60))},hoverinfo:'skip'});}
 const zx=[],zy=[],zc=[];P.x.forEach((x,i)=>{if(Math.abs(x-xf)<=ZM){zx.push(x);zy.push(P.y[i]);zc.push(P.c[i]);}});
 const side=[{type:'scatter',mode:'lines',x:[xf-ZM,xf+ZM],y:[B/2,B/2],
  line:{color:f.il>0.02?'#46ef87':'#32694d',width:1+6*f.il,dash:'dot'},hoverinfo:'skip'},
  {type:'scatter',mode:'markers',x:zx,y:zy,marker:part?{size:4,color:'#888'}:
   {size:6,color:zc,cmin:14,cmax:24.5,colorscale:'Hot',line:{color:'#fff',width:0.3}},hoverinfo:'skip'}];
 if(f.rod){const R=f.rod;
  if(R.r_hot>0)side.push({type:'scatter',mode:'lines',x:[R.x-R.r_hot,R.x+R.r_hot,R.x+R.r_hot,R.x-R.r_hot,R.x-R.r_hot],y:[0,0,B,B,0],
   line:{color:'#6cf',dash:'dash',width:1},hoverinfo:'skip'});
  const re=Math.max(R.r_e,0.01);side.push({type:'scatter',mode:'lines',fill:'toself',fillcolor:R.col,opacity:R.op,
   x:[R.x-re,R.x+re,R.x+re,R.x-re,R.x-re],y:[0,0,B,B,0],line:{color:R.col,width:0.5},hoverinfo:'skip'});}
 const H=C.hist,tx=H.t.map(t=>t+2);
 const history=[
  {type:'scatter',mode:'lines',x:tx,y:H.e,name:'Modal field',line:{color:'#77ccff'}},
  {type:'scatter',mode:'lines',x:tx,y:H.ea,name:'Rod axial field',line:{color:'#a9e3ff',dash:'dot'}},
  {type:'scatter',mode:'lines',x:tx,y:H.ne.map(v=>v===null?null:Math.max(v-14,0)*en/11),name:'Electron density (scaled)',line:{color:'#ffaa44'}},
  {type:'scatter',mode:'lines',x:tx,y:H.tg.map(v=>v===null?null:v*en/15),name:'Gas temperature (scaled)',line:{color:'#ff6666'}}];
 if(C.info.t_bridge)history.push({type:'scatter',mode:'lines',x:[C.info.t_bridge+2,C.info.t_bridge+2],y:[0,en],
  line:{color:'#aaa',dash:'dash',width:1},hoverinfo:'skip'});
 history.push({type:'scatter',mode:'markers',x:[f.t+2],y:[f.e],marker:{color:'#fff',size:9,symbol:'diamond'},hoverinfo:'skip'});
 return {scene,side,history};}
function layouts(){const xf=C.info.x_view,base={paper_bgcolor:'#121d2b',plot_bgcolor:'#121d2b',font:{color:'#d8e4ef',size:11},showlegend:false,uirevision:'keep'};
 return {
  scene:{...base,margin:{l:0,r:25,t:0,b:0},scene:{xaxis:{title:'x (mm)',range:[-A/2,A/2]},yaxis:{title:'z (mm)',range:[-D/2,D/2]},
   zaxis:{title:'y (mm)',range:[0,B]},aspectmode:'data',camera:{eye:{x:1.3,y:1.5,z:0.6}}}},
  side:{...base,margin:{l:48,r:12,t:8,b:43},xaxis:{range:[xf-ZM,xf+ZM],title:'x (mm)',gridcolor:'#344456'},
   yaxis:{range:[0,B],title:'y (mm)',gridcolor:'#344456'}},
  history:{...base,margin:{l:48,r:12,t:8,b:43},xaxis:{type:'log',title:'time + 2 ns (log)',gridcolor:'#344456'},
   yaxis:{title:'kV/cm',gridcolor:'#344456'}}};}
const config={responsive:true,displaylogo:false};
function pause(){if(timer){clearInterval(timer);timer=null;}play.textContent='Play';play.setAttribute('aria-pressed','false');}
function show(k){const f=C.frames[k],ch=C.info.chain,R=f.rod;
 document.getElementById('tt').textContent='t = '+f.t.toFixed(2)+' ns | modal field = '+f.e.toFixed(0)+' kV/cm'+(ch?' | nominal forward power = '+ch.p_fwd_mw.toFixed(2)+' MW':'');
 document.getElementById('ph').textContent=f.phase;
 document.getElementById('st').textContent=R?('axis ne = 1e'+R.ne+' m^-3; Te = '+R.te+' eV; Tg = '+R.tg+' K; axial field = '+R.ea+' kV/cm; display radius = '+(R.r_e*1e3).toFixed(0)+' um; conducting radius = '+(R.rc*1e3).toFixed(0)+' um; R = '+R.rohm+' ohm; X = '+R.xohm+' ohm; rod-loaded Q = '+R.q+'; rod power = '+R.p_kw+' kW'):'No bridged rod in this frame.';
 const t=traces(f),l=layouts();
 Promise.all([Plotly.react('scenePlot',t.scene,l.scene,config),Plotly.react('sidePlot',t.side,l.side,config),
  Plotly.react('historyPlot',t.history,l.history,config)]).catch(err=>{
   pause();document.getElementById('error').textContent='Could not render the simulation: '+(err.message||String(err));});
}
function load(name){pause();C=ALL[name];sl.max=C.frames.length-1;sl.value=0;show(0);}
cs.onchange=()=>load(cs.value);sl.oninput=()=>show(+sl.value);
play.onclick=()=>{if(timer){pause();return;}play.textContent='Pause';play.setAttribute('aria-pressed','true');
 timer=setInterval(()=>{sl.value=(+sl.value+1)%C.frames.length;show(+sl.value);},220);};
if(typeof Plotly==='undefined'){document.getElementById('error').textContent='The Plotly library could not be loaded. Check your connection and reload the page.';play.disabled=true;sl.disabled=true;}
else if(!Object.keys(ALL).length){document.getElementById('error').textContent='No simulation cases were available when this viewer was built.';play.disabled=true;sl.disabled=true;}
else{load(Object.keys(ALL)[0]);}
</script></main></body></html>"""


CASES = (("E200", "chain_E200_cw", "200 kV/cm: continuous RF"),
         ("E200", "chain_E200", "200 kV/cm: RF off at 100 ns"),
         ("E200", "chain_E200_r100", "200 kV/cm: narrow rod, RF off at 100 ns"),
         ("E300", "chain_E300_cw", "300 kV/cm: continuous RF"),
         ("E300", "chain_E300", "300 kV/cm: RF off at 100 ns"),
         ("E100", None, "100 kV/cm: stalled rods"),
         ("E100_f20", None, "100 kV/cm: focus 20 mm off centre"))


def build(cases=CASES, out=OUT):
    out.mkdir(parents=True, exist_ok=True)
    data = {}
    for rc, cc, label in cases:
        if not (ROOT / "rod3d" / rc / "rod3d.npz").exists():
            continue
        d = case_data(rc, cc)
        d["label"] = label
        d["info"]["focus"] = json.loads((ROOT / "rod3d" / rc / "summary.json").read_text(
            encoding="utf-8"))["config"]["x_focus_m"] * 1e3
        data[cc or rc] = d
    html = (HTML.replace("__DATA__", json.dumps(data, separators=(",", ":")))
            .replace("__A__", f"{A_WG * 1e3:.3f}").replace("__B__", f"{B_WG * 1e3:.3f}")
            .replace("__D__", f"{D_CAV * 1e3:.3f}").replace("__ZOOM__", str(ZOOM)).replace("__ZR__", str(Z_R_MM)))
    p = out / "laser_plasma_3d.html"
    p.write_text(html, encoding="utf-8")
    summ = {k: v["info"] for k, v in data.items()}
    (out / "story_cases.json").write_text(json.dumps(summ, indent=1, default=float), encoding="utf-8")
    return p


if __name__ == "__main__":
    print(build())
