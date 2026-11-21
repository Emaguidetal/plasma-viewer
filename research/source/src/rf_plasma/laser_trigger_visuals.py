"""Viewer for the coupled laser-triggered RF breakdown runs (results/coupled3d/laser_trigger_He_1MPa).

Everything drawn comes from coupled3d.run_case output (frames.npz, timeseries.npz, rf_series.npz,
summary.json): the FDTD Ey field, the bulk electron density, the subgrid laser filament
(PPT photoionization + avalanche + He+/He2+ + hydro) and the port/probe RF diagnostics.
"""
import argparse
import base64
import gzip
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation, colors

from .coupled3d import RESULTS

ROOT = RESULTS / "laser_trigger_He_1MPa"
RUNS = (("rfonly", "rf-only", "RF only (no laser)"),
        ("L1", "rf+laser_500ps", "RF + 1 mJ / 500 ps laser"),
        ("L10", "rf+laser_3000ps", "RF + 10 mJ / 3 ns laser"))
GEOMS = {"cavity": "WR-284 iris cavity (standing wave from the start)",
         "guide": "matched WR-284 guide (traveling wave until the plasma reflects it)"}
DATASETS = {"cavity_E8": ("cavity", 8, "Iris cavity · 8 MV/m (above RF self-breakdown)"),
            "guide_E8": ("guide", 8, "Matched guide · 8 MV/m (above RF self-breakdown)"),
            "guide_E6": ("guide", 6, "Matched guide · 6 MV/m (below RF self-breakdown in 24 ns)")}
MAX_FRAMES = 270
NE_LOG = (9.0, 25.0)       # log10 m^-3 colour range for densities


def run_dir(geom, key, root=ROOT, field=8):
    return root / f"{geom}_He_{key}_E{field:g}MVm"


def _b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")


def _q_field(a, scale):
    return np.rint(a / scale * 127).clip(-127, 127).astype(np.int8)


def _q_log(a):
    lg = np.log10(np.maximum(a, 10 ** NE_LOG[0]))
    return np.rint((lg - NE_LOG[0]) / (NE_LOG[1] - NE_LOG[0]) * 255).clip(0, 255).astype(np.uint8)


def _focus_edge(fr):
    fj = int(fr["focus"][1])
    return int(np.argmin(np.abs(fr["fil_y"] / float(fr["dx"]) - 0.5 - fj)))


def _envelope(ey_line, t, f, window):
    """RF amplitude per frame: least-squares fit of a cos + b sin (drive frequency f) over the last
    `window` frames (about one RF period); removes the sampling ripple of a running max."""
    a = np.asarray(ey_line, float)
    out = np.abs(a).copy()
    w = 2 * math.pi * f
    for m in range(2, len(a)):
        i0 = max(0, m - window + 1)
        tt = t[i0:m + 1]
        basis = np.stack([np.cos(w * tt), np.sin(w * tt)], axis=1)
        coef, *_ = np.linalg.lstsq(basis, a[i0:m + 1].reshape(len(tt), -1), rcond=None)
        out[m] = np.hypot(coef[0], coef[1]).reshape(a.shape[1:])
    return out


def _rf_hist(rs, f_drive, n_out=600):
    dt = float(rs["dt"])
    per = int(round(1 / (f_drive * dt)))
    ey = np.asarray(rs["ey_focus"], float)
    nper = len(ey) // per
    pk = np.abs(ey[:nper * per]).reshape(nper, per).max(axis=1)
    t_pk = (np.arange(nper) + 1) * per * dt
    fwd, bwd = np.asarray(rs["fwd"]), np.asarray(rs["bwd"])
    gam = np.abs(bwd) / np.maximum(np.abs(fwd), 1e-30)
    pf = np.abs(fwd) ** 2
    tr = np.asarray(rs["p_out"], float)
    sel = np.linspace(0, len(gam) - 1, min(n_out, len(gam))).astype(int)
    ref = np.median(pf[: max(per, len(pf) // 10)])
    p_ref = float(np.median(tr[: max(per, len(tr) // 10)])) if tr.size else 0.0
    return dict(t_pk=t_pk, pk=pk, t_g=np.asarray(rs["t_phasor"])[sel], gamma=gam[sel],
                trans=(tr[sel] / p_ref) if p_ref > 0 else np.zeros(len(sel)), pf=pf[sel] / max(ref, 1e-30))


def load_geometry(geom, root=ROOT, field=8):
    runs = {}
    for key, tag, label in RUNS:
        d = run_dir(geom, tag, root, field)
        if not (d / "frames.npz").exists():
            continue
        s = json.loads((d / "summary.json").read_text(encoding="utf-8"))
        runs[key] = dict(label=label, fr=dict(np.load(d / "frames.npz")), ts=dict(np.load(d / "timeseries.npz")),
                         rs=dict(np.load(d / "rf_series.npz")), s=s)
    return runs


def _milestones(run):
    ts, s = run["ts"], run["s"]
    t = ts["t"]
    ne = ts["fil_ne_axis_focus"]
    out = {}
    for lab, thr in (("t_fil_1e17", 1e17), ("t_fil_1e20", 1e20), ("t_fil_1e23", 1e23)):
        i = np.where(ne >= thr)[0]
        out[lab] = float(t[i[0]]) if i.size else None
    pk = _rf_hist(run["rs"], s["f_drive_hz"])
    i = np.where(pk["pk"] < 0.5 * pk["pk"][0])[0]
    out["t_focus_field_half"] = float(pk["t_pk"][i[0]]) if i.size else None
    out["fil_ne_axis_max"] = s["peaks"]["fil_ne_axis_max"]
    # skip the first ns: the initial Te relaxes from its seed value and would dominate the peak
    late = t > 1e-9
    out["fil_te_axis_max_eV"] = float(ts["fil_te_axis_focus"][late].max()) if late.any() else s["peaks"]["fil_te_axis_max"]
    out["fil_tg_axis_max_K"] = s["peaks"]["fil_tg_axis_max"]
    out["bulk_ne_max_end"] = s["peaks"]["bulk_ne_max_end"]
    out["laser_delivered_J"] = s["laser"]["delivered_J"]
    out["laser_absorbed_J"] = s["laser"]["absorbed_J"]
    out["gamma_start"] = s["rf"].get("gamma_start")
    out["gamma_end"] = s["rf"].get("gamma_end")
    out["focus_peak_start_V_m"] = s["rf"].get("focus_peak_start_V_m")
    out["focus_peak_end_V_m"] = s["rf"].get("focus_peak_end_V_m")
    out["rf_balance_residual"] = s["rf"].get("balance_residual")
    out["filament_ledger_residual"] = s["ledger_filament_J"].get("closure_residual_rel")
    out["status"] = s["status"]
    return out


def geometry_data(geom, runs, title=None, field=8):
    base = runs["rfonly"]["fr"]
    lz = runs["L10"]["fr"]
    n = min(len(r["fr"]["t"]) for r in runs.values())
    stride = max(1, math.ceil(n / MAX_FRAMES))
    sel = np.arange(0, n, stride)
    dx = float(base["dx"]) * 1e3
    fi, fj, fk = (int(v) for v in base["focus"])
    ins = base["node_inside_xz"]
    nxn, nzn = ins.shape
    k0 = 16 if geom == "guide" else max(0, int(base["meas_k"]) - 4)
    k1 = nzn - 16 if geom == "guide" else nzn
    xs = slice(None, None, 2)
    zs = slice(k0, k1, 2)
    scale = float(max(np.abs(runs[k]["fr"]["ey_xz"][sel][:, :, k0:k1]).max() for k in ("rfonly", "L10")))
    e_f = _focus_edge(base)
    period_frames = max(1, int(round(1 / (runs["rfonly"]["s"]["f_drive_hz"] * (base["t"][1] - base["t"][0])))))
    blocked = ~ins[xs, zs]
    iris = [int(k) for k in range(k0, k1) if (~ins[:, k]).any() and ins[:, k].any()]
    d = dict(geom=geom, title=title or GEOMS[geom], field=field, t=np.round(base["t"][sel] * 1e9, 4).tolist(),
             nx=int(blocked.shape[0]), nz=int(blocked.shape[1]), dxs=2 * dx, z0=k0 * dx, z1=(k1 - 1) * dx,
             width=(nxn - 1) * dx, height=(base["ne_yz"].shape[1] - 1) * dx,
             focus=[fi * dx, (fj + 0.5) * dx, fk * dx], scale=scale / 1e6, ne_log=NE_LOG,
             blocked=_b64(blocked.astype(np.uint8)), iris=[k * dx for k in iris],
             meas=float(base["meas_k"]) * dx, trans=float(base["trans_k"]) * dx if int(base["trans_k"]) >= 0 else None,
             fil_y=np.round(base["fil_y"] * 1e3, 3).tolist(),
             beam_w_um=np.round(lz["beam_w"] * 1e6, 2).tolist(), runs={}, period_frames=period_frames,
             f_drive=float(runs["rfonly"]["s"]["f_drive_hz"]))
    d["eyz_ny"] = int(base["ey_yz"].shape[1])
    lut = plt.get_cmap("RdBu_r")(np.linspace(0, 1, 255))[:, :3]
    d["lut"] = np.rint(lut * 255).astype(int).tolist()
    d["lut_ne"] = np.rint(plt.get_cmap("inferno")(np.linspace(0, 1, 256))[:, :3] * 255).astype(int).tolist()
    for key, r in runs.items():
        fr, ts, s = r["fr"], r["ts"], r["s"]
        rh = _rf_hist(r["rs"], s["f_drive_hz"])
        blk = np.minimum(sel * int(s["case"]["frame_every_blocks"]), len(ts["t"]) - 1)
        lp = fr["laser_power"][:n, 0]
        lp_max = float(lp.max())
        env = _envelope(fr["ey_xz"][:n, fi, k0:k1], fr["t"][:n], s["f_drive_hz"], period_frames)[sel][:, ::2]
        rr = dict(label=r["label"],
                  fil_log=np.round(np.log10(np.maximum(fr["fil_ne_axis"][sel], 1.0)), 2).ravel().tolist(),
                  fil_te=np.round(fr["fil_te_axis"][sel], 2).ravel().tolist(),
                  r50_um=np.round(fr["fil_r50"][sel][:, e_f] * 1e6, 2).tolist(),
                  x2p=np.round(fr["fil_x2p_frac"][sel][:, e_f], 5).tolist(),
                  laser=np.round(lp[sel] / lp_max, 4).tolist() if lp_max else [0.0] * len(sel),
                  laser_peak_W=float(lp.max()),
                  em=np.round(np.asarray(fr["em_energy"])[sel] / max(float(fr["em_energy"][0]), 1e-30), 4).tolist(),
                  bulk_max=np.round(np.log10(np.maximum(ts["bulk_ne_max"][blk], 1.0)), 3).tolist(),
                  fil_focus=np.round(np.log10(np.maximum(ts["fil_ne_axis_focus"][blk], 1.0)), 3).tolist(),
                  te_focus=np.round(ts["fil_te_axis_focus"][blk], 3).tolist(),
                  tg_focus=np.round(ts["fil_tg_axis_focus"][blk], 1).tolist(),
                  env=_b64(np.rint(env / scale * 255).clip(0, 255).astype(np.uint8)),
                  rad_r=np.round(fr["fil_rb_focus"][sel] * 1e6, 3).ravel().tolist(),
                  rad_ne=np.round(np.log10(np.maximum(fr["fil_ne_r_focus"][sel], 1.0)), 3).ravel().tolist(),
                  rad_te=np.round(fr["fil_te_r_focus"][sel], 3).ravel().tolist(),
                  rad_tg=np.round(fr["fil_tg_r_focus"][sel], 1).ravel().tolist(),
                  rad_x2p=np.round(np.log10(np.maximum(fr["fil_x2p_r_focus"][sel], 1.0)), 3).ravel().tolist(),
                  rings=int(fr["fil_ne_r_focus"].shape[1]),
                  t_pk=np.round(rh["t_pk"] * 1e9, 4).tolist(),
                  pk=np.round(rh["pk"] / max(rh["pk"][0], 1e-30), 4).tolist(),
                  t_g=np.round(rh["t_g"] * 1e9, 4).tolist(), gamma=np.round(rh["gamma"], 4).tolist(),
                  trans=np.round(rh["trans"], 4).tolist(), milestones=_milestones(r))
        if key in ("rfonly", "L10"):
            rr["ey"] = _b64(_q_field(fr["ey_xz"][sel][:, xs, zs], scale))
            rr["ne"] = _b64(_q_log(fr["ne_xz"][sel][:, xs, zs]))
            rr["eyz"] = _b64(_q_field(fr["ey_yz"][sel][:, :, zs], scale))
        d["runs"][key] = rr
    return d


VIEWER = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Laser-triggered RF breakdown in 1 MPa helium · coupled 3D physics</title>
<style>
:root{color-scheme:dark;font-family:Segoe UI,Arial,sans-serif;background:#0a1422;color:#deebf6}
body{max-width:1500px;margin:0 auto;padding:20px}h1{font-size:25px;margin:4px 0 8px}p{line-height:1.5}
.sub{color:#aac4d9;margin-top:0}.notice{border-left:3px solid #ffc96c;padding:8px 14px;background:#282719;margin:8px 0}
.controls{display:flex;gap:12px;align-items:center;flex-wrap:wrap;background:#15263b;padding:10px;border-radius:10px;position:sticky;top:0;z-index:2}
button,select{background:#1b4964;color:white;border:1px solid #477890;border-radius:5px;padding:7px 12px;cursor:pointer}
input[type=range]{flex:1;min-width:220px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:12px}
.card{background:#101f31;border:1px solid #294158;border-radius:12px;padding:12px;min-width:0}
h2{font-size:16px;margin:2px 0 8px}canvas{width:100%;height:auto;display:block;border-radius:5px}
.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:8px;margin:8px 0}.metric{background:#172c42;padding:8px;border-radius:7px;font-size:13px}
.value{display:block;color:#86e8ea;font-size:18px;margin-top:3px}.small{font-size:13px;color:#b0c9db}
table{border-collapse:collapse;font-size:13px;width:100%}td,th{border-bottom:1px solid #294158;padding:4px 6px;text-align:right}th:first-child,td:first-child{text-align:left}
@media(max-width:900px){.grid{grid-template-columns:1fr}.metrics{grid-template-columns:1fr 1fr}}
@media(max-width:600px){body{padding:8px}h1{font-size:20px}.controls{gap:6px;padding:6px}.controls label{font-size:13px}
 select,button{padding:5px 8px;max-width:100%}input[type=range]{min-width:0;flex-basis:100%;order:5}#table{overflow-x:auto}
 .metric{font-size:11px}.value{font-size:14px}}
</style>
<!--DATA--></head><body>
<p id="loading" style="color:#ffd27a">Loading simulation data…</p>
<h1>Laser-triggered RF breakdown: how the plasma becomes a reflector</h1>
<p class="sub">Coupled 3D solve: Yee FDTD (WR-284, 3 GHz) ⇄ bulk electron kinetics ⇄ laser filament
(PPT photoionization, avalanche, electron energy, He⁺ → He₂⁺ conversion, recombination, radial hydrodynamics) · helium 1 MPa · RF 8 MV/m at focus · 532 nm, 20 µm waist ·
<a href="https://emaguidetal.github.io/plasma-viewer/scans.html" style="color:#7cc4ff"><b>Parameter scans and comparison with experiment →</b></a></p>
<div class="notice"><b>Nothing here is prescribed.</b> Densities, temperatures and fields are solver outputs. What the coupled runs show:
<ol style="margin:6px 0 2px 18px;padding:0">
<li><b>The laser alone does not make the plasma.</b> At 1 MPa helium, 532 nm ionization needs 11 photons. The 10 mJ / 3 ns pulse deposits only ~10⁻¹² J and leaves a micrometre-thin seed channel. In the 6 MV/m case (RF below self-breakdown), that channel stalls at ~3×10¹⁷ m⁻³ and |Γ| stays at 2×10⁻⁴: <i>no reflector</i>.</li>
<li><b>Laser + strong RF together make a dense channel early.</b> At 8 MV/m the RF-heated seed electrons absorb the laser and avalanche. The channel reaches 10²³ m⁻³ about 5 ns before the RF-only breakdown, peaks at ~6×10²⁴ m⁻³, Tgas ~2×10⁴ K, and locally shorts E<sub>y</sub> at the focus (see the envelope dip).</li>
<li><b>The reflector is the RF-driven volume breakdown, not the thin channel.</b> The channel is only 5–13 µm in radius, so its line resistance (~MΩ/m) is far too high to load a 3 GHz mode. The wave is turned around by the 10²⁰–10²¹ m⁻³ plasma that the RF avalanche builds over the gas volume. In the matched guide |Γ| rises from 0.002 to 0.50 and transmission falls to ~5%: ~25% of the power is reflected and ~70% absorbed. That is because at 1 MPa ν<sub>m</sub> ≈ 2×10¹³ s⁻¹ ≫ ω, so the plasma is a lossy conductor. With or without the laser, the global RF response (|Γ|, transmission, stored energy) is the same to within ~2%.</li></ol></div>
<div class="controls"><label>Case <select id="geom"><option value="cavity_E8">Iris cavity · 8 MV/m</option><option value="guide_E8">Matched guide · 8 MV/m</option><option value="guide_E6">Matched guide · 6 MV/m (sub-threshold)</option></select></label>
<button id="play">Pause</button><input id="frame" type="range" min="0" max="1" value="0" aria-label="time"><strong id="clock"></strong>
<label>Speed <select id="speed"><option value=".5">0.5×</option><option selected value="1">1×</option><option value="2">2×</option><option value="4">4×</option></select></label></div>
<div class="metrics">
<div class="metric">Laser (10 mJ / 3 ns)<span class="value" id="mLaser" style="color:#ff8080"></span></div>
<div class="metric">Filament axis n<sub>e</sub> at focus (laser / RF only)<span class="value" id="mNe"></span></div>
<div class="metric">Filament T<sub>e</sub>, T<sub>gas</sub>, He₂⁺ ion fraction (laser run)<span class="value" id="mTe"></span></div>
<div class="metric" id="mRfBox">RF at focus (laser / RF only)<span class="value" id="mRf"></span></div></div>
<div class="grid">
<section class="card"><h2 id="h0">RF only: E<sub>y</sub> (x–z mid-plane) and bulk n<sub>e</sub></h2><canvas id="c0" width="680" height="400"></canvas></section>
<section class="card"><h2 id="h1">RF + 10 mJ / 3 ns laser: E<sub>y</sub> and plasma</h2><canvas id="c1" width="680" height="400"></canvas></section>
<section class="card"><h2>3D: cut planes, laser beam (along −y, parallel to E<sub>y</sub>) and computed filament · drag to orbit</h2><canvas id="three" width="680" height="470" style="touch-action:none"></canvas>
<div class="small">This run uses the parallel geometry: the laser enters through the broad wall along E. The experiment uses a laser perpendicular to E (along x); see the <a href="https://emaguidetal.github.io/plasma-viewer/dynamics/orientation3d.html">3D orientation comparison</a>. x–z plane: E<sub>y</sub> of the laser run at y = mid-height. y–z plane: E<sub>y</sub> through the filament. The filament line is coloured by its computed axis n<sub>e</sub> at each height, and its width shows the computed half-electron radius. The red cone is the laser beam with its radius magnified ×40.</div></section>
<section class="card"><h2>Standing-wave envelope |E<sub>y</sub>| along the centreline (RF amplitude fitted over the last period)</h2><canvas id="env" width="680" height="300"></canvas>
<div class="small" id="envNote"></div>
<h2 style="margin-top:10px">Inside the filament at the focus: computed radial profiles (laser run)</h2><canvas id="rad" width="680" height="250"></canvas></section>
<section class="card" style="grid-column:1/-1"><h2>Time histories (all three runs)</h2><canvas id="hist" width="1400" height="330"></canvas></section>
</div>
<div class="card" style="margin-top:12px"><h2>Milestones from the solver</h2><div id="table"></div>
<h2 style="margin-top:12px">Physics included and not included</h2><ul class="small">
<li><b>Included and computed:</b> FDTD Maxwell solve with the plasma current fed back as a Drude term on every RF step. Bulk electron energy and ionization use BOLSIG-type Boltzmann swarm tables for helium. In the subgrid laser filament (48 radial rings per 1 mm cell along y), the solver computes Gaussian beam propagation and absorption, PPT multiphoton/tunnel ionization, inverse bremsstrahlung, RF Joule heating, avalanche, electron energy, He⁺ + 2He → He₂⁺ + He conversion, dissociative and three-body recombination, electron–gas heating and Lagrangian radial expansion.</li>
<li><b>Not in this coupled run:</b> multiple ionization (He²⁺), superelastic (excited-state → electron) energy return and resonance radiation transport. Excitation energy is a local sink here. These processes are implemented and tested separately in <code>helium_cr.py</code> (low-density demonstration only). Dense-helium He₂*/He₂⁺ data at 1 MPa are not available with traceable sources.</li>
<li>PEC copper walls, 1 mm mesh. The optical port is not meshed. The seed density is 10⁹ m⁻³ everywhere, roughly natural background ionization. 8 MV/m is just above the RF-only self-breakdown threshold of this model within 24 ns, and 6 MV/m is below it. A ring-to-ring diffusion term is included, but no ionization-wave or streamer model spreads the dense channel outward into the surrounding gas. A real channel could widen and reflect more than computed here.</li>
<li>The field animation is sampled every ~90 ps against a 333 ps RF period, so the wave phase advances visibly between frames.</li></ul></div>
<script>
let DS={},D,G='cavity_E8',index=0,playing=true,last=0,yaw=-.55,pitch=.6;
const dec8=s=>Int8Array.from(atob(s),c=>c.charCodeAt(0)),decU=s=>Uint8Array.from(atob(s),c=>c.charCodeAt(0));
function prep(d){if(!d||d._p)return d;d._p=1;d.B=decU(d.blocked);for(const k in d.runs){const r=d.runs[k];if(r.ey){r.EY=dec8(r.ey);r.NE=decU(r.ne);r.EYZ=dec8(r.eyz)}r.ENV=decU(r.env)}return d}
function lutc(q,a=1){const v=D.lut[Math.max(0,Math.min(254,q+127))];return `rgba(${v[0]},${v[1]},${v[2]},${a})`}
function nec(lg,a=1){const f=Math.max(0,Math.min(1,(lg-D.ne_log[0])/(D.ne_log[1]-D.ne_log[0])));const v=D.lut_ne[Math.round(f*255)];return `rgba(${v[0]},${v[1]},${v[2]},${a})`}
const off=document.createElement('canvas'),oc=off.getContext('2d');
function slice(cid,key){
 const r=D.runs[key],cv=document.getElementById(cid),ctx=cv.getContext('2d'),nx=D.nx,nz=D.nz;off.width=nz;off.height=nx;
 const im=oc.createImageData(nz,nx),o=index*nx*nz;
 for(let i=0;i<nx;i++)for(let k=0;k<nz;k++){const j=4*((nx-1-i)*nz+k),p=i*nz+k;let c;
  if(D.B[p]){c=[120,130,140]}else{const v=D.lut[Math.max(0,Math.min(254,r.EY[o+p]+127))];c=v.slice();const q=r.NE[o+p]/255;
   if(q>.32){const w=Math.min(1,(q-.32)/.45),n=D.lut_ne[Math.min(255,Math.round(110+145*w))];c=[c[0]*(1-w)+n[0]*w,c[1]*(1-w)+n[1]*w,c[2]*(1-w)+n[2]*w]}}
  im.data[j]=c[0];im.data[j+1]=c[1];im.data[j+2]=c[2];im.data[j+3]=255}
 oc.putImageData(im,0,0);ctx.fillStyle='#101f31';ctx.fillRect(0,0,680,400);ctx.imageSmoothingEnabled=false;
 const L=56,T=14,W=600,H=Math.min(300,W*D.width/(D.z1-D.z0)*1.6);ctx.drawImage(off,L,T,W,H);
 const zp=z=>L+(z-D.z0)/(D.z1-D.z0)*W,xp=x=>T+H-x/D.width*H;
 ctx.strokeStyle='#b9ccdb';ctx.strokeRect(L,T,W,H);
 const fz=zp(D.focus[2]),fx=xp(D.focus[0]),lg=r.fil_log[index*D.fil_y.length+fe()],
   rad=3+Math.max(0,lg-D.ne_log[0])/(D.ne_log[1]-D.ne_log[0])*12;
 ctx.beginPath();ctx.arc(fz,fx,rad,0,7);ctx.fillStyle=nec(lg,.95);ctx.fill();ctx.strokeStyle='#fff';ctx.lineWidth=1;ctx.stroke();
 if(r.laser[index]>.02){ctx.save();ctx.shadowColor='#ff2020';ctx.shadowBlur=20;ctx.strokeStyle=`rgba(255,60,60,${.3+.7*r.laser[index]})`;ctx.lineWidth=3;ctx.beginPath();ctx.arc(fz,fx,rad+6,0,7);ctx.stroke();ctx.restore();
  ctx.fillStyle='#ff9a9a';ctx.font='bold 13px Segoe UI';ctx.fillText('⊙ laser (into the page, −y)',fz+rad+10,fx-8)}
 ctx.font='13px Segoe UI';ctx.fillStyle='#d8e8f4';
 for(let z=Math.ceil(D.z0/20)*20;z<=D.z1;z+=20)ctx.fillText(String(z),zp(z)-9,T+H+16);
 ctx.fillText('z (mm) → from the RF port',L+W/2-80,T+H+33);ctx.fillText('x',18,T+H/2);ctx.fillText(D.width.toFixed(0),24,T+10);ctx.fillText('0',38,T+H);
 ctx.fillStyle='#9fb3c4';ctx.fillText(D.geom==='guide'?'matched load (CPML) →':'PEC short',L+W-(D.geom==='guide'?140:62),T-2+H+48);
 for(let j=0;j<160;j++){ctx.fillStyle=lutc(Math.round((j/159*2-1)*127));ctx.fillRect(L+j,T+H+56,1,9)}
 ctx.fillStyle='#d8e8f4';ctx.font='12px Segoe UI';ctx.fillText(`Ey ±${D.scale.toFixed(1)} MV/m`,L+165,T+H+65);
 for(let j=0;j<160;j++){ctx.fillStyle=nec(D.ne_log[0]+(j/159)*(D.ne_log[1]-D.ne_log[0]));ctx.fillRect(L+300+j,T+H+56,1,9)}
 ctx.fillText(`ne 10^${D.ne_log[0]}–10^${D.ne_log[1]} m⁻³ (bulk overlay; dot = filament axis)`,L+300,T+H+78);
}
function fe(){let b=0,e=1e9;D.fil_y.forEach((y,i)=>{const d=Math.abs(y-D.focus[1]);if(d<e){e=d;b=i}});return b}
function project(z,x,y){const a=z-(D.focus[2]),b=x-D.width/2,h=y-D.height/2;const u=Math.cos(yaw)*a-Math.sin(yaw)*b,v=Math.sin(yaw)*a+Math.cos(yaw)*b;const s=560/(D.z1-D.z0+40);return [340+s*u,250+s*(Math.sin(pitch)*v-Math.cos(pitch)*h),Math.cos(pitch)*v+Math.sin(pitch)*h]}
function poly(ctx,pts,fill,stroke){ctx.beginPath();pts.forEach((v,i)=>{const p=project(...v);i?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1])});ctx.closePath();if(fill){ctx.fillStyle=fill;ctx.fill()}if(stroke){ctx.strokeStyle=stroke;ctx.stroke()}}
function three(){
 const r=D.runs.L10,ctx=document.getElementById('three').getContext('2d');ctx.fillStyle='#101f31';ctx.fillRect(0,0,680,470);
 const faces=[],st=2,nx=D.nx,nz=D.nz,ny=D.eyz_ny,o=index*nx*nz,oy=index*ny*nz,dz=D.dxs,a=D.width,b=D.height;
 for(let i=0;i<nx-1;i+=st)for(let k=0;k<nz-1;k+=st){const z=D.z0+k*dz,x=i*dz;const p=i*nz+k;if(D.B[p])continue;
  faces.push({pts:[[z,x,b/2],[z+st*dz,x,b/2],[z+st*dz,x+st*dz,b/2],[z,x+st*dz,b/2]],col:lutc(r.EY[o+p],.85),d:project(z,x,b/2)[2]})}
 for(let j=0;j<ny-1;j+=st)for(let k=0;k<nz-1;k+=st){const z=D.z0+k*dz,y=(j+.5)*dz/2;
  faces.push({pts:[[z,a/2,y],[z+st*dz,a/2,y],[z+st*dz,a/2,y+st*dz/2],[z,a/2,y+st*dz/2]],col:lutc(r.EYZ[oy+j*nz+k],.8),d:project(z,a/2,y)[2]})}
 faces.sort((p,q)=>p.d-q.d);for(const f of faces)poly(ctx,f.pts,f.col,null);
 ctx.lineWidth=1.2;for(const z of [D.z0,D.z1])poly(ctx,[[z,0,0],[z,a,0],[z,a,b],[z,0,b]],null,'#8baabc');
 for(const x of [0,a])for(const y of [0,b]){const p=project(D.z0,x,y),q=project(D.z1,x,y);ctx.strokeStyle='#8baabc';ctx.beginPath();ctx.moveTo(p[0],p[1]);ctx.lineTo(q[0],q[1]);ctx.stroke()}
 for(const zi of D.iris)poly(ctx,[[zi,0,0],[zi,a,0],[zi,a,b],[zi,0,b]],null,'rgba(160,180,195,.25)');
 const [fx,fy,fz]=D.focus,ne=D.fil_y.length;
 if(r.laser[index]>.01){ctx.save();ctx.shadowColor='#ff2020';ctx.shadowBlur=16;for(let e=0;e<ne;e++){const w=D.beam_w_um[e]*40/1000,y=D.fil_y[e];
  poly(ctx,[[fz-w,fx,y-.5],[fz+w,fx,y-.5],[fz+w,fx,y+.5],[fz-w,fx,y+.5]],`rgba(255,50,50,${.15+.6*r.laser[index]})`,null)}
  const p=project(fz,fx,b+14),q=project(fz,fx,b);ctx.strokeStyle='#ff6060';ctx.lineWidth=3;ctx.beginPath();ctx.moveTo(p[0],p[1]);ctx.lineTo(q[0],q[1]);ctx.stroke();ctx.fillStyle='#ff9a9a';ctx.font='bold 13px Segoe UI';ctx.fillText('laser ↓ (−y)',p[0]+6,p[1]);ctx.restore()}
 for(let e=0;e<ne;e++){const lg=r.fil_log[index*ne+e],y=D.fil_y[e],p=project(fz,fx,y-.5),q=project(fz,fx,y+.5);
  ctx.strokeStyle=nec(lg);ctx.lineWidth=2+Math.max(0,(lg-D.ne_log[0]))/(D.ne_log[1]-D.ne_log[0])*9;ctx.beginPath();ctx.moveTo(p[0],p[1]);ctx.lineTo(q[0],q[1]);ctx.stroke()}
 ctx.fillStyle='#d7e5ef';ctx.font='13px Segoe UI';ctx.fillText('3D solution · laser run · '+D.title,14,22);
 const pp=project(D.z0,0,0);ctx.fillText('RF port',pp[0]-10,pp[1]+16);const ps=project(D.z1,a,0);ctx.fillText(D.geom==='guide'?'matched load':'short',ps[0]+4,ps[1]+14);
}
function envelope(){
 const c=document.getElementById('env'),ctx=c.getContext('2d');ctx.fillStyle='#101f31';ctx.fillRect(0,0,680,300);const L=50,T=16,W=610,H=230;
 ctx.strokeStyle='#36556e';ctx.strokeRect(L,T,W,H);const zp=z=>L+(z-D.z0)/(D.z1-D.z0)*W,yp=v=>T+H-v/255*H;
 for(const zi of D.iris){ctx.fillStyle='rgba(160,180,195,.25)';ctx.fillRect(zp(zi)-1,T,2,H)}
 ctx.strokeStyle='rgba(255,90,90,.6)';ctx.setLineDash([4,4]);ctx.beginPath();ctx.moveTo(zp(D.focus[2]),T);ctx.lineTo(zp(D.focus[2]),T+H);ctx.stroke();ctx.setLineDash([]);
 for(const [k,col] of [['rfonly','#efbc62'],['L1','#9be37a'],['L10','#49e1e9']]){const r=D.runs[k];if(!r)continue;ctx.strokeStyle=col;ctx.lineWidth=k==='L1'?1:2;ctx.beginPath();
  for(let q=0;q<D.nz;q++){const x=zp(D.z0+q*D.dxs),y=yp(r.ENV[index*D.nz+q]);q?ctx.lineTo(x,y):ctx.moveTo(x,y)}ctx.stroke()}
 ctx.font='12px Segoe UI';ctx.fillStyle='#c9dce9';ctx.fillText(`0`,L-12,T+H);ctx.fillText(`${D.scale.toFixed(1)} MV/m`,4,T+10);
 for(let z=Math.ceil(D.z0/20)*20;z<=D.z1;z+=20)ctx.fillText(String(z),zp(z)-9,T+H+15);ctx.fillText('z (mm)',L+W-40,T+H+30);
 ctx.fillStyle='#efbc62';ctx.fillText('RF only',L+8,T+H+30);ctx.fillStyle='#9be37a';ctx.fillText('1 mJ',L+70,T+H+30);ctx.fillStyle='#49e1e9';ctx.fillText('10 mJ',L+115,T+H+30);ctx.fillStyle='#ff8080';ctx.fillText('laser focus',zp(D.focus[2])+4,T+12);
}
function radial(){
 const r=D.runs.L10,c=document.getElementById('rad'),ctx=c.getContext('2d');ctx.fillStyle='#101f31';ctx.fillRect(0,0,680,250);
 const n=r.rings,o=index*n,ob=index*(n+1),L=50,T=12,W=600,H=190;ctx.strokeStyle='#36556e';ctx.strokeRect(L,T,W,H);
 const rmax=Math.max(...r.rad_r.slice(ob,ob+n+1)),rc=i=>.5*(r.rad_r[ob+i]+r.rad_r[ob+i+1]);
 const xr=v=>L+Math.sqrt(v/rmax)*W;
 const curves=[[r.rad_ne,'#ffb347',v=>(v-9)/(26-9),'ne log10 9–26'],[r.rad_x2p,'#d68bff',v=>(v-9)/(26-9),'He₂⁺ log10 9–26'],
  [r.rad_te,'#49e1e9',v=>v/Math.max(10,Math.max(...r.rad_te)),'Te (eV, 0–max)'],[r.rad_tg,'#ff6b6b',v=>(v-300)/Math.max(100,Math.max(...r.rad_tg)-300),'Tgas (K, 300–max)']];
 for(const [arr,col,f] of curves){ctx.strokeStyle=col;ctx.lineWidth=2;ctx.beginPath();for(let i=0;i<n;i++){const x=xr(rc(i)),y=T+H-Math.max(0,Math.min(1,f(arr[o+i])))*H;i?ctx.lineTo(x,y):ctx.moveTo(x,y)}ctx.stroke()}
 ctx.font='12px Segoe UI';let x=L;for(const [,col,,lab] of curves){ctx.fillStyle=col;ctx.fillText(lab,x,T+H+32);x+=ctx.measureText(lab).width+18}
 ctx.fillStyle='#c9dce9';for(const f of [0,.25,.5,1]){const v=f*f*rmax;ctx.fillText(v<1000?v.toFixed(0)+' µm':(v/1000).toFixed(2)+' mm',xr(v)-14,T+H+15)}
 ctx.fillText(`max Te ${Math.max(...r.rad_te.slice(o,o+n)).toFixed(2)} eV · max Tgas ${Math.max(...r.rad_tg.slice(o,o+n)).toFixed(0)} K · outer radius (√r axis)`,L+250,T+H+32-16);
}
function lineplot(ctx,x0,y0,w,h,t,series,ylab,ymin,ymax,log){
 ctx.strokeStyle='#36556e';ctx.strokeRect(x0,y0,w,h);const tm=D.t.at(-1),yp=v=>y0+h-(Math.max(ymin,Math.min(ymax,v))-ymin)/(ymax-ymin)*h;
 for(const [tt,v,col,lw] of series){if(!v)continue;ctx.strokeStyle=col;ctx.lineWidth=lw||1.6;ctx.beginPath();v.forEach((q,i)=>{const x=x0+tt[i]/tm*w,y=yp(q);i?ctx.lineTo(x,y):ctx.moveTo(x,y)});ctx.stroke()}
 const xc=x0+D.t[index]/tm*w;ctx.strokeStyle='#fff';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(xc,y0);ctx.lineTo(xc,y0+h);ctx.stroke();
 ctx.fillStyle='#c9dce9';ctx.font='12px Segoe UI';ctx.fillText(ylab,x0,y0-4);ctx.fillText(log?'1e'+ymax:String(ymax),x0+3,y0+12);ctx.fillText(log?'1e'+ymin:String(ymin),x0+3,y0+h-3);ctx.fillText('t (ns) → '+tm.toFixed(0),x0+w-70,y0+h+14);
}
function hist(){
 const ctx=document.getElementById('hist').getContext('2d');ctx.fillStyle='#101f31';ctx.fillRect(0,0,1400,330);const R=D.runs,t=D.t,cols={rfonly:'#efbc62',L1:'#9be37a',L10:'#49e1e9'};
 const s1=[];for(const k of ['rfonly','L1','L10'])if(R[k]){s1.push([t,R[k].fil_focus,cols[k],2]);if(k==='rfonly')s1.push([t,R[k].bulk_max,'#a08050',1])}
 s1.push([t,D.runs.L10.laser.map(v=>9+v*16),'#ff5a5a',1.4]);
 lineplot(ctx,30,24,420,250,t,s1,'filament axis ne at focus (m⁻³, log) · red: laser power · brown: RF-only bulk max',9,25,true);
 if(D.geom==='guide'){lineplot(ctx,490,24,420,250,null,[['rfonly','#efbc62'],['L1','#9be37a'],['L10','#49e1e9']].filter(([k])=>R[k]).flatMap(([k,c])=>[[R[k].t_g,R[k].gamma,c,2],[R[k].t_g,R[k].trans,c,1]]),'|Γ| at port (thick) · transmitted power / initial (thin)',0,1.2)}
 else{lineplot(ctx,490,24,420,250,null,[['rfonly','#efbc62'],['L1','#9be37a'],['L10','#49e1e9']].filter(([k])=>R[k]).flatMap(([k,c])=>[[R[k].t_pk,R[k].pk,c,2],[t,R[k].em,c,1]]),'focus |Ey| peak / initial (thick) · stored EM energy / initial (thin)',0,1.2)}
 lineplot(ctx,950,24,420,250,null,[[t,R.L10.te_focus.map(v=>v/10),'#49e1e9',2],[t,R.rfonly.te_focus.map(v=>v/10),'#efbc62',1.4],[t,R.L10.tg_focus.map(v=>v/1e4),'#ff6b6b',2],[t,R.L10.x2p.map(v=>Math.max(0,Math.log10(Math.max(v,1e-12))+12)/12),'#d68bff',2]],'Te/10 eV (cyan, gold) · Tgas/10⁴ K (red) · log He₂⁺ ion frac 1e-12→1 (violet)',0,1.2);
 ctx.font='13px Segoe UI';let x=30;for(const [k,l] of [['rfonly','RF only'],['L1','RF + 1 mJ/500 ps'],['L10','RF + 10 mJ/3 ns']]){ctx.fillStyle=cols[k];ctx.fillText('■ '+l,x,318);x+=180}
}
const fmt=v=>v==null?'—':Math.abs(v)>=1e4||Math.abs(v)<1e-2?v.toExponential(2):v.toPrecision(3);
function table(){const R=D.runs,rows=[['first filament ne ≥ 1e17 m⁻³ (ns)','t_fil_1e17',1e9],['first filament ne ≥ 1e20 (ns)','t_fil_1e20',1e9],['first filament ne ≥ 1e23 (ns)','t_fil_1e23',1e9],['focus field < 50 % (ns)','t_focus_field_half',1e9],['peak filament axis ne (m⁻³)','fil_ne_axis_max',1],['peak filament Te after 1 ns (eV)','fil_te_axis_max_eV',1],['peak filament Tgas (K)','fil_tg_axis_max_K',1],['bulk ne max at end (m⁻³)','bulk_ne_max_end',1],['laser absorbed (J)','laser_absorbed_J',1],['|Γ| start → end','gamma',1],['RF energy-balance residual','rf_balance_residual',1],['filament energy-ledger residual','filament_ledger_residual',1]];
 let h='<table><tr><th>quantity</th>'+Object.keys(R).map(k=>`<th>${R[k].label}</th>`).join('')+'</tr>';
 for(const [lab,key,s] of rows){h+=`<tr><td>${lab}</td>`;for(const k in R){const m=R[k].milestones;h+='<td>'+(key==='gamma'?`${fmt(m.gamma_start)} → ${fmt(m.gamma_end)}`:fmt(m[key]==null?null:m[key]*s))+'</td>'}h+='</tr>'}
 document.getElementById('table').innerHTML=h+'</table>'}
function metrics(){const n=D.fil_y.length,e=fe(),L=D.runs.L10,R0=D.runs.rfonly;
 document.getElementById('mLaser').textContent=L.laser[index]>.01?`on · ${(L.laser[index]*L.laser_peak_W/1e6).toFixed(2)} MW`:(D.t[index]<5?'armed (fires 5 ns)':'off');
 document.getElementById('mNe').textContent=`${Math.pow(10,L.fil_focus[index]).toExponential(1)} / ${Math.pow(10,R0.fil_focus[index]).toExponential(1)} m⁻³`;
 document.getElementById('mTe').textContent=`${L.te_focus[index].toFixed(2)} eV · ${L.tg_focus[index].toFixed(0)} K · ${L.x2p[index].toExponential(1)}`;
 const pk=(r)=>{let i=r.t_pk.findIndex(v=>v>=D.t[index]);if(i<0)i=r.t_pk.length-1;return r.pk[i]},gm=(r)=>{let i=r.t_g.findIndex(v=>v>=D.t[index]);if(i<0)i=r.t_g.length-1;return r.gamma[i]};
 document.getElementById('mRf').textContent=D.geom==='guide'?`|Γ| ${gm(L).toFixed(2)} / ${gm(R0).toFixed(2)}`:`${(100*pk(L)).toFixed(0)} % / ${(100*pk(R0)).toFixed(0)} %`;
 document.getElementById('envNote').textContent=D.geom==='guide'?'In a matched guide the envelope starts flat (a traveling wave). As the RF avalanche fills the gas volume (about 18 ns on), the field drops along the plasma and a reflected standing-wave component builds up upstream. The laser channel only adds a narrow local dip at the focus; it does not change the global response.':'The cavity already contains a standing wave. The RF-driven volume plasma loads the cavity and the stored field collapses. The laser channel at the field maximum shorts Ey only locally (the narrow dip), so the global decay is set by the bulk plasma.'}
function draw(){const sl=document.getElementById('frame');sl.max=D.t.length-1;sl.value=index;document.getElementById('clock').textContent=D.t[index].toFixed(2)+' ns';
 slice('c0','rfonly');slice('c1','L10');three();envelope();radial();hist();metrics()}
function setGeom(g){G=g;D=prep(DS[g]);index=Math.min(index,D.t.length-1);table();draw()}
document.getElementById('geom').onchange=e=>setGeom(e.target.value);
document.getElementById('play').onclick=()=>{playing=!playing;document.getElementById('play').textContent=playing?'Pause':'Play'};
document.getElementById('frame').oninput=e=>{index=+e.target.value;playing=false;document.getElementById('play').textContent='Play';draw()};
let drag=null;const tc=document.getElementById('three');tc.onpointerdown=e=>{drag=[e.clientX,e.clientY];tc.setPointerCapture(e.pointerId)};tc.onpointerup=()=>drag=null;
tc.onpointermove=e=>{if(drag){yaw+=(e.clientX-drag[0])*.008;pitch=Math.max(.1,Math.min(1.45,pitch+(e.clientY-drag[1])*.008));drag=[e.clientX,e.clientY];three()}};
function tick(now){if(playing&&now-last>1000/(12*+document.getElementById('speed').value)){index=(index+1)%D.t.length;draw();last=now}requestAnimationFrame(tick)}
async function gunzip(b64){const b=Uint8Array.from(atob(b64),c=>c.charCodeAt(0));
 return JSON.parse(await new Response(new Blob([b]).stream().pipeThrough(new DecompressionStream('gzip'))).text())}
async function boot(){
 if(window.DATA_GZ){for(const k in DATA_GZ)DS[k]=await gunzip(DATA_GZ[k])}
 else DS={cavity_E8:window.DATA_cavity_E8,guide_E8:window.DATA_guide_E8,guide_E6:window.DATA_guide_E6};
 for(const o of document.getElementById('geom').options)if(!DS[o.value])o.disabled=true;if(!DS[G]){G=Object.keys(DS).find(k=>DS[k]);document.getElementById('geom').value=G}
 document.getElementById('loading').remove();setGeom(G);requestAnimationFrame(tick)}
boot().catch(e=>{document.getElementById('loading').textContent='Could not load the simulation data ('+e+'). Open this file in a current Chrome, Edge, Firefox or Safari browser.'});
window.viewerState=()=>({geom:G,frame:index,frames:D.t.length,time_ns:D.t[index]});
</script></body></html>"""


def overview(out, key, runs, d):
    geom = d["geom"]
    fr0, fr1 = runs["rfonly"]["fr"], runs["L10"]["fr"]
    t = fr0["t"]
    ins = fr0["node_inside_xz"]
    k0 = 16 if geom == "guide" else max(0, int(fr0["meas_k"]) - 4)
    k1 = ins.shape[1] - 16 if geom == "guide" else ins.shape[1]
    dx = float(fr0["dx"]) * 1e3
    ext = [k0 * dx, (k1 - 1) * dx, 0, (ins.shape[0] - 1) * dx]
    times = [5e-9, 10e-9, 16e-9, 20e-9, 23.5e-9]
    sc = d["scale"] * 1e6
    fig, ax = plt.subplots(3, len(times), figsize=(19, 8.6), constrained_layout=True)
    for c, tt in enumerate(times):
        m = int(np.argmin(np.abs(t - tt)))
        for r, (fr, lab) in enumerate(((fr0, "RF only"), (fr1, "RF + 10 mJ laser"))):
            a = ax[r, c]
            im = a.imshow(fr["ey_xz"][m][:, k0:k1] / 1e6, origin="lower", extent=ext, cmap="RdBu_r",
                          vmin=-sc / 1e6, vmax=sc / 1e6, aspect="auto")
            ne = fr["ne_xz"][m][:, k0:k1]
            if ne.max() > 1e16:
                a.contour(np.linspace(ext[0], ext[1], ne.shape[1]), np.linspace(0, ext[3], ne.shape[0]),
                          np.log10(np.maximum(ne, 1)), levels=[16, 18, 20], colors=["#ffe066", "#ff9f1c", "#ffffff"],
                          linewidths=0.8)
            e = _focus_edge(fr)
            lg = math.log10(max(fr["fil_ne_axis"][m][e], 1))
            a.scatter([float(fr["focus"][2]) * dx], [float(fr["focus"][0]) * dx], s=20 + 12 * max(0, lg - 9),
                      c=[lg], cmap="inferno", vmin=9, vmax=25, edgecolors="w")
            a.set_title(f"{lab} · t = {t[m] * 1e9:.1f} ns\nfilament ne = 1e{lg:.1f} m⁻³", fontsize=9)
            if c == 0:
                a.set_ylabel("x (mm)")
        a = ax[2, c]
        m = int(np.argmin(np.abs(t - tt)))
        for rk, col in (("rfonly", "#c98a00"), ("L1", "#4caf50"), ("L10", "#0097a7")):
            if rk in runs:
                fr = runs[rk]["fr"]
                env = _envelope(fr["ey_xz"][:, int(fr["focus"][0]), k0:k1], fr["t"], d["f_drive"], d["period_frames"])
                a.plot(np.arange(k0, k1) * dx, env[m] / 1e6, color=col, label=runs[rk]["label"])
        a.axvline(float(fr0["focus"][2]) * dx, color="r", ls="--", lw=0.8)
        a.set_xlabel("z (mm)")
        if c == 0:
            a.set_ylabel("|Ey| envelope (MV/m)")
            a.legend(fontsize=7)
    fig.colorbar(im, ax=ax[:2, :], shrink=0.6, label="Ey (MV/m)")
    fig.suptitle(f"Coupled 3D solve, {d['title']}, He 1 MPa: "
                 "contours show bulk ne = 1e16/1e18/1e20 m⁻³; dot = computed filament axis ne", fontsize=11)
    fig.savefig(out / f"overview_{key}.png", dpi=110)
    plt.close(fig)


def histories_png(out, data):
    keys = [k for k in DATASETS if k in data]
    fig, ax = plt.subplots(len(keys), 3, figsize=(17, 4.3 * len(keys)), constrained_layout=True, squeeze=False)
    cols = {"rfonly": "#c98a00", "L1": "#4caf50", "L10": "#0097a7"}
    for row, key in enumerate(keys):
        d = data[key]
        geom = d["geom"]
        t = np.asarray(d["t"])
        for k, r in d["runs"].items():
            ax[row, 0].semilogy(t, 10 ** np.asarray(r["fil_focus"]), color=cols[k], lw=2, label=f"{r['label']}: filament axis")
            if k == "rfonly":
                ax[row, 0].semilogy(t, 10 ** np.asarray(r["bulk_max"]), color=cols[k], lw=1, ls=":", label="RF only: bulk max")
            if geom == "guide":
                ax[row, 1].plot(r["t_g"], r["gamma"], color=cols[k], lw=2, label=f"{r['label']} |Γ|")
                ax[row, 1].plot(r["t_g"], r["trans"], color=cols[k], lw=1, ls="--")
            else:
                ax[row, 1].plot(r["t_pk"], r["pk"], color=cols[k], lw=2, label=f"{r['label']}")
            ax[row, 2].plot(t, r["te_focus"], color=cols[k], lw=1.6, label=f"{r['label']} Te (eV)")
        lz = d["runs"]["L10"]
        ax2 = ax[row, 0].twinx()
        ax2.plot(t, lz["laser"], color="r", lw=1)
        ax2.set_ylabel("laser power / peak", color="r")
        ax[row, 0].set_ylim(1e8, 1e26)
        ax[row, 0].set_title(f"{d['title']}\nelectron density at the focus", fontsize=9)
        ax[row, 1].set_title(f"{geom}: " + ("|Γ| at port (solid), transmitted/initial (dashed)" if geom == "guide"
                                            else "focus |Ey| peak / initial"))
        ax[row, 1].set_ylim(0, 1.25)
        ax3 = ax[row, 2].twinx()
        ax3.semilogy(t, np.maximum(lz["x2p"], 1e-14), color="purple", lw=1.5)
        ax3.set_ylabel("He₂⁺ molecular-ion fraction at focus (10 mJ)", color="purple")
        ax[row, 2].set_title(f"{geom}: filament Te (right: 10 mJ He₂⁺ ion fraction)")
        for a in ax[row]:
            a.set_xlabel("t (ns)")
            a.grid(alpha=0.3)
        ax[row, 0].legend(fontsize=7, loc="upper left")
        ax[row, 1].legend(fontsize=7)
        ax[row, 2].legend(fontsize=7)
    fig.savefig(out / "histories.png", dpi=110)
    plt.close(fig)


def gif(out, key, runs, d, n_frames=120):
    geom = d["geom"]
    fr0, fr1 = runs["rfonly"]["fr"], runs["L10"]["fr"]
    t = fr0["t"]
    ins = fr0["node_inside_xz"]
    k0 = 16 if geom == "guide" else max(0, int(fr0["meas_k"]) - 4)
    k1 = ins.shape[1] - 16 if geom == "guide" else ins.shape[1]
    dx = float(fr0["dx"]) * 1e3
    ext = [k0 * dx, (k1 - 1) * dx, 0, (ins.shape[0] - 1) * dx]
    sc = d["scale"]
    idx = np.linspace(0, len(t) - 1, n_frames).astype(int)
    fig, ax = plt.subplots(2, 2, figsize=(12.5, 7.2), dpi=64, constrained_layout=True)
    ims, dots, cont = [], [], []
    ne_cmap = colors.ListedColormap(plt.get_cmap("inferno")(np.linspace(0.35, 1, 128)))
    for a, fr, lab in ((ax[0, 0], fr0, "RF only"), (ax[0, 1], fr1, "RF + 10 mJ / 3 ns laser")):
        im = a.imshow(fr["ey_xz"][0][:, k0:k1] / 1e6, origin="lower", extent=ext, cmap="RdBu_r", vmin=-sc, vmax=sc,
                      aspect="auto")
        ov = a.imshow(np.ma.masked_less(np.log10(np.maximum(fr["ne_xz"][0][:, k0:k1], 1)), 14), origin="lower",
                      extent=ext, cmap=ne_cmap, vmin=14, vmax=22, aspect="auto", alpha=0.85)
        dot = a.scatter([float(fr["focus"][2]) * dx], [float(fr["focus"][0]) * dx], s=30, c=[9], cmap="inferno",
                        vmin=9, vmax=25, edgecolors="w", zorder=5)
        a.set_title(lab)
        a.set_xlabel("z (mm)")
        a.set_ylabel("x (mm)")
        ims.append((im, ov, fr))
        dots.append(dot)
    envs = {}
    for rk, col in (("rfonly", "#c98a00"), ("L1", "#4caf50"), ("L10", "#0097a7")):
        if rk in runs:
            fr = runs[rk]["fr"]
            envs[rk] = (_envelope(fr["ey_xz"][:, int(fr["focus"][0]), k0:k1], fr["t"], d["f_drive"], d["period_frames"]) / 1e6,
                        ax[1, 0].plot([], [], color=col, label=runs[rk]["label"])[0])
    zz = np.arange(k0, k1) * dx
    ax[1, 0].set_xlim(zz[0], zz[-1])
    ax[1, 0].set_ylim(0, sc * 1.05)
    ax[1, 0].axvline(float(fr0["focus"][2]) * dx, color="r", ls="--", lw=0.8)
    ax[1, 0].set_xlabel("z (mm)")
    ax[1, 0].set_ylabel("|Ey| envelope (MV/m)")
    ax[1, 0].legend(fontsize=8, loc="upper left")
    ax[1, 0].set_title("standing-wave envelope on the centreline")
    tt = d["t"]
    a = ax[1, 1]
    for rk, col in (("rfonly", "#c98a00"), ("L1", "#4caf50"), ("L10", "#0097a7")):
        if rk in d["runs"]:
            a.semilogy(tt, 10 ** np.asarray(d["runs"][rk]["fil_focus"]), color=col, lw=2)
    a.semilogy(tt, 10 ** (9 + 16 * np.asarray(d["runs"]["L10"]["laser"])), color="r", lw=1, label="laser power (scaled)")
    a.set_ylim(1e8, 1e26)
    a.set_xlabel("t (ns)")
    a.set_ylabel("filament axis ne at focus (m⁻³)")
    a.legend(fontsize=8, loc="upper left")
    cur = a.axvline(0, color="k")
    title = fig.suptitle("")

    def upd(n):
        m = idx[n]
        for (im, ov, fr), dot in zip(ims, dots):
            im.set_data(fr["ey_xz"][m][:, k0:k1] / 1e6)
            ov.set_data(np.ma.masked_less(np.log10(np.maximum(fr["ne_xz"][m][:, k0:k1], 1)), 14))
            lg = math.log10(max(fr["fil_ne_axis"][m][_focus_edge(fr)], 1))
            dot.set_array(np.array([lg]))
            dot.set_sizes([30 + 14 * max(0, lg - 9)])
        for env, ln in envs.values():
            ln.set_data(zz, env[m])
        cur.set_xdata([t[m] * 1e9])
        lp = runs["L10"]["fr"]["laser_power"][m, 0]
        title.set_text(f"{d['title']} · He 1 MPa · t = {t[m] * 1e9:5.2f} ns"
                       + ("   LASER ON" if lp > 0.02 * runs['L10']['fr']['laser_power'][:, 0].max() else ""))
        return []

    anim = animation.FuncAnimation(fig, upd, frames=len(idx), blit=False)
    anim.save(out / f"laser_trigger_{key}.gif", writer=animation.PillowWriter(fps=10))
    plt.close(fig)


def make(root=ROOT, out=None, gifs=True):
    out = out or root / "viewer"
    out.mkdir(parents=True, exist_ok=True)
    data, report = {}, {}
    for key, (geom, field, title) in DATASETS.items():
        runs = load_geometry(geom, root, field)
        if not {"rfonly", "L10"} <= runs.keys():
            continue
        d = geometry_data(geom, runs, title, field)
        data[key] = d
        (out / f"data_{key}.js").write_text(f"window.DATA_{key}=" + json.dumps(d, separators=(",", ":")) + ";",
                                            encoding="utf-8")
        overview(out, key, runs, d)
        if gifs:
            gif(out, key, runs, d)
        report[key] = {k: r["milestones"] for k, r in d["runs"].items()}
    tags = "".join(f'<script src="data_{k}.js"></script>' for k in data)
    (out / "laser_trigger.html").write_text(VIEWER.replace("<!--DATA-->", tags), encoding="utf-8")
    gz = {k: base64.b64encode(gzip.compress(json.dumps(d, separators=(",", ":")).encode(), 9, mtime=0)).decode()
          for k, d in data.items()}
    (out / "laser_trigger_standalone.html").write_text(
        VIEWER.replace("<!--DATA-->", "<script>window.DATA_GZ=" + json.dumps(gz) + ";</script>"), encoding="utf-8")
    histories_png(out, data)
    (out / "milestones.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m rf_plasma.laser_trigger_visuals")
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--no-gif", action="store_true")
    a = ap.parse_args(argv)
    print(json.dumps(make(a.root, gifs=not a.no_gif), indent=2, default=float))


if __name__ == "__main__":
    main()
