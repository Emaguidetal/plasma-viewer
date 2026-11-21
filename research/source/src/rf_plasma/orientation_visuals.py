"""3D viewer comparing a laser perpendicular and parallel to the WR-284 RF field.

Reads matched `cavity_dynamics` runs (same field, mesh, seed kernel and coupling
block; only the laser axis differs) and writes an orbitable 3D page, a static
summary figure and a metrics JSON.
"""
from __future__ import annotations

import base64
import gzip
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
GEOMETRIES = (("perp", "Perpendicular: laser along x, across E (experiment)"),
              ("par", "Parallel: laser along y, through the broad wall, along E"))
NE_FLOOR = 8.0


def _load(run):
    with np.load(run / "frames.npz") as f:
        frames = {k: f[k] for k in f.files}
    return frames, json.loads((run / "summary.json").read_text())


def _b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def _q_signed(a, scale):
    return np.clip(np.rint(127 * np.asarray(a) / scale), -127, 127).astype(np.int8)


def _q_log(a, hi):
    lg = np.log10(np.maximum(np.asarray(a, dtype=float), 10 ** NE_FLOOR))
    return np.clip(np.rint(255 * (lg - NE_FLOOR) / (hi - NE_FLOOR)), 0, 255).astype(np.uint8)


def _first_time(t, y, level):
    hit = np.nonzero(np.asarray(y) >= level)[0]
    return float(t[hit[0]]) if hit.size else None


def metrics(frames, summary):
    t = frames["time_s"]
    line = frames["ne_line"]
    peak = line.max(axis=1)
    focus = np.abs(frames["field_focus_v_m"])
    early = focus[t < summary["laser_time_s"]].max()
    plane = np.abs(frames["ey_xy"])
    empty = plane[t < summary["laser_time_s"]].max()
    f, i, j = np.unravel_index(plane.argmax(), plane.shape)
    return dict(laser_axis=summary["laser_axis"], seed_electrons=summary["seed_expected_electrons"],
                final_electrons=float(frames["electrons"][-1]), peak_line_ne_m3=float(peak.max()),
                t_line_ne_1e19_ns=None if _first_time(t, peak, 1e19) is None else 1e9 * _first_time(t, peak, 1e19),
                t_line_ne_1e21_ns=None if _first_time(t, peak, 1e21) is None else 1e9 * _first_time(t, peak, 1e21),
                final_focus_field_over_pre_laser=float(focus[-20:].max() / early),
                ionized_length_mm_at_1e19=float((line[-1] >= 1e19).sum()),
                max_field_focus_plane_kv_cm=float(plane.max() / 1e5),
                field_enhancement_over_empty=float(plane.max() / empty),
                max_field_time_ns=float(1e9 * t[f]),
                max_field_location_mm=dict(x=float(i), y_edge_centre=float(j + 0.5)),
                rf_net_work_j=summary["rf_net_work_j"], gas_heat_j=summary["gas_heat_j"],
                rf_energy_residual_rel=summary["rf_energy_residual_rel"],
                plasma_energy_residual_rel=summary["plasma_energy_residual_rel"])


def build(root: Path, stride=2):
    out = root / "orientation"
    cases, table = [], {}
    for field in (50, 200):
        runs = {g: _load(out / f"E{field}_{g}") for g, _ in GEOMETRIES}
        scale = max(np.abs(f["ey_xy"][f["time_s"] < s["laser_time_s"]]).max() for f, s in runs.values())
        ne_hi = max(np.log10(max(f["ne_xy"].max(), f["ne_line"].max(), 1e9)) for f, _ in runs.values())
        entry = dict(field=field, scale_v_m=float(scale), ne_log=[NE_FLOOR, float(ne_hi)], geoms={})
        for g, label in GEOMETRIES:
            f, s = runs[g]
            sel = np.arange(0, len(f["time_s"]), stride)
            table[f"E{field}_{g}"] = metrics(f, s)
            entry["geoms"][g] = dict(
                label=label, summary=dict((k, s[k]) for k in (
                    "laser_axis", "laser_time_s", "laser_waist_m", "laser_energy_j", "seed_expected_electrons",
                    "rf_net_work_j", "gas_heat_j", "rf_energy_residual_rel", "plasma_energy_residual_rel")),
                t_ns=np.round(1e9 * f["time_s"][sel], 4).tolist(),
                electrons=[float(f"{v:.4g}") for v in f["electrons"][sel]],
                focus_kv_cm=np.round(f["field_focus_v_m"][sel] / 1e5, 3).tolist(),
                plane_max_kv_cm=np.round(np.abs(f["ey_xy"][sel]).max(axis=(1, 2)) / 1e5, 1).tolist(),
                laser_w=[float(f"{v:.4g}") for v in f["laser_power_w"][sel]],
                lat_ex=_b64(_q_signed(f["ex_lattice"][sel], scale)),
                lat_ey=_b64(_q_signed(f["ey_lattice"][sel], scale)),
                lat_ez=_b64(_q_signed(f["ez_lattice"][sel], scale)),
                ey_xy=_b64(_q_signed(f["ey_xy"][sel], scale)),
                ne_xy=_b64(_q_log(f["ne_xy"][sel], ne_hi)),
                ne_line=_b64(_q_log(f["ne_line"][sel], ne_hi)),
                ey_line=_b64(_q_signed(f["ey_line"][sel], scale)),
                seed_line=_b64(_q_log(f["seed_density_line_m3"], ne_hi)))
        cases.append(entry)
    f0 = runs["perp"][0]
    geom = dict(node_shape=f0["node_shape"].tolist(), focus=f0["focus_index"].tolist(),
                gas_k=f0["gas_k_range"].tolist(), dx_mm=1.0,
                lat_x=np.round(1e3 * f0["lattice_x_m"], 3).tolist(),
                lat_y=np.round(1e3 * f0["lattice_y_m"], 3).tolist(),
                lat_z=np.round(1e3 * f0["lattice_z_m"], 3).tolist(),
                active_xy=_b64(f0["active_xy"].astype(np.uint8)),
                ey_xy_shape=list(f0["ey_xy"].shape[1:]), ne_xy_shape=list(f0["ne_xy"].shape[1:]))
    data = dict(geometry=geom, cases=cases, metrics=table)
    serial = json.dumps(data, separators=(",", ":"), allow_nan=False)
    payload = base64.b64encode(gzip.compress(serial.encode())).decode()
    (root / "orientation3d.html").write_text(PAGE.replace("__DATA__", payload), encoding="utf-8")
    (root / "orientation_metrics.json").write_text(json.dumps(dict(
        scope="same 3D Maxwell/CR cavity model, seed kernel, mesh and coupling; only the laser axis differs",
        limitations=["plasma can develop only on laser-seeded nodes (no transport); rod growth along E "
                     "away from the seed line is not modeled in 3D",
                     "1 mm mesh: the 20 um seed is node-averaged; field enhancement at sub-mm tips is unresolved, "
                     "so the reported enhancement is a mesh-limited lower bound",
                     "seed kernel is empirical, not a measured multiphoton yield"],
        cases=table), indent=2, allow_nan=False), encoding="utf-8")
    figure(root, out)
    return table


def figure(root, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    colors = dict(perp="C0", par="C3")
    for row, field in enumerate((50, 200)):
        for g, label in GEOMETRIES:
            f, s = _load(out / f"E{field}_{g}")
            t = 1e9 * f["time_s"]
            ax[row, 0].semilogy(t, np.maximum(f["electrons"], 1e-6), color=colors[g], label=label.split(":")[0])
            ax[row, 1].plot(t, f["field_focus_v_m"] / 1e5, color=colors[g], lw=0.6, label=label.split(":")[0])
            pos = np.arange(f["ne_line"].shape[1]) - (f["focus_index"][0] if g == "perp" else f["focus_index"][1] + 0.5)
            ax[row, 2].semilogy(pos, np.maximum(f["ne_line"][-1], 1e8), color=colors[g],
                                label=f"{label.split(':')[0]} (along {'x' if g == 'perp' else 'y'})")
        ax[row, 0].set(title=f"{field} kV/cm: total electrons", xlabel="RF window time (ns)", ylabel="electrons")
        ax[row, 1].set(title=f"{field} kV/cm: Ey at the focus", xlabel="RF window time (ns)", ylabel="Ey (kV/cm)")
        ax[row, 2].set(title=f"{field} kV/cm: ne along the laser line, final frame",
                       xlabel="distance from focus along the beam (mm)", ylabel="ne (m$^{-3}$)")
        for a in ax[row]:
            a.grid(alpha=0.25)
            a.legend(fontsize=8)
    fig.suptitle("Same 3D cavity model; only the laser axis differs. He 200 psig, WR-284, 1 mm mesh")
    fig.tight_layout()
    fig.savefig(root / "orientation_comparison.png", dpi=140)
    plt.close(fig)


PAGE = r"""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Laser orientation vs RF field in 3D: WR-284 cavity</title>
<style>
body{margin:0;background:#0f1522;color:#e5edf7;font:15px/1.45 system-ui}main{max-width:1500px;margin:auto;padding:18px}
a{color:#7dcfff}h1{font-size:24px;margin:6px 0}.notice{background:#263042;padding:10px 12px;border-left:4px solid #e7b650;margin:8px 0}
.controls{display:flex;gap:12px;align-items:center;flex-wrap:wrap;padding:10px 0}input[type=range]{width:360px}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:12px}.panel{background:#1b2333;border-radius:6px;padding:10px}
.panel h2{font-size:16px;margin:0 0 6px}canvas{width:100%;background:#0b1020;border-radius:4px;touch-action:none}
.small{font-size:12.5px;color:#b9c8dc}button,select{background:#39445b;color:#fff;padding:6px 9px;border:1px solid #62728b;border-radius:4px}
table{border-collapse:collapse;font-size:13px;width:100%}td,th{border-bottom:1px solid #34405a;padding:4px 6px;text-align:right}td:first-child,th:first-child{text-align:left}
.key span{display:inline-block;margin-right:14px}.sw{display:inline-block;width:12px;height:12px;border-radius:2px;vertical-align:-1px;margin-right:4px}
@media(max-width:1050px){.grid{grid-template-columns:1fr}}
</style><main>
<p><a href="../">Project viewer</a> · <a href="dynamics.html">Full RF/optics/plasma/PIC dynamics</a> · <a href="orientation_comparison.png">Static comparison figure</a> · <a href="orientation_metrics.json">Metrics</a></p>
<h1>Laser perpendicular vs parallel to the RF field, in 3D</h1>
<div class="notice"><b>Geometry.</b> WR-284: broad dimension <b>a = 72 mm along x</b>, narrow gap <b>b = 34 mm along y</b>, cavity axis along z.
In the TE<sub>10</sub> mode the RF electric field points <b>only along y</b>. It runs across the short 34 mm gap, from one broad wall to the other, and is strongest at the centre: |E| ∝ sin(πx/a).
It forms a standing wave along z. Red arrows point to +y and blue to −y; arrow length is the local |E|.
<b>Perpendicular (experiment, Magid et al. IPSTA slides 9 and 13):</b> the laser enters through the narrow side walls and crosses the guide along x at mid-height. The camera looks down through the broad wall.
<b>Parallel:</b> the laser enters through a hole in the broad wall and runs along y, along E, through the field maximum.</div>
<div class="controls"><label>RF field <select id="case"></select></label><button id="play">Play</button>
<input id="time" type="range" min="0" value="0"><b id="clock"></b>
<label><input type="checkbox" id="lat" checked> 3D field arrows</label><label><input type="checkbox" id="walls" checked> broad walls</label></div>
<div class="key small"><span><i class="sw" style="background:#ff5a4e"></i>E<sub>y</sub> &gt; 0</span><span><i class="sw" style="background:#4f8cff"></i>E<sub>y</sub> &lt; 0</span>
<span><i class="sw" style="background:#4be38a"></i>laser beam (drawn wide; real waist 20 µm)</span><span><i class="sw" style="background:linear-gradient(90deg,#5a2d82,#ff8c1a,#fff6a0)"></i>computed n<sub>e</sub> on the laser line (log)</span>
<span>Drag a 3D view to orbit; both views turn together.</span></div>
<div class="grid">
<section class="panel"><h2 id="h0"></h2><canvas id="v0" width="700" height="470"></canvas></section>
<section class="panel"><h2 id="h1"></h2><canvas id="v1" width="700" height="470"></canvas></section>
<section class="panel"><h2>Front view at the focus plane (as on slide 9): E<sub>y</sub> across the gap</h2><canvas id="f0" width="700" height="360"></canvas><p class="small" id="n0"></p></section>
<section class="panel"><h2>Front view at the focus plane: parallel geometry</h2><canvas id="f1" width="700" height="360"></canvas><p class="small" id="n1"></p></section>
<section class="panel"><h2>Electron density along each laser line, now</h2><canvas id="line" width="700" height="320"></canvas><p class="small">Distance from the focus along the beam. Perpendicular: along x, where the RF field falls as sin(πx/a) away from the centre. Parallel: along y, where the RF field is the same along the whole gap. The dashed curves are the injected laser seed.</p></section>
<section class="panel"><h2>Time histories</h2><canvas id="hist" width="700" height="320"></canvas><p class="small">Solid: total electrons (log). Thin: E<sub>y</sub> at the focus (linear, right axis). Laser peak at 2 ns.</p></section>
<section class="panel" style="grid-column:1/-1"><h2>Matched comparison (same field, mesh, seed kernel and coupling; only the laser axis differs)</h2><div id="table"></div></section>
</div>
<div class="notice small"><b>What this does and does not show.</b> The RF field and plasma are computed together by a 3D Maxwell (Yee FDTD) / collisional-radiative fluid model of the WR-284 cavity at 200 psig He. The plasma current feeds back on the field.
The plasma can only grow on mesh nodes that the laser seeded. Diffusion, drift and streamer growth of a rod <i>along E away from the laser line</i> are not modelled in 3D.
The experiment shows this growth (front-view frames, slides 13 and 14), and the axisymmetric channel-growth model covers it separately.
The 1 mm mesh averages the 20 µm seed; field enhancement at sub-mm channel tips is unresolved. The seed kernel is empirical. This is not a calibrated reproduction of the measured images.</div>
</main><script id="payload" type="text/plain">__DATA__</script><script>
let D,C=0,I=0,playing=false,last=0,yaw=-0.75,pitch=0.42,drag=null;
const el=id=>document.getElementById(id),G=['perp','par'];
const b2=(s,T)=>{const b=Uint8Array.from(atob(s),c=>c.charCodeAt(0));return new T(b.buffer)};
function ne_rgb(q){q=Math.max(0,Math.min(1,q));const a=[[90,45,130],[200,60,90],[255,140,26],[255,246,160]],x=q*3,i=Math.min(2,Math.floor(x)),f=x-i;return a[i].map((v,k)=>Math.round(v+(a[i+1][k]-v)*f))}
function prep(){D.cases.forEach(c=>G.forEach(g=>{const r=c.geoms[g];for(const k of ['lat_ex','lat_ey','lat_ez','ey_xy','ey_line'])r[k]=b2(r[k],Int8Array);for(const k of ['ne_xy','ne_line','seed_line'])r[k]=b2(r[k],Uint8Array)}));D.geometry.active_xy=b2(D.geometry.active_xy,Uint8Array)}
function frameCount(){return D.cases[C].geoms.perp.t_ns.length}
function dims(){const g=D.geometry,[nx,ny,nz]=g.node_shape;return {a:nx-1,b:ny-1,z0:g.gas_k[0],z1:g.gas_k[1],fx:g.focus[0],fy:g.focus[1]+0.5,fz:g.focus[2]}}
function P(x,y,z,W,H){const d=dims(),X=x-d.a/2,Y=y-d.b/2,Z=z-d.fz,u=Math.cos(yaw)*X+Math.sin(yaw)*Z,w=-Math.sin(yaw)*X+Math.cos(yaw)*Z,s=Math.min(W,H)/Math.max(d.a*1.15,(d.z1-d.z0)*0.8);
return [W/2+s*u,H/2+s*(-Math.cos(pitch)*Y+Math.sin(pitch)*w)*1.0,Math.cos(pitch)*w+Math.sin(pitch)*Y,s]}
function seg(g,p,q,col,w){g.strokeStyle=col;g.lineWidth=w;g.beginPath();g.moveTo(p[0],p[1]);g.lineTo(q[0],q[1]);g.stroke()}
function view3d(k){const geo=G[k],c=el('v'+k),g=c.getContext('2d'),W=c.width,H=c.height,d=dims(),cs=D.cases[C],r=cs.geoms[geo],L=D.geometry,pr=(x,y,z)=>P(x,y,z,W,H);
g.fillStyle='#0b1020';g.fillRect(0,0,W,H);
const corners=[];for(const x of [0,d.a])for(const y of [0,d.b])for(const z of [d.z0,d.z1])corners.push([x,y,z]);
if(el('walls').checked)for(const y of [0,d.b]){const q=[[0,y,d.z0],[d.a,y,d.z0],[d.a,y,d.z1],[0,y,d.z1]].map(v=>pr(...v));g.fillStyle=y?'rgba(150,170,200,.07)':'rgba(150,170,200,.12)';g.beginPath();q.forEach((p,i)=>i?g.lineTo(p[0],p[1]):g.moveTo(p[0],p[1]));g.closePath();g.fill()}
g.lineWidth=1;g.strokeStyle='#7e93ad';for(let i=0;i<8;i++)for(let j=i+1;j<8;j++){const a=corners[i],b=corners[j],n=(a[0]!==b[0])+(a[1]!==b[1])+(a[2]!==b[2]);if(n===1)seg(g,pr(...a),pr(...b),'#6f86a3',1)}
const items=[],nx=L.lat_x.length,ny=L.lat_y.length,nz=L.lat_z.length,off=I*nx*ny*nz,len=3.6;
if(el('lat').checked)for(let i=0;i<nx;i++)for(let j=0;j<ny;j++)for(let m=0;m<nz;m++){const o=off+(i*ny+j)*nz+m,ex=r.lat_ex[o]/127,ey=r.lat_ey[o]/127,ez=r.lat_ez[o]/127,mag=Math.hypot(ex,ey,ez);if(mag<0.04)continue;
const x=L.lat_x[i],y=L.lat_y[j],z=L.lat_z[m],p=pr(x-ex*len/2,y-ey*len/2,z-ez*len/2),q=pr(x+ex*len/2,y+ey*len/2,z+ez*len/2);items.push({d:(p[2]+q[2])/2,f:()=>{const col=ey>=0?`rgba(255,90,78,${0.25+0.75*Math.min(1,mag)})`:`rgba(79,140,255,${0.25+0.75*Math.min(1,mag)})`;seg(g,p,q,col,1.4);
const ang=Math.atan2(q[1]-p[1],q[0]-p[0]),h=Math.min(6,3+4*mag);g.fillStyle=col;g.beginPath();g.moveTo(q[0],q[1]);g.lineTo(q[0]-h*Math.cos(ang-0.45),q[1]-h*Math.sin(ang-0.45));g.lineTo(q[0]-h*Math.cos(ang+0.45),q[1]-h*Math.sin(ang+0.45));g.closePath();g.fill()}})}
const beam=geo==='perp'?[[-14,d.fy,d.fz],[d.a+14,d.fy,d.fz]]:[[d.fx,d.b+16,d.fz],[d.fx,0,d.fz]],lp=r.laser_w[I]/Math.max(...r.laser_w),bp=pr(...beam[0]),bq=pr(...beam[1]);
items.push({d:(bp[2]+bq[2])/2-1e-3,f:()=>{seg(g,bp,bq,`rgba(75,227,138,${0.22+0.6*lp})`,7);g.fillStyle='#8ff5b8';g.font='bold 13px system-ui';g.fillText(geo==='perp'?'laser → (+x), ⟂ E':'laser ↓ (−y), ∥ E',bp[0]+4,bp[1]-8)}});
const nl=r.ne_line.length/frameCount(),lo=I*nl;for(let s=0;s<nl;s++){const q8=r.ne_line[lo+s];if(q8<3)continue;const x0=geo==='perp'?s:d.fx,y0=geo==='perp'?d.fy:s,p=pr(geo==='perp'?x0-0.5:x0,geo==='perp'?y0:y0-0.5,d.fz),q=pr(geo==='perp'?x0+0.5:x0,geo==='perp'?y0:y0+0.5,d.fz),rgb=ne_rgb(q8/255);
items.push({d:(p[2]+q[2])/2-2e-3,f:()=>seg(g,p,q,`rgb(${rgb})`,2+7*q8/255)})}
items.sort((u,v)=>v.d-u.d).forEach(it=>it.f());
const o=pr(0,0,d.z0),ax=[[[14,0,0],'x (a = 72 mm)'],[[0,14,0],'y (b = 34 mm), E'],[[0,0,14],'z (cavity axis)']];ax.forEach(([v,t])=>{const q=pr(v[0],v[1],d.z0+v[2]);seg(g,o,q,'#e8eef8',1.6);g.fillStyle='#e8eef8';g.font='12px system-ui';g.fillText(t,q[0]+3,q[1]-3)});
const f=cs.geoms[geo];g.fillStyle='#cfdaea';g.font='12px system-ui';g.fillText(`t = ${f.t_ns[I].toFixed(2)} ns   N_e = ${f.electrons[I].toExponential(2)}   Ey(focus) = ${f.focus_kv_cm[I].toFixed(1)} kV/cm`,10,H-10)}
function front(k){const geo=G[k],c=el('f'+k),g=c.getContext('2d'),W=c.width,H=c.height,d=dims(),r=D.cases[C].geoms[geo],[ex,ey]=D.geometry.ey_xy_shape,[nx,ny]=D.geometry.ne_xy_shape;
g.fillStyle='#0b1020';g.fillRect(0,0,W,H);const s=Math.min((W-80)/d.a,(H-90)/d.b),X0=50,Y0=25,px=x=>X0+s*x,py=y=>Y0+s*(d.b-y);
const img=g.createImageData(ex,ey),oe=I*ex*ey;for(let i=0;i<ex;i++)for(let j=0;j<ey;j++){const v=r.ey_xy[oe+i*ey+j]/127,k4=4*((ey-1-j)*ex+i);img.data[k4]=v>0?255*v:0;img.data[k4+1]=30*Math.abs(v);img.data[k4+2]=v<0?255*-v:0;img.data[k4+3]=255}
const tc=document.createElement('canvas');tc.width=ex;tc.height=ey;tc.getContext('2d').putImageData(img,0,0);g.imageSmoothingEnabled=true;g.drawImage(tc,px(0),py(d.b),s*d.a,s*d.b);
const on=I*nx*ny;for(let i=0;i<nx;i++)for(let j=0;j<ny;j++){const q8=r.ne_xy[on+i*ny+j];if(q8<3)continue;const rgb=ne_rgb(q8/255);g.fillStyle=`rgba(${rgb},${0.35+0.65*q8/255})`;g.fillRect(px(i-0.5),py(j+0.5),s,s)}
for(let i=4;i<ex-1;i+=6)for(let j=3;j<ey;j+=7){const v=r.ey_xy[oe+i*ey+j]/127;if(Math.abs(v)<0.03)continue;const x=px(i),y0=py(j+0.5),h=0.5*6*s*v;seg(g,[x,y0+h],[x,y0-h],'#ffffff',1.2);g.fillStyle='#fff';g.beginPath();g.moveTo(x,y0-h);g.lineTo(x-3,y0-h+(v>0?5:-5));g.lineTo(x+3,y0-h+(v>0?5:-5));g.fill()}
g.strokeStyle='#ffe27a';g.lineWidth=2.2;g.beginPath();for(let i=0;i<ex;i++){const v=Math.abs(r.ey_xy[oe+i*ey+Math.floor(ey/2)])/127,x=px(i),y=py(d.b*0.08)-v*s*d.b*0.8;i?g.lineTo(x,y):g.moveTo(x,y)}g.stroke();
g.strokeStyle='#cfd9e8';g.lineWidth=2;g.strokeRect(px(0),py(d.b),s*d.a,s*d.b);
const lp=r.laser_w[I]/Math.max(...r.laser_w);g.strokeStyle=`rgba(75,227,138,${0.25+0.6*lp})`;g.lineWidth=6;g.beginPath();if(geo==='perp'){g.moveTo(px(-8),py(d.fy));g.lineTo(px(d.a+8),py(d.fy))}else{g.moveTo(px(d.fx),py(d.b+5));g.lineTo(px(d.fx),py(0))}g.stroke();
g.fillStyle='#dfe8f5';g.font='12px system-ui';g.fillText('x: broad dimension a = 72 mm →',px(0),H-12);g.save();g.translate(18,py(0));g.rotate(-Math.PI/2);g.fillText('y: gap b = 34 mm (E) →',0,0);g.restore();
g.fillText('yellow curve: |Ey| along x at mid-gap (computed)',px(d.a)-290,py(d.b)-6);
el('n'+k).textContent=`${D.cases[C].geoms[geo].label}. Colour: signed Ey (red +y, blue −y), scaled to the empty-cavity peak ±${(D.cases[C].scale_v_m/1e5).toFixed(0)} kV/cm and saturating above it. Max |Ey| in this plane now: ${r.plane_max_kv_cm[I].toFixed(0)} kV/cm. White arrows: Ey. Orange: computed ne > 10^8.3 m⁻³ in this plane.`}
function axes(g,W,H,xl,yl){g.strokeStyle='#5f6f88';g.lineWidth=1;g.beginPath();g.moveTo(55,10);g.lineTo(55,H-38);g.lineTo(W-50,H-38);g.stroke();g.fillStyle='#dbe5f3';g.font='12px system-ui';g.fillText(xl,60,H-10);g.save();g.translate(14,H-45);g.rotate(-Math.PI/2);g.fillText(yl,0,0);g.restore()}
function lineplot(){const c=el('line'),g=c.getContext('2d'),W=c.width,H=c.height,d=dims(),cs=D.cases[C];g.fillStyle='#0b1020';g.fillRect(0,0,W,H);axes(g,W,H,'distance from focus along the beam (mm)','log10 ne (m⁻³)');
const span=40,X=v=>55+(W-110)*(v+span)/(2*span),lg=cs.ne_log,Y=q=>H-38-(H-50)*q/255;const cols={perp:'#7ec8ff',par:'#ff7a6b'};
G.forEach(geo=>{const r=cs.geoms[geo],nl=r.seed_line.length,f0=geo==='perp'?d.fx:d.fy;[[r.ne_line.subarray(I*nl,(I+1)*nl),[]],[r.seed_line,[5,4]]].forEach(([a,dash])=>{g.setLineDash(dash);g.strokeStyle=cols[geo];g.lineWidth=dash.length?1:2;g.beginPath();let st=0;for(let s=0;s<nl;s++){const x=s-f0;if(Math.abs(x)>span)continue;const p=[X(x),Y(a[s])];st++?g.lineTo(...p):g.moveTo(...p)}g.stroke()});g.setLineDash([])});
g.fillStyle='#dbe5f3';for(let q=0;q<=4;q++){const v=lg[0]+(lg[1]-lg[0])*q/4;g.fillText(v.toFixed(0),24,Y(255*q/4)+4)}for(let v=-span;v<=span;v+=10)g.fillText(v,X(v)-6,H-24);
g.fillStyle=cols.perp;g.fillText('perpendicular (along x)',W-200,22);g.fillStyle=cols.par;g.fillText('parallel (along y; walls at ±17 mm)',W-200,38)}
function hist(){const c=el('hist'),g=c.getContext('2d'),W=c.width,H=c.height,cs=D.cases[C];g.fillStyle='#0b1020';g.fillRect(0,0,W,H);axes(g,W,H,'RF window time (ns)','log10 electrons');
const t=cs.geoms.perp.t_ns,X=v=>55+(W-110)*v/t[t.length-1];let lo=1e9,hi=-1e9,em=1;G.forEach(geo=>{cs.geoms[geo].electrons.forEach(v=>{const l=Math.log10(Math.max(v,1e-6));lo=Math.min(lo,l);hi=Math.max(hi,l)});em=Math.max(em,...cs.geoms[geo].focus_kv_cm.map(Math.abs))});
const Y=l=>H-38-(H-50)*(l-lo)/(hi-lo||1),YE=e=>H-38-(H-50)*(e/em+1)/2,cols={perp:'#7ec8ff',par:'#ff7a6b'};
G.forEach(geo=>{const r=cs.geoms[geo];g.strokeStyle=cols[geo]+'66';g.lineWidth=0.7;g.beginPath();r.focus_kv_cm.forEach((e,i)=>i?g.lineTo(X(r.t_ns[i]),YE(e)):g.moveTo(X(r.t_ns[i]),YE(e)));g.stroke();
g.strokeStyle=cols[geo];g.lineWidth=2.2;g.beginPath();r.electrons.forEach((v,i)=>{const p=[X(r.t_ns[i]),Y(Math.log10(Math.max(v,1e-6)))];i?g.lineTo(...p):g.moveTo(...p)});g.stroke()});
g.strokeStyle='#fff';g.setLineDash([3,3]);g.beginPath();g.moveTo(X(t[I]),10);g.lineTo(X(t[I]),H-38);g.stroke();g.setLineDash([]);
g.fillStyle='#dbe5f3';for(let q=0;q<=4;q++){const l=lo+(hi-lo)*q/4;g.fillText(l.toFixed(1),22,Y(l)+4)}g.fillText(`±${em.toFixed(0)} kV/cm`,W-48,18);for(let v=0;v<=t[t.length-1];v+=3)g.fillText(v,X(v)-4,H-24);
g.fillStyle=cols.perp;g.fillText('perpendicular',W-200,22);g.fillStyle=cols.par;g.fillText('parallel',W-200,38)}
function table(){const m=D.metrics,f=D.cases[C].field,rows=[['laser axis','laser_axis',null],['injected seed electrons','seed_electrons',3],['total electrons, final','final_electrons',3],['peak ne on the laser line (m⁻³)','peak_line_ne_m3',3],['time line ne ≥ 1e19 m⁻³ (ns)','t_line_ne_1e19_ns',2],['time line ne ≥ 1e21 m⁻³ (ns)','t_line_ne_1e21_ns',2],['length with ne ≥ 1e19 m⁻³ at end (mm)','ionized_length_mm_at_1e19',0],['max |Ey| in the focus plane (kV/cm)','max_field_focus_plane_kv_cm',0],['field enhancement over empty cavity','field_enhancement_over_empty',2],['time of max field (ns)','max_field_time_ns',2],['late focus |Ey| / pre-laser','final_focus_field_over_pre_laser',3],['net RF work on plasma (J)','rf_net_work_j',3],['gas heating (J)','gas_heat_j',3],['plasma energy-ledger residual','plasma_energy_residual_rel',2]];
const A=m[`E${f}_perp`],B=m[`E${f}_par`],fmt=(v,p)=>v===null||v===undefined?'not reached':p===null?v:(Math.abs(v)>=1e4||Math.abs(v)<1e-2&&v!==0)?v.toExponential(p):(+v).toFixed(p);
el('table').innerHTML=`<table><tr><th>${f} kV/cm (empty-cavity focus peak), He 200 psig</th><th>perpendicular (experiment)</th><th>parallel</th></tr>`+rows.map(([n,k,p])=>`<tr><td>${n}</td><td>${fmt(A[k],p)}</td><td>${fmt(B[k],p)}</td></tr>`).join('')+'</table>'}
function draw(){el('time').max=frameCount()-1;I=Math.min(I,frameCount()-1);el('time').value=I;const r=D.cases[C].geoms.perp;el('clock').textContent=`t = ${r.t_ns[I].toFixed(2)} ns (laser peak at 2 ns)`;
G.forEach((geo,k)=>{el('h'+k).textContent='3D: '+D.cases[C].geoms[geo].label;view3d(k);front(k)});lineplot();hist()}
async function init(){const bytes=Uint8Array.from(atob(el('payload').textContent.trim()),x=>x.charCodeAt(0));D=await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).json();prep();
D.cases.forEach((c,i)=>el('case').add(new Option(`${c.field} kV/cm`,i)));C=D.cases.length-1;el('case').value=C;I=Math.floor(frameCount()*0.3);
el('case').onchange=e=>{C=+e.target.value;table();draw()};el('time').oninput=e=>{I=+e.target.value;draw()};el('lat').onchange=draw;el('walls').onchange=draw;
el('play').onclick=()=>{playing=!playing;el('play').textContent=playing?'Pause':'Play'};
['v0','v1'].forEach(id=>{const c=el(id);c.onpointerdown=e=>{drag=[e.clientX,e.clientY];c.setPointerCapture(e.pointerId)};c.onpointerup=()=>drag=null;c.onpointermove=e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*0.008;pitch=Math.max(-1.45,Math.min(1.45,pitch+(e.clientY-drag[1])*0.008));drag=[e.clientX,e.clientY];view3d(0);view3d(1)}});
table();draw();requestAnimationFrame(tick)}
function tick(t){if(playing&&t-last>110){I=(I+1)%frameCount();draw();last=t}requestAnimationFrame(tick)}
init().catch(e=>{el('clock').textContent='Viewer error: '+e.message;console.error(e)});
</script></html>"""


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=ROOT / "results" / "dynamics")
    a = p.parse_args()
    print(json.dumps(build(a.root), indent=2))
