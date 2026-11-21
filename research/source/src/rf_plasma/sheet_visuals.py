"""Self-contained interactive 2D/3D viewer and GIF for the prescribed gas sheet."""

import argparse
import base64
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from .drude import te10_slab_scattering
from .sheet_cavity import C0, laser_front_fraction, laser_intensity


VIEWER = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Helium gas sheet · 3D standing-wave cavity</title>
<style>
:root{color-scheme:dark;font-family:Segoe UI,Arial,sans-serif;background:#0a1422;color:#deebf6}
body{max-width:1480px;margin:0 auto;padding:22px}h1{font-size:26px;margin:5px 0 9px}
p{line-height:1.55}.sub{color:#aac4d9;margin-top:0}.tag{color:#69e3ea}
.controls{display:flex;gap:14px;align-items:center;flex-wrap:wrap;background:#15263b;padding:12px;border-radius:10px;position:sticky;top:0;z-index:2}
button,select{background:#1b4964;color:white;border:1px solid #477890;border-radius:5px;padding:8px 14px;cursor:pointer}
input[type=range]{flex:1;min-width:220px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-top:14px}
.card{background:#101f31;border:1px solid #294158;border-radius:12px;padding:12px;min-width:0}
h2{font-size:17px;margin:2px 0 8px}canvas{width:100%;height:auto;display:block;border-radius:5px}
.metrics{display:grid;grid-template-columns:1fr 1fr;gap:10px;margin:10px 0}.metric{background:#172c42;padding:10px;border-radius:7px}.value{display:block;color:#86e8ea;font-size:20px;margin-top:4px}
.small{font-size:13px;color:#b0c9db}.notice{border-left:3px solid #ffc96c;padding:8px 14px;background:#282719}.links a{color:#8ed8ff}
@media(max-width:850px){.grid{grid-template-columns:1fr}body{padding:10px}}
</style></head><body>
<h1>Standing wave → ionized-helium reflector</h1>
<p class="sub">Full 3D Yee FDTD · WR-284 iris cavity · 1 MPa absolute · weak RF probe · <span class="tag">no metal inserted at the sheet</span></p>
<div class="notice"><b>Laser-induced trigger.</b> A short laser pulse (red) enters through an optical port in the narrow wall at x = 0 and crosses the guide at the speed of light. Gas behind the pulse starts ionizing, building a 6 mm ionized-helium sheet in the pulse's path. The <i>timing and path</i> follow the laser. The <i>growth</i> to 1×10²² m⁻³ (0.6 ns rise) is prescribed, not calculated breakdown physics. The 3D RF field responds self-consistently. Playback slows ×6 while the pulse crosses the guide.</div>
<div class="controls"><button id="play" aria-label="Play or pause animation">Pause</button>
<input id="frame" aria-label="Simulation time" type="range" min="0" max="1" value="0">
<strong id="clock"></strong><label>Playback <select id="speed"><option value=".5">0.5×</option><option selected value="1">1×</option><option value="2">2×</option></select></label>
<button id="reset">Reset</button></div>
<div class="grid">
<section class="card"><h2>No sheet: prefilled standing wave</h2><canvas id="control" width="660" height="510"></canvas></section>
<section class="card"><h2>Laser fired: plasma sheet forms in the same standing wave</h2><canvas id="sheet" width="660" height="510"></canvas></section>
<section class="card"><h2>3D cavity, laser pulse and field slices · drag to orbit</h2><canvas id="three" width="660" height="470" style="touch-action:none"></canvas><div class="small">Colors are the simulated E<sub>y</sub>, not an invented standing-wave drawing. Red is the laser pulse, cyan is the ionized gas it leaves behind, and the wireframe is the cavity wall.</div></section>
<section class="card"><h2>What makes the sheet an effective reflector?</h2>
<div class="metrics">
<div class="metric">Laser pulse<span class="value" id="laser" style="color:#ff7a7a"></span></div>
<div class="metric">Mean sheet electron density<span class="value" id="density"></span></div>
<div class="metric">Isolated-slab reflected power R<span class="value" id="reflection"></span></div>
<div class="metric">Isolated-slab transmitted power T<span class="value" id="transmission"></span></div>
<div class="metric">Isolated-slab absorbed power A<span class="value" id="absorption"></span></div>
</div><canvas id="history" width="660" height="225"></canvas>
<p class="small">R, T and A above are an independent analytic <b>matched-guide slab reference</b> evaluated at the x-averaged density, not the cavity's port S-parameters. The cavity already has a PEC short. Downstream stored energy is the direct shielding diagnostic here.</p>
</section></div>
<div class="card" style="margin-top:14px"><h2>Read this animation correctly</h2>
<ul><li>RF phase is slowed for viewing. Time labels are physical nanoseconds. Both 2D panels use the same fixed color scale.</li>
<li>The laser pulse (illustrative 30 ps FWHM) crosses the 72 mm guide in 0.24 ns, so ionization starts at x = 0 first. The optical port is drawn but not meshed, because the RF wall stays PEC. The optical field is not part of the RF solve; it only sets where and when ionization begins.</li>
<li>The cyan band is the 6 mm ionized-helium sheet at the cavity center. At its final density the collision frequency is much greater than the RF frequency: conductivity and impedance mismatch cause strong reflection even though Re(ε) remains positive.</li>
<li>The 1 mm mesh only marginally resolves the approximately 2 mm field-attenuation length. This is a qualitative 3D cavity demonstration; use the independently checked slab calculation for a reflection reference.</li>
<li>Creating this full-area sheet requires at least <span id="ionenergy"></span> J of ionization energy, before heating and losses. It is not claimed to be achievable with the earlier 30 mJ focused pulse; a sheet-forming (cylindrically focused) beam of that energy is assumed.</li>
<li>Multiple ionization, superelastic exchange and nonlocal reduced radiation transport are tested separately at low density. Missing dense-helium molecular and pressure-broadening data prevent a defensible 1 MPa chemical-lifetime prediction.</li></ul>
<p class="links"><a href="standing_wave.gif">2D animation GIF</a> · <a href="overview.png">Overview figure</a> · <a href="atomic_physics.png">Restricted atomic-physics example</a> · <a href="summary.json">Simulation inputs and metrics</a></p></div>
<script>
const D=__DATA__;
function decode(s){return Int8Array.from(atob(s),x=>x.charCodeAt(0))}
const C=decode(D.control),S=decode(D.sheet),Y=decode(D.yz),nx=D.nx,nz=D.nz,ny=D.ny;
const slider=document.getElementById('frame');slider.max=D.t.length-1;
const lut=D.lut;let index=0,playing=true,yaw=-.48,pitch=.65,last=0;
const cv=document.createElement('canvas');cv.width=nz-D.kmin;cv.height=nx;const ic=cv.getContext('2d');
function rgba(q){return lut[Math.max(0,Math.min(254,q+127))]}
function color(q,alpha=1){const v=rgba(q);return `rgba(${v[0]},${v[1]},${v[2]},${alpha})`}
function field(canvas,arr,label){
 const ctx=canvas.getContext('2d'),im=ic.createImageData(cv.width,nx),offset=index*nx*nz;
 for(let i=0;i<nx;i++)for(let k=D.kmin;k<nz;k++){const a=rgba(arr[offset+i*nz+k]),j=4*((nx-1-i)*cv.width+k-D.kmin);im.data[j]=a[0];im.data[j+1]=a[1];im.data[j+2]=a[2];im.data[j+3]=255}
 ic.putImageData(im,0,0);ctx.fillStyle='#101f31';ctx.fillRect(0,0,660,510);ctx.imageSmoothingEnabled=false;ctx.drawImage(cv,62,18,560,360);
 const z=x=>62+(x-D.zmin)/(D.zmax-D.zmin)*560,xp=x=>378-x/D.width*360;
 if(canvas.id==='sheet'){const zl=z(D.sheetZ-.5*D.thickness),zw=z(D.sheetZ+.5*D.thickness)-zl;ctx.fillStyle='#00dbed';for(let i=0;i<nx;i++){const f=D.fx[index*nx+i];if(f>0){ctx.globalAlpha=.08+.4*f;const y0=xp(Math.min(D.width,i+.5)),y1=xp(Math.max(0,i-.5));ctx.fillRect(zl,y0,zw,y1-y0)}}ctx.globalAlpha=1;laser2d(ctx,z,xp)}
 else ctx.setLineDash([5,5]);
 ctx.strokeStyle='#23e6ee';ctx.lineWidth=2;ctx.strokeRect(z(D.sheetZ-.5*D.thickness),18,z(D.sheetZ+.5*D.thickness)-z(D.sheetZ-.5*D.thickness),360);
 ctx.setLineDash([]);
 ctx.strokeStyle='#b9ccdb';ctx.strokeRect(62,18,560,360);ctx.lineWidth=5;ctx.beginPath();ctx.moveTo(z(D.iris),18);ctx.lineTo(z(D.iris),xp((D.width+20)/2));ctx.moveTo(z(D.iris),xp((D.width-20)/2));ctx.lineTo(z(D.iris),378);ctx.stroke();ctx.lineWidth=1;
 ctx.fillStyle='#d8e8f4';ctx.font='14px Segoe UI';for(let zz=40;zz<=130;zz+=15)ctx.fillText(String(zz),z(zz)-9,400);
 ctx.fillText('z (mm) → propagation from feed',190,426);ctx.fillText('x',22,185);ctx.fillText('(mm)',9,205);
 ctx.fillText('72',30,29);ctx.fillText('0',38,378);ctx.fillText('iris',z(D.iris)-10,448);ctx.fillText(canvas.id==='sheet'?'gas sheet':'reference plane',z(D.sheetZ)-35,448);ctx.fillText('PEC short',z(D.zmax)-67,448);
 if(canvas.id==='sheet'){ctx.fillStyle='#ff6b6b';ctx.fillText('▲ laser enters through an optical port in the x = 0 wall (red tick)',62,467)}else{ctx.fillStyle='#9fb3c4';ctx.fillText('no laser fired in this reference run',62,467)}
 for(let j=0;j<180;j++){ctx.fillStyle=color(Math.round((j/179*2-1)*127));ctx.fillRect(430+j,477,1,10)}
 ctx.fillStyle='#d8e8f4';ctx.font='12px Segoe UI';ctx.fillText('Ey · fixed scale shared by both panels',62,488);ctx.fillText(`−${D.scale.toFixed(2)}`,420,501);ctx.fillText(`+${D.scale.toFixed(2)} MV/m`,557,501);
}
function laserActive(){const x=D.laserX[index];return x>-2*D.L&&x<D.width+2*D.L}
function laser2d(ctx,z,xp){
 const zl=z(D.sheetZ-.5*D.thickness),zr=z(D.sheetZ+.5*D.thickness),zc=(zl+zr)/2;
 ctx.strokeStyle='#ff3b3b';ctx.lineWidth=4;ctx.beginPath();ctx.moveTo(zl-4,378);ctx.lineTo(zr+4,378);ctx.stroke();ctx.lineWidth=1;
 if(!laserActive())return;const xl=D.laserX[index];
 ctx.save();ctx.strokeStyle='rgba(255,80,80,.55)';ctx.setLineDash([4,4]);ctx.beginPath();ctx.moveTo(zc,396);ctx.lineTo(zc,xp(Math.max(0,Math.min(D.width,xl))));ctx.stroke();ctx.setLineDash([]);
 ctx.shadowColor='#ff2020';ctx.shadowBlur=18;
 for(let s=-1.6*D.L;s<=1.6*D.L;s+=.5){const x=xl+s;if(x<0||x>D.width)continue;ctx.fillStyle=`rgba(255,${60+120*Math.exp(-4*Math.LN2*(s/D.L)**2)|0},60,${.95*Math.exp(-4*Math.LN2*(s/D.L)**2)})`;ctx.fillRect(zl-3,xp(x+.5),zr-zl+6,xp(x)-xp(x+.5)+1)}
 ctx.restore();ctx.fillStyle='#ffb0b0';ctx.font='bold 14px Segoe UI';const ly=xp(Math.max(0,Math.min(D.width,xl)));ctx.fillText('laser pulse ↑',zr+10,Math.max(34,Math.min(370,ly+5)));
}
function project(z,x,y){const a=z-(D.zmin+D.zmax)/2,b=x-D.width/2,h=y-D.height/2;const u=Math.cos(yaw)*a-Math.sin(yaw)*b,v=Math.sin(yaw)*a+Math.cos(yaw)*b;return [330+3.9*u,240+3.9*(Math.sin(pitch)*v-Math.cos(pitch)*h),Math.cos(pitch)*v+Math.sin(pitch)*h]}
function poly(ctx,points,fill,stroke){ctx.beginPath();points.forEach((v,i)=>{const p=project(...v);if(i)ctx.lineTo(p[0],p[1]);else ctx.moveTo(p[0],p[1])});ctx.closePath();if(fill){ctx.fillStyle=fill;ctx.fill()}if(stroke){ctx.strokeStyle=stroke;ctx.stroke()}}
function three(){
 const ctx=document.getElementById('three').getContext('2d');ctx.fillStyle='#101f31';ctx.fillRect(0,0,660,470);const faces=[],step=3;
 for(let i=0;i<nx-step;i+=step)for(let k=D.kmin;k<nz-step;k+=step){const pts=[[k,i,D.height/2],[k+step,i,D.height/2],[k+step,i+step,D.height/2],[k,i+step,D.height/2]];faces.push({pts,col:color(S[index*nx*nz+i*nz+k],.90),depth:project(k+step/2,i+step/2,D.height/2)[2]})}
 for(let j=0;j<ny-step;j+=step)for(let k=D.kmin;k<nz-step;k+=step){const pts=[[k,D.width/2,j+.5],[k+step,D.width/2,j+.5],[k+step,D.width/2,j+step+.5],[k,D.width/2,j+step+.5]];faces.push({pts,col:color(Y[index*ny*nz+j*nz+k],.85),depth:project(k+step/2,D.width/2,j+step/2)[2]})}
 faces.sort((a,b)=>a.depth-b.depth);for(const f of faces)poly(ctx,f.pts,f.col,null);
 const z0=D.zmin,z1=D.zmax,a=D.width,b=D.height;
 ctx.lineWidth=1.5;for(const z of [z0,z1])poly(ctx,[[z,0,0],[z,a,0],[z,a,b],[z,0,b]],null,'#8baabc');
 for(const x of [0,a])for(const y of [0,b]){let p=project(z0,x,y),q=project(z1,x,y);ctx.strokeStyle='#8baabc';ctx.beginPath();ctx.moveTo(p[0],p[1]);ctx.lineTo(q[0],q[1]);ctx.stroke()}
 for(const [x0,x1] of [[0,(a-20)/2],[(a+20)/2,a]])poly(ctx,[[D.iris,x0,0],[D.iris,x1,0],[D.iris,x1,b],[D.iris,x0,b]],'rgba(100,125,145,.65)','#9db5c4');
 const zs=D.sheetZ;for(let i=0;i<nx-1;i+=2){const f=D.fx[index*nx+i];if(f>0)poly(ctx,[[zs,i,0],[zs,Math.min(a,i+2),0],[zs,Math.min(a,i+2),b],[zs,i,b]],`rgba(0,225,230,${.08+.5*f})`,null)}
 for(const zz of [zs-D.thickness/2,zs+D.thickness/2])poly(ctx,[[zz,0,0],[zz,a,0],[zz,a,b],[zz,0,b]],null,'rgba(26,216,226,.7)');
 const port=[[zs-4,0,b/2-4],[zs+4,0,b/2-4],[zs+4,0,b/2+4],[zs-4,0,b/2+4]];poly(ctx,port,'rgba(255,60,60,.8)','#ff8080');
 if(laserActive()){const xl=D.laserX[index],p=project(zs,-35,b/2),q=project(zs,Math.max(-35,Math.min(a,xl)),b/2);ctx.save();ctx.strokeStyle='rgba(255,70,70,.7)';ctx.lineWidth=3;ctx.setLineDash([6,4]);ctx.beginPath();ctx.moveTo(p[0],p[1]);ctx.lineTo(q[0],q[1]);ctx.stroke();ctx.setLineDash([]);
  const x0=Math.max(0,xl-D.L/2),x1=Math.min(a,xl+D.L/2);if(x1>x0){ctx.shadowColor='#ff2020';ctx.shadowBlur=20;poly(ctx,[[zs,x0,0],[zs,x1,0],[zs,x1,b],[zs,x0,b]],'rgba(255,70,50,.85)','#ffd0c0')}
  ctx.restore();const r=project(zs,Math.max(-20,Math.min(a,xl)),b);ctx.fillStyle='#ff9090';ctx.font='bold 14px Segoe UI';ctx.fillText('laser pulse',r[0]+6,r[1]-6)}
 ctx.font='14px Segoe UI';ctx.fillStyle='#d7e5ef';ctx.fillText('3D solution: x–z and y–z cut planes',16,28);ctx.fillText('RF cavity walls (wireframe), not a new metal barrier',16,448);
 for(const [text,zz,xx,yy] of [['feed',z0,0,0],['short',z1,a,b],['ionized gas',D.sheetZ,a,b],['optical port',D.sheetZ,0,0]]){const p=project(zz,xx,yy);ctx.fillStyle=text==='ionized gas'?'#65f4eb':text==='optical port'?'#ff8a8a':'#d7e5ef';ctx.fillText(text,p[0]+5,p[1]+(text==='optical port'?16:-8))}
}
function history(){
 const ctx=document.getElementById('history').getContext('2d');ctx.clearRect(0,0,660,225);ctx.strokeStyle='#36556e';ctx.strokeRect(55,20,580,155);ctx.font='13px Segoe UI';ctx.fillStyle='#c9dce9';ctx.fillText('Downstream stored energy / initial control value',56,14);
 for(const [data,col] of [[D.uc,'#efbc62'],[D.us,'#49e1e9']]){ctx.strokeStyle=col;ctx.lineWidth=2;ctx.beginPath();data.forEach((v,i)=>{let x=55+D.t[i]/D.t.at(-1)*580,y=175-Math.min(1.15,v)/1.15*155;if(i)ctx.lineTo(x,y);else ctx.moveTo(x,y)});ctx.stroke()}
 ctx.strokeStyle='#ff5a5a';ctx.lineWidth=2;ctx.beginPath();D.li.forEach((v,i)=>{let x=55+D.t[i]/D.t.at(-1)*580,y=175-v/1.15*155;if(i)ctx.lineTo(x,y);else ctx.moveTo(x,y)});ctx.stroke();
 const x=55+D.t[index]/D.t.at(-1)*580;ctx.strokeStyle='white';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(x,20);ctx.lineTo(x,175);ctx.stroke();ctx.fillStyle='#efbc62';ctx.fillText('no sheet',60,200);ctx.fillStyle='#49e1e9';ctx.fillText('laser + gas sheet',150,200);ctx.fillStyle='#ff6b6b';ctx.fillText('laser pulse in guide',290,200);ctx.fillStyle='#c9dce9';ctx.fillText('time (ns) →',540,200);
}
function draw(){slider.value=index;document.getElementById('clock').textContent=D.t[index].toFixed(3)+' ns'+(laserActive()?'  · slow motion':'');field(document.getElementById('control'),C);field(document.getElementById('sheet'),S);three();history();document.getElementById('density').textContent=(D.fraction[index]*D.density).toExponential(2)+' m⁻³';
 const xl=D.laserX[index];document.getElementById('laser').textContent=xl<-1.5*D.L?`armed · fires at ${D.trigger.toFixed(2)} ns`:xl<=D.width+1.5*D.L?`crossing guide · x ≈ ${Math.max(0,Math.min(D.width,xl)).toFixed(0)} mm`:(D.fraction[index]<.999?'passed · plasma growing':'passed · sheet fully formed');
 for(const [id,key] of [['reflection','R'],['transmission','T'],['absorption','A']])document.getElementById(id).textContent=(100*D[key][index]).toPrecision(3)+' %'}
document.getElementById('ionenergy').textContent=D.ionenergy.toFixed(3);
function toggle(){playing=!playing;document.getElementById('play').textContent=playing?'Pause':'Play'}
document.getElementById('play').onclick=toggle;slider.oninput=()=>{index=+slider.value;if(playing)toggle();draw()};
document.getElementById('reset').onclick=()=>{index=0;draw()};
let drag=null;const threeCanvas=document.getElementById('three');threeCanvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];threeCanvas.setPointerCapture(e.pointerId)};threeCanvas.onpointerup=()=>drag=null;threeCanvas.onpointermove=e=>{if(drag){yaw+=(e.clientX-drag[0])*.008;pitch=Math.max(.15,Math.min(1.4,pitch+(e.clientY-drag[1])*.008));drag=[e.clientX,e.clientY];three()}};
function tick(now){if(playing&&now-last>(laserActive()?6:1)*1000/(24*+document.getElementById('speed').value)){index=(index+1)%D.t.length;draw();last=now}requestAnimationFrame(tick)}draw();requestAnimationFrame(tick);
window.viewerState=()=>({frame:index,playing,frames:D.t.length,time_ns:D.t[index],nx,nz});
</script></body></html>"""


def _encoded(a, scale):
    return base64.b64encode(np.rint(a / scale * 127).clip(-127, 127).astype(np.int8).tobytes()).decode("ascii")


def _laser_frames(meta, t, nx):
    """Per-frame local ionized fraction along x, mean fraction, pulse position (mm) and envelope."""
    c = meta["case"]
    x = np.arange(nx) * meta["grid"]["dx_m"]
    fx = laser_front_fraction(t[:, None], x[None, :], c["trigger_s"], c["rise_s"])
    xg = np.linspace(0, meta["width_m"], 145)
    li = laser_intensity(t[:, None], xg[None, :], c["trigger_s"], c["laser_fwhm_s"]).max(axis=1)
    return fx, fx.mean(axis=1), C0 * (t - c["trigger_s"]) * 1e3, li


def make(output):
    meta = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    co, sh = (np.load(output / f"{name}_fields.npz") for name in ("control", "sheet"))
    hist = np.load(output / "histories.npz")
    t = sh["t"]
    fx, frac, laser_x, li = _laser_frames(meta, t, sh["ey_xz"].shape[1])
    ex0, ex1, eyz = co["ey_xz"], sh["ey_xz"], sh["ey_yz"]
    scale = float(max(np.abs(ex0).max(), np.abs(ex1).max(), np.abs(eyz).max()))
    lut = np.rint(plt.get_cmap("RdBu_r")(np.linspace(0, 1, 255))[:, :3] * 255).astype(int).tolist()
    rta = []
    for ff in frac:
        r, tr = te10_slab_scattering(meta["frequency_hz"], meta["width_m"],
                                    meta["sheet_thickness_m"], meta["case"]["density_m3"] * ff,
                                    meta["collision_frequency_s"])
        rta.append([abs(r)**2, abs(tr)**2, 1 - abs(r)**2 - abs(tr)**2])
    rta = np.asarray(rta)
    normal = float(hist["control_downstream"][0])
    data = dict(control=_encoded(ex0, scale), sheet=_encoded(ex1, scale), yz=_encoded(eyz, scale),
                nx=ex0.shape[1], nz=ex0.shape[2], ny=eyz.shape[1], kmin=40,
                zmin=40, zmax=meta["cavity_end_m"] * 1e3, iris=meta["cavity_entrance_m"] * 1e3,
                sheetZ=meta["sheet_center_m"] * 1e3, thickness=meta["sheet_thickness_m"] * 1e3,
                width=meta["width_m"] * 1e3, height=meta["height_m"] * 1e3,
                density=meta["case"]["density_m3"], fraction=frac.tolist(),
                fx=np.round(fx, 3).ravel().tolist(), laserX=np.round(laser_x, 2).tolist(),
                li=np.round(li, 3).tolist(), L=C0 * meta["case"]["laser_fwhm_s"] * 1e3,
                trigger=meta["case"]["trigger_s"] * 1e9,
                t=(t * 1e9).tolist(), scale=scale / 1e6, lut=lut,
                R=rta[:, 0].tolist(), T=rta[:, 1].tolist(), A=rta[:, 2].tolist(),
                uc=(hist["control_downstream"] / normal).tolist(),
                us=(hist["sheet_downstream"] / normal).tolist(),
                ionenergy=meta["minimum_sheet_ionization_energy_J"])
    (output / "standing_wave.html").write_text(VIEWER.replace("__DATA__", json.dumps(data, separators=(",", ":"))),
                                              encoding="utf-8")
    overview(output, meta, co, sh, hist)
    gif(output, meta, co, sh, scale)
    atomic_figure(output)
    return dict(frames=len(t), html_bytes=(output / "standing_wave.html").stat().st_size,
                color_scale_V_m=scale)


def overview(output, meta, co, sh, hist):
    fig, axs = plt.subplots(2, 2, figsize=(13, 8), layout="constrained")
    extent = (40, 132, 0, 72)
    scale = max(abs(co["ey_xz"]).max(), abs(sh["ey_xz"]).max()) / 1e6
    for ax, arr, title in zip(axs[0], (co["ey_xz"], sh["ey_xz"]), ("No laser: standing wave", "Laser-induced gas sheet: downstream shielding")):
        im = ax.imshow(arr[-1, :, 40:] / 1e6, origin="lower", aspect="auto", extent=extent,
                       cmap="RdBu_r", vmin=-scale, vmax=scale)
        if title.startswith("Laser"):
            ax.axvspan((meta["sheet_center_m"] - .003) * 1e3,
                       (meta["sheet_center_m"] + .003) * 1e3, color="cyan", alpha=.18)
        else:
            ax.axvline(meta["sheet_center_m"] * 1e3, color="teal", ls=":", alpha=.5)
        ax.set(xlabel="z (mm)", ylabel="x (mm)", title=title)
    fig.colorbar(im, ax=axs[0], label="Instantaneous Ey (MV/m), common scale")
    ts = hist["sheet_t"] * 1e9
    for mode, label in (("control", "No laser"), ("sheet", "Laser-induced gas sheet")):
        axs[1, 0].semilogy(ts, hist[f"{mode}_downstream"] / hist["control_downstream"][0], label=label)
    axs[1, 0].axvspan(meta["case"]["trigger_s"] * 1e9, (meta["case"]["trigger_s"] + meta["laser"]["crossing_time_s"]) * 1e9,
                      color="red", alpha=.35, label="Laser pulse crosses guide")
    axs[1, 0].set(xlabel="Time (ns)", ylabel="Downstream EM energy / initial control")
    axs[1, 0].legend()
    scan = np.load(output / "analytic_density_scan.npy")
    for j, label in enumerate(("Reflected R", "Transmitted T", "Absorbed A"), 1):
        axs[1, 1].semilogx(scan[:, 0], scan[:, j], label=label)
    axs[1, 1].axvline(meta["case"]["density_m3"], color="gray", ls="--")
    axs[1, 1].set(xlabel="Prescribed electron density (m$^{-3}$)", ylabel="Power fraction",
                   title="Independent analytic matched-guide slab (not cavity S)")
    axs[1, 1].legend()
    fig.suptitle("3D WR-284 cavity, 1 MPa helium: laser-induced 6 mm ionized-gas sheet (laser-timed, prescribed growth), not a metal insert")
    fig.savefig(output / "overview.png", dpi=170)
    plt.close(fig)


def gif(output, meta, co, sh, scale):
    font = ImageFont.truetype(r"C:\Windows\Fonts\segoeui.ttf", 20)
    small = ImageFont.truetype(r"C:\Windows\Fonts\segoeui.ttf", 16)
    cmap = plt.get_cmap("RdBu_r")
    images = []
    durations = []
    t = sh["t"]
    nx = sh["ey_xz"].shape[1]
    fx, frac, laser_x, _ = _laser_frames(meta, t, nx)
    pulse_mm = C0 * meta["case"]["laser_fwhm_s"] * 1e3
    width_mm = meta["width_m"] * 1e3
    crossing = (laser_x > -2 * pulse_mm) & (laser_x < width_mm + 2 * pulse_mm)
    # Every third stored frame, plus every frame while the laser crosses (GIF is a preview, not raw data).
    for i in sorted(set(range(0, len(t), 3)) | set(np.flatnonzero(crossing).tolist())):
        im = Image.new("RGB", (1100, 610), "#0b1725")
        d = ImageDraw.Draw(im, "RGBA")
        d.text((25, 16), "3D FDTD: laser-induced ionized-helium sheet in a WR-284 standing wave", fill="white", font=font)
        d.text((25, 48), "1 MPa He | laser sets ionization onset; growth to 1e22 m^-3 prescribed | no metal inserted", fill="#b9d1e4", font=small)
        for col, (a, title) in enumerate(((co["ey_xz"], "No laser (reference)"), (sh["ey_xz"], "Laser fired"))):
            arr = np.uint8(cmap((a[i, :, 40:][::-1] / scale + 1) / 2)[:, :, :3] * 255)
            tile = Image.fromarray(arr).resize((510, 370), resample=Image.Resampling.NEAREST)
            x = 25 + 550 * col
            im.paste(tile, (x, 126))
            d.text((x, 92), title, fill="white", font=font)
            zs = (meta["sheet_center_m"] * 1e3 - 40) / 92 * 510 + x
            half = 3 / 92 * 510
            row = lambda xm: 496 - xm / width_mm * 370
            if col == 1:
                for k in range(nx):
                    if fx[i, k] > 0:
                        d.rectangle((zs - half, row(min(width_mm, k + .5)), zs + half, row(max(0, k - .5))),
                                    fill=(0, 225, 235, int(40 + 110 * fx[i, k])))
                d.rectangle((zs - half - 3, 494, zs + half + 3, 499), fill="#ff3b3b")
                if crossing[i]:
                    for s in np.arange(-1.6 * pulse_mm, 1.6 * pulse_mm, .5):
                        xm = laser_x[i] + s
                        if 0 <= xm <= width_mm:
                            w = math.exp(-4 * math.log(2) * (s / pulse_mm) ** 2)
                            d.rectangle((zs - half - 3, row(xm + .5), zs + half + 3, row(xm)),
                                        fill=(255, int(60 + 120 * w), 60, int(245 * w)))
                    ly = row(min(width_mm, max(0, laser_x[i])))
                    d.text((zs + half + 8, max(130, min(480, ly - 10))), "laser pulse ^", fill="#ffb0b0", font=small)
            d.rectangle((zs - half, 126, zs + half, 496), outline="#00e8ed", width=2)
            d.text((x, 503), "z: 40 -> 132 mm | x: 0 -> 72 mm | Ey common scale", fill="#c5dbe8", font=small)
        state = ("laser crossing the guide" if crossing[i] else
                 "laser armed" if laser_x[i] < 0 else "laser passed")
        d.text((25, 545), f"Time: {t[i]*1e9:.3f} ns  |  {state}  |  mean sheet ne = {frac[i]*meta['case']['density_m3']:.2e} m^-3",
               fill="#75edf0", font=font)
        d.text((25, 577), "Laser enters at x = 0 (red port), ionizes the gas behind it; the plasma then reflects the RF and shields the downstream cavity.", fill="#c5dbe8", font=small)
        images.append(im)
        durations.append(220 if crossing[i] else 70)
    images[0].save(output / "standing_wave.gif", save_all=True, append_images=images[1:],
                   duration=durations, loop=0, optimize=False)


def atomic_figure(output):
    data = np.load(output / "atomic_demonstration" / "restricted_cr.npz")
    h, t = data["states"], data["t"] * 1e6
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for i, label in enumerate(("He ground", "He 2 singlet S", "He 2 singlet P", "He+", "He++ bare ion")):
        axs[0, 0].semilogy(t, np.maximum(h[i].mean(axis=0), 1e-10), label=label)
    axs[0, 0].set(xlabel="Time (microseconds)", ylabel="Nuclei fraction")
    axs[0, 0].legend(fontsize=8)
    te = h[5] / (1.5 * (h[3] + 2 * h[4]))
    axs[0, 1].plot(t, te.mean(axis=0))
    axs[0, 1].set(xlabel="Time (microseconds)", ylabel="Mean electron temperature (eV)")
    axs[1, 0].plot(t, h[6].mean(axis=0))
    axs[1, 0].set(xlabel="Time (microseconds)", ylabel="Escaped radiation (eV per nucleus)")
    axs[1, 1].plot(data["z"] * 1e3, data["escape"], marker=".")
    axs[1, 1].set(xlabel="z (mm)", ylabel="Resonance escape probability per emission")
    fig.suptitle("Separate restricted atomic-He example at 0.414 Pa, NOT the 1 MPa sheet lifetime\n"
                 "Sequential ionization + superelastic exchange + nonlocal resonance redistribution")
    fig.savefig(output / "atomic_physics.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    print(json.dumps(make(args.output.resolve()), indent=2))
