"""Standalone, provenance-labelled RF/optics/fluid/afterglow/PIC dynamics viewer."""
from __future__ import annotations

import base64
import gzip
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def read_arrays(path):
    with np.load(path) as data:
        result = {k: data[k].tolist() for k in data.files}
        if any(not np.isfinite(data[k]).all() for k in data.files):
            raise ValueError(f"non-finite viewer data in {path}")
    return result


def build(root, pic_path):
    cases = []
    for label in ("E50", "E200"):
        cavity = root / "cavity" / ("E200_b1" if label == "E200" else "E50_verified")
        gas = root / "afterglow" / ("E200_b1_n160" if label == "E200" else "E50_n160_verified")
        cases.append(dict(label=label,
                          summary=json.loads((cavity / "summary.json").read_text()),
                          rf=read_arrays(cavity / "frames.npz"),
                          gas=read_arrays(gas / "afterglow.npz"),
                          gas_summary=json.loads((gas / "summary.json").read_text())))
    pic = json.loads(pic_path.read_text())
    if not pic.get("cases"):
        raise ValueError("no executed PIC cases for the viewer")
    for case in pic["cases"]:
        if not case.get("particles"):
            raise ValueError(f"PIC case has no actual particle frames: {case.get('label')}")
        if not case.get("eedf_probability") or not case.get("energy_bins_ev"):
            raise ValueError(f"PIC case has no EEDF data: {case.get('label')}")
    data = dict(cases=cases, pic=pic)
    serial = json.dumps(data, separators=(",", ":"), allow_nan=False)
    (root / "dynamics_data.json").write_text(serial, encoding="utf-8")
    payload = base64.b64encode(gzip.compress(serial.encode())).decode()
    page = PAGE.replace("__DATA__", payload)
    (root / "dynamics.html").write_text(page, encoding="utf-8")
    figures(data, root)
    validation_figures(root, pic)
    return data


def validation_figures(root, pic):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cavity = root / "cavity"
    base = read_arrays(cavity / "E200_verified" / "frames.npz")
    fine = read_arrays(cavity / "E200_b2" / "frames.npz")
    bs = json.loads((cavity / "E200_verified" / "summary.json").read_text())
    fs = json.loads((cavity / "E200_b2" / "summary.json").read_text())
    empty = json.loads((cavity / "E200_empty" / "summary.json").read_text())
    single = read_arrays(cavity / "E200_b1" / "frames.npz")
    single_summary = json.loads((cavity / "E200_b1" / "summary.json").read_text())
    t = np.asarray(base["probe_time_s"])
    common = t <= fine["probe_time_s"][-1]
    difference = np.asarray(base["probe_ey_v_m"])[common] - np.interp(
        t[common], fine["probe_time_s"], fine["probe_ey_v_m"])
    rf_nrmse = float(np.linalg.norm(difference) / np.linalg.norm(np.asarray(base["probe_ey_v_m"])[common]))
    common_end = min(base["time_s"][-1], fine["time_s"][-1])
    n_base = np.interp(common_end, base["time_s"], base["electrons"])
    n_fine = np.interp(common_end, fine["time_s"], fine["electrons"])
    rf_metrics = dict(focus_field_nrmse_b4_vs_b2=rf_nrmse,
                      electron_count_relative_difference=float(abs(n_base-n_fine)/n_fine),
                      net_work_relative_difference=abs(bs["rf_net_work_j"]-fs["rf_net_work_j"])/abs(fs["rf_net_work_j"]),
                      gas_heat_relative_difference=abs(bs["gas_heat_j"]-fs["gas_heat_j"])/abs(fs["gas_heat_j"]),
                      empty_control_residual_j=empty["rf_energy_residual_j"],
                      empty_control_residual_relative_em=empty["rf_energy_residual_rel"],
                      empty_control_seed_count=empty["seed_expected_electrons"],
                      baseline_subtracted_residual_relative_work=abs(
                          bs["rf_energy_residual_j"]-empty["rf_energy_residual_j"])/abs(bs["rf_net_work_j"]),
                      caveat="RF staggered-energy/flux diagnostic is not an exact discrete identity; weak loading is unresolved")
    tf = np.asarray(fine["probe_time_s"])
    common_single = tf <= single["probe_time_s"][-1]
    ey_single = np.interp(tf[common_single], single["probe_time_s"], single["probe_ey_v_m"])
    rf_metrics["focus_field_nrmse_b2_vs_b1"] = float(np.linalg.norm(
        np.asarray(fine["probe_ey_v_m"])[common_single]-ey_single)/np.linalg.norm(ey_single))
    rf_metrics["primary_strong_case"] = "E200_b1 (one FDTD step per kinetic coupling)"
    rf_metrics["primary_summary"] = single_summary
    meshes = []
    for n, name in ((160, "E200_b1_n160"), (320, "E200_b1_n320"), (640, "E200_b1_n640")):
        path = root / "afterglow" / name
        meshes.append((n, read_arrays(path / "afterglow.npz"), json.loads((path / "summary.json").read_text())))
    gas_metrics = dict(meshes=[s for _, _, s in meshes])
    peaks = [s["maximum_velocity_m_s"] for _, _, s in meshes]
    gas_metrics["peak_velocity_relative_difference_320_to_640"] = abs(peaks[1]-peaks[2])/peaks[2]
    gas_metrics["peak_velocity_converged_5_percent"] = gas_metrics["peak_velocity_relative_difference_320_to_640"] < 0.05
    high = meshes[-1][1]
    r = np.asarray(high["r_m"])
    weights = r
    p0 = high["pressure"][0][-1]
    gas_metrics["final_pressure_perturbation_l1_relative_320_to_640"] = float(
        np.sum(np.abs(np.interp(r, meshes[1][1]["r_m"], meshes[1][1]["pressure"][-1])-high["pressure"][-1])*weights)
        / np.sum(np.abs(np.asarray(high["pressure"][-1])-p0)*weights))
    ss = json.loads((root / "afterglow" / "E200_b1_source_refined" / "summary.json").read_text())
    ts = json.loads((root / "afterglow" / "E200_b1_trap10" / "summary.json").read_text())
    gas_metrics["half_source_step_peak_velocity_relative_difference"] = abs(peaks[0]-ss["maximum_velocity_m_s"])/ss["maximum_velocity_m_s"]
    gas_metrics["trapping_radius_5_to_10_mm_radiative_energy_relative_difference"] = abs(
        meshes[0][2]["radiative_energy_j_per_m"]-ts["radiative_energy_j_per_m"])/ts["radiative_energy_j_per_m"]
    gas_metrics["caveat"] = "velocity peak is mesh-dependent; 1-mm projected source is not an optical-filament resolution study"
    pic_metrics = [dict(label=c["label"], mean_energy_nrmse=c["validation"]["fluid_energy_nrmse_last_cycle"],
                        ionization_events=c["validation"]["ionization_events"],
                        warnings=c["validation"]["validation_warnings"]) for c in pic["cases"]]
    metrics = dict(scope="restricted numerical diagnostics, not complete or experimental validation",
                   rf=rf_metrics, gas=gas_metrics, pic=pic_metrics,
                   shock=json.loads((root / "gas_validation.json").read_text()))
    (root / "validation.json").write_text(json.dumps(metrics, indent=2, allow_nan=False), encoding="utf-8")
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    ax[0, 0].plot(t[common]*1e9, difference/1e5)
    ax[0, 0].set(title=f"RF: b4/b2 {100*rf_nrmse:.2f}%, b2/b1 {100*rf_metrics['focus_field_nrmse_b2_vs_b1']:.2f}%",
                 xlabel="RF window time (ns)", ylabel="Ey b4 - b2 (kV/cm)")
    for label, result in (("block 4", base), ("block 2", fine), ("block 1 (primary)", single)):
        ax[0, 1].semilogy(np.asarray(result["time_s"])*1e9, result["electrons"], label=label)
    ax[0, 1].set(title="Perpendicular-source avalanche refinement", xlabel="RF window time (ns)", ylabel="total electrons")
    for n, result, summary in meshes:
        rr = np.asarray(result["r_m"])*1e3
        ax[0, 2].plot(rr, np.asarray(result["pressure"][-1])/p0, label=f"{n} cells")
        ax[1, 0].plot(rr, np.asarray(result["gas_temperature"][-1])-300, label=f"{n} cells")
        ax[1, 1].plot(rr, result["velocity"][-1], label=f"{n} cells")
        ax[1, 2].plot(np.asarray(result["time_s"])*1e6, result["radiative_energy_j_per_m"], label=f"{n} cells")
    ax[0, 2].set(title="Local pressure at 1 us", xlabel="radius (mm)", ylabel="p / ambient p", xlim=(0, 2.5))
    ax[1, 0].set(title="Local gas heating at 1 us", xlabel="radius (mm)", ylabel="Tg - 300 K", xlim=(0, 2.5))
    ax[1, 1].set(title="Velocity peak NOT mesh-converged", xlabel="radius (mm)", ylabel="outward velocity (m/s)", xlim=(0, 2.5))
    ax[1, 2].set(title="Fixed-radius CR radiation; uncertain rates", xlabel="local RF-off time (us)", ylabel="modeled radiative energy (J/m)")
    for a in ax.flat:
        a.grid(alpha=0.25)
        if a.lines and len(a.lines) > 1:
            a.legend(fontsize=8)
    fig.suptitle("Validation distinguishes conservation, convergence and missing physics")
    fig.tight_layout()
    fig.savefig(root / "validation.png", dpi=140)
    plt.close(fig)


def figures(data, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    for case in data["cases"]:
        rf, gas = case["rf"], case["gas"]
        t = np.asarray(rf["time_s"]) * 1e9
        probe_t = np.asarray(rf["probe_time_s"]) * 1e9
        ax[0, 0].plot(probe_t, np.asarray(rf["probe_ey_v_m"]) / 1e5, label=case["label"])
        ax[0, 0].set_xlim(0.5, 3.5)
        ax[0, 1].semilogy(t, np.maximum(rf["electrons"], 1e-30), label=case["label"])
        ax[0, 2].plot(t, rf["gas_heat_j"], label=f"{case['label']} heavy gas")
        ax[0, 2].plot(t, rf["sink_j"], "--", label=f"{case['label']} CR sink")
        r = np.asarray(gas["r_m"]) * 1e3
        p0 = gas["pressure"][0][-1]
        for i in (0, len(gas["time_s"]) // 2, len(gas["time_s"]) - 1):
            ax[1, 0].plot(r, np.asarray(gas["pressure"][i]) / p0,
                          label=f"{case['label']} {gas['time_s'][i]*1e6:.2f} us")
        ax[1, 1].plot(np.asarray(gas["time_s"]) * 1e6,
                      np.max(np.abs(gas["velocity"]), axis=1), label=case["label"])
        ax[1, 2].plot(np.asarray(gas["time_s"]) * 1e6,
                      gas["sink_energy_j_per_m"], label=f"{case['label']} total sink")
        ax[1, 2].plot(np.asarray(gas["time_s"]) * 1e6,
                      gas["radiative_energy_j_per_m"], "--", label=f"{case['label']} modeled radiation")
    labels = [("RF window time (ns)", "Ey at focus (kV/cm)", "Actual transient Yee field"),
              ("RF window time (ns)", "total electrons", "Empirical saturated perpendicular seed"),
              ("RF window time (ns)", "energy (J)", "Gas heat and CR energy sink (not all photons)"),
              ("radius (mm)", "pressure / ambient", "Local RF-off expansion: coarse 3D source"),
              ("local afterglow time (us)", "maximum |velocity| (m/s)", "Computed gas motion"),
              ("local afterglow time (us)", "energy sink (J/m)", "Radiation separated from unresolved disposal")]
    for a, (x, y, title) in zip(ax.flat, labels):
        a.set(xlabel=x, ylabel=y, title=title)
        a.grid(alpha=0.25)
        a.legend(fontsize=7)
    fig.suptitle("Restricted multiphysics study: 1.48 MPa He; not experimental validation")
    fig.tight_layout()
    fig.savefig(out / "dynamics_overview.png", dpi=140)
    plt.close(fig)
    fig, a = plt.subplots(figsize=(12, 4.5))
    a.set(xlim=(0, 1), ylim=(0, 1))
    a.axis("off")
    for x, y, text, color in (
            (0.14, 0.76, "Prescribed Gaussian\nEmpirical seed along x", "#fff2cd"),
            (0.46, 0.76, "3D Maxwell + CR fluid\nRF window: 0-18 ns", "#deecff"),
            (0.80, 0.76, "Conservative projection\nAssumed RF quench", "#deecff"),
            (0.80, 0.25, "Local Euler + CR gas\nAfterglow: 0-1 us", "#deecff"),
            (0.14, 0.25, "Independent real PIC\nImposed 50/100 kV/cm", "#def4e6"),
            (0.46, 0.25, "Weighted EEDF / particles\nPassed AND failed checks", "#def4e6")):
        a.text(x, y, text, ha="center", va="center", fontsize=12,
               bbox=dict(boxstyle="round,pad=0.6", facecolor=color, edgecolor="#465466"))
    for start, end in (((0.26, 0.76), (0.33, 0.76)), ((0.60, 0.76), (0.67, 0.76)),
                       ((0.80, 0.63), (0.80, 0.39)), ((0.27, 0.25), (0.32, 0.25))):
        a.annotate("", end, start, arrowprops=dict(arrowstyle="->", lw=1.8))
    a.text(0.03, 0.48, "PIC is not driven by exported cavity fields; its clock is independent.", fontsize=11)
    fig.suptitle("Schematic: three separate calculations, not full-device PIC or the measured discharge")
    fig.tight_layout()
    fig.savefig(out / "workflow.png", dpi=150)
    plt.close(fig)


PAGE = r"""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>RF, optics, plasma, gas and PIC dynamics</title>
<style>
body{margin:0;background:#101522;color:#e5edf7;font:15px/1.45 system-ui}main{max-width:1450px;margin:auto;padding:20px}
a{color:#7dcfff}h1{font-size:25px}.notice{background:#283040;padding:12px;border-left:4px solid #e7b650}
.controls{display:flex;gap:12px;align-items:center;flex-wrap:wrap;padding:15px 0}input[type=range]{width:330px}
.grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px}.panel{background:#1c2434;border-radius:6px;padding:12px}
.panel h2{font-size:17px;margin:0 0 6px}canvas{width:100%;height:250px;background:#0c1120}.small{font-size:12px;color:#b6c6db}
button,select{background:#39445b;color:white;padding:7px;border:1px solid #62728b;border-radius:4px}
@media(max-width:1000px){.grid{grid-template-columns:repeat(2,minmax(0,1fr))}}@media(max-width:650px){.grid{grid-template-columns:1fr}}
</style><main>
<p><a href="../">Project viewer</a> · <a href="dynamics_overview.png">Overview</a> · <a href="validation.png">Refinement diagnostics</a> · <a href="shock_validation.png">Sod shock benchmark</a> · <a href="pic_diagnostics.png">PIC checks</a> · <a href="dynamics_data.json">Viewer data</a> · <a href="validation.json">Validation metrics</a></p>
<h1>RF → electron avalanche → gas expansion and radiation</h1>
<div class="notice"><b>Three distinct calculations, not one full-device PIC simulation.</b>
The short RF stage is a transient 3D Maxwell/CR-fluid calculation with a number-conserving empirical seed along x, perpendicular to Ey.
Optics is a prescribed vacuum Gaussian; physical narrow-beam attenuation and refraction remain unverified. The subsequent local cylindrical afterglow assumes RF quench and couples gas motion to CR chemistry.
PIC displays real particles from a separate imposed-field, fixed-neutral local calculation. Its clock is independent.</div>
<div class="controls"><select id="case"></select><select id="stage"><option value="rf">3D RF / laser / fluid stage</option><option value="gas">Local RF-off afterglow</option></select>
<button id="play">Play</button><input id="time" type="range" min="0" value="0"><b id="clock"></b></div>
<div class="grid">
<section class="panel"><h2>RF electric field Ey</h2><canvas id="rf" width="500" height="300"></canvas><p class="small" id="rfnote"></p></section>
<section class="panel"><h2>Fluid electron density</h2><canvas id="ne" width="500" height="300"></canvas><p class="small" id="nenote"></p></section>
<section class="panel"><h2>Laser focal optics</h2><canvas id="optics" width="500" height="300"></canvas><p class="small" id="optnote"></p></section>
<section class="panel"><h2>Radiation / optical emission</h2><canvas id="uv" width="500" height="300"></canvas><p class="small" id="uvnote"></p></section>
<section class="panel"><h2>Heavy gas: expansion and compression</h2><canvas id="gas" width="500" height="300"></canvas><p class="small" id="gasnote"></p></section>
<section class="panel"><h2>Real local PIC macroparticles</h2><select id="piccase"></select><button id="picplay">Play PIC</button>
<input id="pictime" type="range" min="0" value="0" style="width:100%">
<canvas id="pic" width="500" height="300" style="height:180px"></canvas>
<canvas id="eedf" width="500" height="180" style="height:150px"></canvas><p class="small" id="picnote"></p><p id="picwarning" style="color:#ffbd80"></p></section>
</div><p id="stats"></p>
<div class="notice">An emission estimate is not a calibrated camera image. UV panels use lumped excited-He/excimer escape rates.
Visible 400–700 nm free-free estimates assume Maxwellian electrons and optical thinness; helium visible lines and pressure quenching are not resolved.
Radiative energy is distinguished from the total CR sink, which includes phenomenological recombination disposal. Escape uses a fixed trapping radius and uncertain excimer rate.
The 1-mm RF/source mesh does not resolve the 20-um seed or certify shock resolution. Gas-flow refinement and a Sod shock benchmark are reported separately.
The seed normalization, Gaunt factor and uncertain CR coefficients are assumptions, not fitted experimental conclusions.</div>
</main><script id="payload" type="text/plain">__DATA__</script><script>
let D,C=0,S='rf',I=0,P=0,PI=0,playing=false,picplaying=false,last=0;
const el=id=>document.getElementById(id),max=a=>a.reduce((x,y)=>Math.max(x,y),-Infinity);
function note(id,s){el(id).textContent=s}
function clear(id,text){let c=el(id),x=c.getContext('2d');x.fillStyle='#0c1120';x.fillRect(0,0,c.width,c.height);x.fillStyle='#bbcadd';x.font='15px system-ui';x.fillText(text,20,140)}
function heat(id,a,log,lo,hi,label){const c=el(id),ctx=c.getContext('2d'),nx=a.length,nz=a[0].length;
let temp=document.createElement('canvas');temp.width=nx;temp.height=nz;let tx=temp.getContext('2d'),img=tx.createImageData(nx,nz);
const rgb=q=>log?[255*Math.min(1,2*q),230*Math.max(0,2*q-1),120*(1-q)]:[255*q,100*(1-Math.abs(2*q-1)),255*(1-q)];
for(let i=0;i<nx;i++)for(let j=0;j<nz;j++){let value=log?Math.log10(Math.max(a[i][j],1e-40)):a[i][j],q=Math.max(0,Math.min(1,(value-lo)/(hi-lo||1))),k=4*((nz-1-j)*nx+i),color=rgb(q);
for(let m=0;m<3;m++)img.data[k+m]=color[m];img.data[k+3]=255}
tx.putImageData(img,0,0);ctx.clearRect(0,0,c.width,c.height);ctx.imageSmoothingEnabled=false;ctx.drawImage(temp,40,12,c.width-55,c.height-95);
for(let i=0;i<256;i++){ctx.fillStyle=`rgb(${rgb(i/255).join(',')})`;ctx.fillRect(40+i*(c.width-55)/256,c.height-60,(c.width-55)/256+1,9)}
ctx.fillStyle='#e0e8f6';ctx.font='12px system-ui';ctx.fillText(`${log?'log10 ':''}${lo.toPrecision(3)}`,40,c.height-65);ctx.textAlign='right';ctx.fillText(`${log?'log10 ':''}${hi.toPrecision(3)}`,c.width-15,c.height-65);ctx.textAlign='left';ctx.fillText(label,20,c.height-14)}
function radial(r,v){let n=65,a=[];for(let i=0;i<n;i++){a[i]=[];for(let j=0;j<n;j++){let rr=Math.hypot((i-32)/32,(j-32)/32)*r[r.length-1],k=r.findIndex(x=>x>=rr);a[i][j]=k<0?0:v[k]}}return a}
function lines(id,x,series,label){let c=el(id),g=c.getContext('2d'),W=c.width,H=c.height;g.clearRect(0,0,W,H);let ym=Math.max(1e-12,max(series.flatMap(s=>s.y))),colors=['#ffc65f','#7ed5ff','#9aeeaa'];
g.strokeStyle='#64718a';g.beginPath();g.moveTo(40,12);g.lineTo(40,H-40);g.lineTo(W-12,H-40);g.stroke();
series.forEach((s,k)=>{g.strokeStyle=colors[k];g.beginPath();s.y.forEach((y,i)=>{let xx=40+(W-55)*(x[i]-x[0])/(x[x.length-1]-x[0]||1),yy=H-40-(H-55)*y/ym;i?g.lineTo(xx,yy):g.moveTo(xx,yy)});g.stroke();g.fillStyle=colors[k];g.font='12px system-ui';g.fillText(s.name,55,25+16*k)});
g.fillStyle='#dde7f5';g.fillText(label,40,H-12);g.fillText(ym.toPrecision(3),0,20)}
function draw(){let c=D.cases[C],f=c.rf,g=c.gas,b=c.summary;el('time').max=(S==='rf'?f.time_s:g.time_s).length-1;I=Math.min(I,+el('time').max);el('time').value=I;
if(S==='rf'){let peak=c.scales.rf;heat('rf',f.ey_xz[I],false,-peak,peak,`x: 0–72 mm; z: 0–132 mm; colour ±${(peak/1e5).toFixed(1)} kV/cm`);
heat('ne',f.ne_xz[I],true,8,Math.log10(c.scales.ne),'ne [m^-3]; shared scales across both cases');heat('uv',f.uv_xz[I],true,0,Math.log10(c.scales.uv),'lumped modeled UV power [W/m^3]; shared scales');
let w0=b.laser_waist_m,zr=Math.PI*w0*w0/532e-9,a=[],power=f.laser_power_w[I];for(let ix=0;ix<81;ix++){a[ix]=[];let x=(ix-40)*.000375,w=w0*Math.sqrt(1+(x/zr)**2);for(let j=0;j<61;j++){let z=(j-30)*5e-6;a[ix][j]=2*power/(Math.PI*w*w)*Math.exp(-2*z*z/(w*w))}}
heat('optics',a,true,8,Math.log10(c.scales.optics),'I [W/m^2]; x: ±15 mm; signed transverse z: ±150 um');
lines('gas',f.time_s.map(x=>x*1e9),[{name:'maximum displayed bulk Tg / 300 K',y:f.tg_xz.map(a=>Math.max(300,max(a.flat()))/300)}],'RF window time [ns]; bulk gas is isochoric here');
note('clock',`RF window: ${(f.time_s[I]*1e9).toFixed(2)} ns; laser at 2 ns`);note('rfnote','Actual Yee-field snapshots, no imposed 3.9× field multiplier.');
note('nenote',`Fluid ne is not a PIC density map. Maps average ~2 mm and sample ${(1e12*(f.time_s[1]-f.time_s[0])).toFixed(1)} ps; dense focus probe is exported. Sub-mm structure unresolved.`);
note('optnote',`Prescribed Gaussian, peak I ${(c.scales.optics).toExponential(3)} W/m², zR ${(zr*1e3).toFixed(2)} mm. Bulk-grid IB depth ${b.optical_depth_during_pulse.toExponential(2)} is NOT a physical beam-attenuation/refraction validation.`);
note('uvnote',`CR UV/excimer diagnostic; conditional visible continuum maximum ${max(f.visible_xz[I].flat()).toExponential(2)} W/m³.`);
note('gasnote',`Displayed Tg max: ${max(f.tg_xz[I].flat()).toFixed(1)} K. Gas-motion calculation is in the afterglow stage.`);
note('stats',`N_e = ${f.electrons[I].toExponential(3)}; seed = ${f.seed_electrons[I].toFixed(1)} electrons; RF stored energy = ${f.em_energy_j[I].toExponential(3)} J.`)
}else{clear('rf','RF is not solved in this local afterglow');clear('optics','Laser off; no optical propagation solve');
heat('ne',radial(g.r_m,g.electron_density[I]),true,8,Math.log10(c.scales.ngas),'axisymmetric reconstruction: radial domain 0–5 mm');
heat('uv',radial(g.r_m,g.uv_power_density[I]),true,0,Math.log10(c.scales.uvgas),'local lumped UV power [W/m^3], log10');
let p0=g.pressure[0][g.r_m.length-1],rho0=g.rho[0][g.r_m.length-1],cs=Math.sqrt(5/3*p0/rho0);
lines('gas',g.r_m.map(x=>x*1000),[{name:'p / ambient p',y:g.pressure[I].map(x=>x/p0)},{name:'rho / initial rho',y:g.rho[I].map(x=>x/rho0)},{name:'|u| / ambient sound speed',y:g.velocity[I].map(x=>Math.abs(x)/cs)}],'radius [mm]; ratios');
note('clock',`Local RF-off afterglow: ${(g.time_s[I]*1e6).toFixed(3)} us`);note('rfnote','This is a separate assumed-quench scenario, not continuation of the driven cavity.');
note('nenote','CR species advect with the gas; no post-handoff drift diffusion or sheath model.');note('optnote','No artificial optical field is animated after the pulse.');
note('uvnote',`Modeled radiation ${g.radiative_energy_j_per_m[I].toExponential(3)} J/m; total sink ${g.sink_energy_j_per_m[I].toExponential(3)} J/m includes unresolved disposal.`);
note('gasnote',`Tg max ${max(g.gas_temperature[I]).toFixed(1)} K; |u| max ${max(g.velocity[I].map(Math.abs)).toFixed(1)} m/s. Velocity peaks are NOT mesh-converged (160/320/640-cell check). A compression front is not automatically a shock.`);
if(!c.gas_summary.plasma_energy_signal_resolved)note('gasnote','Gas-heating signal is below the conservative floating-point resolution bound against ambient energy; tiny computed motion must not be interpreted as a measured expansion or shock.');
note('stats',`Handoff at ${(c.gas_summary.handoff_time_s*1e9).toFixed(2)} ns; local mass residual ${c.gas_summary.mass_residual_rel.toExponential(2)}, energy residual ${c.gas_summary.energy_residual_rel.toExponential(2)}.`)}
drawPIC()}
function drawPIC(){let p=D.pic.cases[P],frames=p.particles;el('pictime').max=frames.length-1;PI=Math.min(PI,frames.length-1);el('pictime').value=PI;
let a=frames[PI],c=el('pic'),g=c.getContext('2d'),emax=p.scales.energy,zmin=p.scales.zmin,zmax=p.scales.zmax;g.clearRect(0,0,c.width,c.height);
g.fillStyle='#7fd4ff';for(let i=0;i<a.z_m.length;i++){let x=40+(c.width-60)*(a.z_m[i]-zmin)/(zmax-zmin||1e-9),y=c.height-40-(c.height-60)*a.energy_ev[i]/emax;g.beginPath();g.arc(x,y,2,0,6.3);g.fill()}
g.fillStyle='#e7edf7';g.font='12px system-ui';g.fillText(`position ${(zmin*1e6).toFixed(1)}–${(zmax*1e6).toFixed(1)} um; energy 0–${emax.toFixed(1)} eV`,15,c.height-12);
let k=p.time_s.reduce((best,t,i)=>Math.abs(t-a.time_s)<Math.abs(p.time_s[best]-a.time_s)?i:best,0);
let bins=p.energy_bins_ev,centres=bins.slice(1).map((x,i)=>.5*(x+bins[i])),pdf=p.eedf_probability[k].map((x,i)=>x/(bins[i+1]-bins[i]));
lines('eedf',centres,[{name:'PIC electron energy probability density',y:pdf}],'energy [eV]; probability density [1/eV]');
note('picnote',`PIC local clock ${(a.time_s*1e9).toFixed(3)} ns. ${a.z_m.length} surviving REAL macroparticles in a fixed-ID sample; weights are counts per m², not 3D totals. Mean energy ${p.mean_energy_ev[k].toFixed(3)} eV. Separate imposed-field swarm, not full-cavity PIC.`);
note('picwarning',`Validation limits: ${(p.validation.validation_warnings||[]).join('; ')}. Mean-energy agreement does not validate drift, exact MCC work balance or the whole cavity.`)}
async function init(){let bytes=Uint8Array.from(atob(el('payload').textContent.trim()),x=>x.charCodeAt(0));D=await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).json();
let scales={rf:1,ne:1e9,uv:10,ngas:1e9,uvgas:10,optics:1e9};
D.cases.forEach(c=>{scales.rf=Math.max(scales.rf,max(c.rf.ey_xz.flat(2).map(Math.abs)));scales.ne=Math.max(scales.ne,max(c.rf.ne_xz.flat(2)));scales.uv=Math.max(scales.uv,max(c.rf.uv_xz.flat(2)));scales.ngas=Math.max(scales.ngas,max(c.gas.electron_density.flat()));scales.uvgas=Math.max(scales.uvgas,max(c.gas.uv_power_density.flat()));let b=c.summary;scales.optics=Math.max(scales.optics,2*b.laser_energy_j*Math.sqrt(4*Math.log(2)/Math.PI)/(b.laser_fwhm_s*Math.PI*b.laser_waist_m**2))});
D.cases.forEach((c,i)=>{c.scales=scales;el('case').add(new Option(`${c.summary.field_peak_kv_cm} kV/cm empty-cavity focus peak, 200 psig He`,i))});
D.pic.cases.forEach((c,i)=>{c.scales={energy:Math.max(40,max(c.particles.flatMap(a=>a.energy_ev))),zmin:-max(c.particles.flatMap(a=>a.z_m.map(x=>-x))),zmax:max(c.particles.flatMap(a=>a.z_m))};el('piccase').add(new Option(c.label||`PIC ${i+1}`,i))});
el('case').onchange=e=>{C=+e.target.value;I=0;draw()};el('stage').onchange=e=>{S=e.target.value;I=0;draw()};el('time').oninput=e=>{I=+e.target.value;draw()};
el('piccase').onchange=e=>{P=+e.target.value;PI=0;drawPIC()};el('pictime').oninput=e=>{PI=+e.target.value;drawPIC()};
el('play').onclick=()=>{playing=!playing;el('play').textContent=playing?'Pause':'Play'};el('picplay').onclick=()=>{picplaying=!picplaying;el('picplay').textContent=picplaying?'Pause PIC':'Play PIC'};
draw();requestAnimationFrame(tick)}
function tick(t){if(t-last>90){if(playing){I=(I+1)%(+el('time').max+1);draw()}if(picplaying){PI=(PI+1)%(+el('pictime').max+1);drawPIC()}last=t}requestAnimationFrame(tick)}
init().catch(e=>{note('clock',`Viewer error: ${e.message}`);console.error(e)});
</script></html>"""


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT / "results" / "dynamics")
    parser.add_argument("--pic", type=Path, required=True)
    args = parser.parse_args()
    build(args.root, args.pic)
