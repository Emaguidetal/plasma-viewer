"""Bounded real WarpX 26.11 local He PIC/MCC dynamics. No cavity kinetic model.

Run: D:\\mm\\wx\\python.exe pic_dynamics.py suite
Each case is a fresh process; no dependencies are installed. Outputs stay in pic.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from datetime import datetime, timezone

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import numpy as np

HERE = Path(__file__).resolve().parent
QE, ME, C, KB, EPS0 = 1.602176634e-19, 9.1093837015e-31, 299792458., 1.380649e-23, 8.8541878128e-12
MHE = 4.002602 * 1.66053906660e-27
PRESSURE = 200 * 6894.757 + 101325


def save(path, obj):
    path.write_text(json.dumps(obj, indent=1, allow_nan=False), encoding="utf-8")


def cat(arrays):
    return np.concatenate([np.asarray(a) for a in arrays]) if arrays else np.array([])


def run_case(a):
    from pywarpx import picmi, callbacks
    from pywarpx.particle_containers import ParticleContainerWrapper, ParticleBoundaryBufferWrapper
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(a.seed)
    n = a.pressure / (KB * 300)
    steps = round(a.periods / a.frequency / a.dt) if a.steps is None else a.steps
    if steps > 700000 or a.particles > 8000:
        raise ValueError("Per-case bound exceeded")
    gap = a.gap
    z = rng.uniform(.45 * gap, .55 * gap, a.particles)
    u = rng.normal(0, np.sqrt(QE * .5 / ME), (3, a.particles))
    if a.mode == "acceleration":
        u[:] = 0
    elif a.mode == "relaxation":
        u[:] = 0
        u[2] = 1.e6
    # 1D weights are particles per unit transverse area, NOT a physical 3D count.
    weight = 1.e12 * (.1 * gap) / a.particles
    dist = picmi.ParticleListDistribution(x=[0.] * a.particles, y=[0.] * a.particles,
        z=z.tolist(), ux=u[0].tolist(), uy=u[1].tolist(), uz=u[2].tolist(),
        weight=[weight] * a.particles)
    field = f"{-a.field * gap:.16e}"
    if a.mode == "rf":
        field += f"*sin(2*pi*{a.frequency:.16e}*t)"
    grid = picmi.Cartesian1DGrid(number_of_cells=[a.cells], lower_bound=[0.], upper_bound=[gap],
        warpx_max_grid_size=a.cells, lower_boundary_conditions=["dirichlet"],
        upper_boundary_conditions=["dirichlet"], lower_boundary_conditions_particles=["absorbing"],
        upper_boundary_conditions_particles=["absorbing"], warpx_potential_lo_z="0.",
        warpx_potential_hi_z=field)
    electrons = picmi.Species(particle_type="electron", name="electrons", initial_distribution=dist,
        warpx_save_particles_at_zlo=True, warpx_save_particles_at_zhi=True)
    ions = picmi.Species(particle_type="He", name="he_ions", charge="q_e", mass=MHE,
        warpx_save_particles_at_zlo=True, warpx_save_particles_at_zhi=True)
    manifest = json.loads((HERE / "xs_he/manifest.json").read_text())
    processes = {}
    xs_meta = []
    for p in manifest["processes"]:
        path = HERE / "xs_he" / p["file"]
        data = np.loadtxt(path)
        xs_meta.append(dict(name=p["name"], threshold_ev=p["energy_eV"],
            range_ev=[float(data[0, 0]), float(data[-1, 0])],
            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        entry = {"cross_section": str(path)}
        if p["energy_eV"] is not None:
            entry["energy"] = p["energy_eV"]
        if p["name"] == "ionization":
            entry["species"] = ions
        processes[p["name"]] = entry
    if a.mode == "relaxation":
        # Deliberately artificial constant cross section ONLY for analytic relaxation test.
        path = out / "benchmark_elastic.dat"
        np.savetxt(path, [[0., 1.e-21], [1.e4, 1.e-21]])
        processes = {"elastic": {"cross_section": str(path)}}
    collisions = [] if a.mode == "acceleration" else [picmi.MCCCollisions(name="he_mcc",
        species=electrons, background_density=n, background_temperature=300.,
        background_mass=MHE, scattering_processes=processes)]
    sim = picmi.Simulation(solver=picmi.ElectrostaticSolver(grid=grid), time_step_size=a.dt,
        max_steps=steps, warpx_collisions=collisions, warpx_collisions_split_momentum_push=False,
        warpx_random_seed=a.seed, warpx_used_inputs_file=str(out / "warpx_used_inputs"),
        warpx_amrex_the_arena_init_size=16000000, verbose=0)
    sim.add_species(electrons, layout=None)
    sim.add_species(ions, layout=None)
    # C++-side species totals at high cadence: window-averaged drift and step-resolved
    # field work are computed from these, not from the 65 noisy snapshot frames.
    for kind in ("ParticleMomentum", "ParticleEnergy", "ParticleNumber"):
        sim.add_diagnostic(picmi.ReducedDiagnostic(diag_type=kind, name=kind,
            period=a.reduced_period, path=str(out / "reduced") + os.sep))
    config = dict(vars(a), steps_actual=steps, neutral_density_m3=n, temperature_k=300.,
        initial_density_m3=1.e12, initial_macro_weight_per_m2=weight, xs_source=manifest["source"],
        cross_sections=xs_meta, software={k: importlib.metadata.version(k)
            for k in ("pywarpx", "amrex", "numpy", "scipy", "matplotlib")},
        python=sys.executable)
    save(out / "config.json", config)
    start = time.perf_counter()
    sim.initialize_inputs()
    sim.initialize_warpx()
    pc = ParticleContainerWrapper("electrons")
    ip = ParticleContainerWrapper("he_ions")
    buffer = ParticleBoundaryBufferWrapper()
    edges = np.linspace(0., 80., 161)
    targets = np.unique(np.round(np.linspace(0, steps, a.frames)).astype(int))
    records = []
    particles = []
    hists = []
    tracked = None
    min_pos, max_pos = gap, 0.
    target_set = set(targets.tolist())
    def sample():
        nonlocal tracked, min_pos, max_pos
        step = sim.extension.warpx.getistep(lev=0)
        if step not in target_set:
            return
        ux, uy, uz, w, zz, ids = [cat(getattr(pc, method)()) for method in
            ("get_particle_ux", "get_particle_uy", "get_particle_uz", "get_particle_weight",
             "get_particle_z", "get_particle_id")]
        if not len(w):
            raise RuntimeError("All electrons lost")
        u2 = ux*ux + uy*uy + uz*uz
        gam = np.sqrt(1. + u2 / C**2)
        energy = u2 / (gam + 1.) * ME / QE
        vz = uz / gam
        if tracked is None:
            # Stable IDs selected once: actual trajectories, not changing synthetic points.
            tracked = np.sort(rng.choice(ids, min(200, len(ids)), replace=False))
        select = np.flatnonzero(np.isin(ids, tracked))
        select = select[np.argsort(ids[select])]
        particles.append(dict(time_s=step*a.dt, z_m=zz[select].tolist(),
            energy_ev=energy[select].tolist(), weight=w[select].tolist(),
            particle_id=ids[select].astype(int).tolist()))
        ws = float(w.sum())
        hi = np.histogram(energy, bins=edges, weights=w)[0] / ws
        hists.append(hi)
        iw = cat(ip.get_particle_weight())
        escaped = {}
        for sp in ("electrons", "he_ions"):
            ew, en, count = 0., 0., 0
            for boundary in ("z_lo", "z_hi"):
                boundary_count = buffer.get_particle_boundary_buffer_size(sp, boundary)
                count += boundary_count
                if boundary_count == 0:
                    continue  # buffer is lazily allocated only after the first escape
                bw = cat(buffer.get_particle_boundary_buffer(sp, boundary, "w", 0))
                bu = [cat(buffer.get_particle_boundary_buffer(sp, boundary, k, 0))
                      for k in ("ux", "uy", "uz")]
                b2 = sum(x*x for x in bu)
                mass = ME if sp == "electrons" else MHE
                ew += float(bw.sum())
                en += float(np.sum(bw * b2 / (np.sqrt(1+b2/C**2)+1) * mass/QE))
            escaped[sp] = (count, ew, en)
        min_pos, max_pos = min(min_pos, float(zz.min())), max(max_pos, float(zz.max()))
        # Expected MCC frequencies integrated over the actual EEDF, not fitted rates.
        rates = dict(momentum=0., ionization=0., excitation_energy=0., elastic_energy=0.)
        speed = np.sqrt(u2) / gam
        for p in manifest["processes"] if a.mode == "rf" else []:
            d = np.loadtxt(HERE / "xs_he" / p["file"])
            r = n * speed * np.interp(energy, d[:, 0], d[:, 1])
            average = float(np.sum(w*r)/ws)
            if p["name"] == "elastic":
                rates["momentum"] = average
                rates["elastic_energy"] = float(np.sum(w*r*(2*ME/MHE)*
                    (energy - 1.5*KB*300/QE))/ws)
            elif p["name"] == "ionization":
                rates["ionization"] = average
            else:
                rates["excitation_energy"] += average * p["energy_eV"]
        records.append([step*a.dt, float(np.sum(w*energy)/ws), float(np.sum(w*vz)/ws),
            ws, len(w), len(iw), float(iw.sum()), float(w[energy > edges[-1]].sum()/ws),
            float(w[energy < edges[0]].sum()/ws), *escaped["electrons"], *escaped["he_ions"],
            *rates.values(), float(energy.max()), float(np.sum(w*(energy-np.sum(w*energy)/ws)**2)/ws)])
        if len(w) > 20000 or time.perf_counter()-start > a.timeout:
            raise RuntimeError("Runtime or growing-particle bound exceeded")
    sample()
    callbacks.installafterstep(sample)
    sim.step(steps)
    r = np.asarray(records)
    events = int(r[-1, 5] + r[-1, 12])
    from scipy.stats import chi2
    exposure = float(np.trapezoid(r[:, 4], r[:, 0]))
    upper = float(.5*chi2.ppf(.95, 2*(events+1))/exposure)
    residual = r[:, 3]+r[:, 10]-r[0, 3]-r[:, 6]-r[:, 13]
    validation = dict(status="executed_not_experimentally_validated", ionization_events=events,
        ionization_event_method="created ion macros + scraped ion macros; one product ion per MCC event",
        ionization_rate_validated=False, ionization_poisson_95_upper_s=upper,
        ionization_poisson_relative_se=1/np.sqrt(events) if events else None,
        ionization_rate_event_estimate_s=events/exposure,
        exposure_macro_s=exposure, few_events=events < 100,
        electron_weight_gain=float(r[-1, 3]/r[0, 3]-1),
        initial_macro_count=int(r[0, 4]), final_macro_count=int(r[-1, 4]),
        ion_weight=float(r[-1, 6]), charge_accounting_relative_max=float(np.max(abs(residual))/r[0, 3]),
        escaped_electron_macros=int(r[-1, 9]), escaped_electron_weight=float(r[-1, 10]),
        escaped_electron_energy_ev_per_m2=float(r[-1, 11]),
        escaped_ion_macros=int(r[-1, 12]), escaped_ion_weight=float(r[-1, 13]),
        histogram_overflow_probability=r[:, 7].tolist(), histogram_underflow_probability=r[:, 8].tolist(),
        histogram_normalization_error_max=float(np.max(abs(np.sum(hists, axis=1)+r[:, 7]+r[:, 8]-1))),
        observed_z_range_m=[min_pos, max_pos], maximum_energy_ev=float(r[:, 19].max()),
        runtime_s=time.perf_counter()-start, seconds_per_step=(time.perf_counter()-start)/steps)
    if a.mode == "acceleration":
        expected_v = -QE*a.field/ME*r[:, 0]
        expected_e = .5*ME/QE*expected_v**2
        validation.update(acceleration_velocity_relative_error=float(np.max(abs(r[:, 2]-expected_v))/
            max(abs(expected_v[-1]), 1)), acceleration_energy_relative_error=float(abs(r[-1, 1]/
            expected_e[-1]-1)), acceptance_relative=.01)
        validation["benchmark_pass"] = validation["acceleration_velocity_relative_error"] < .01 and validation["acceleration_energy_relative_error"] < .01
    if a.mode == "relaxation":
        predicted = 1.e6*np.exp(-n*1.e-21*1.e6*r[:, 0])
        validation.update(relaxation_drift_max_error_over_initial=float(np.max(abs(r[:, 2]-predicted))/1.e6),
            acceptance_absolute_over_initial=.10, analytic_note="constant sigma, nearly monoenergetic heavy-target isotropic drift decay; thermal/recoil corrections small")
        validation["benchmark_pass"] = validation["relaxation_drift_max_error_over_initial"] < .10
    case = dict(label=a.label, pressure_pa=a.pressure, e_peak_v_m=a.field, f_hz=a.frequency,
        dt_s=a.dt, seed=a.seed, time_s=r[:, 0].tolist(), mean_energy_ev=r[:, 1].tolist(),
        mean_vz_m_s=r[:, 2].tolist(), relative_electron_count=(r[:, 3]/r[0, 3]).tolist(),
        energy_bins_ev=edges.tolist(), eedf_probability=np.asarray(hists).tolist(),
        particles=particles, validation=validation)
    np.savez_compressed(out / "raw.npz", records=r, energy_bins_ev=edges,
        eedf_probability=hists, particle_z_m=np.array([p["z_m"] for p in particles]),
        particle_energy_ev=np.array([p["energy_ev"] for p in particles]),
        particle_weight=np.array([p["weight"] for p in particles]),
        particle_id=np.array([p["particle_id"] for p in particles]))
    save(out / "case.json", case)
    save(out / "meta.json", validation)
    print(json.dumps(dict(label=a.label, **validation)), flush=True)
    sim.finalize()


def suite(a):
    root = HERE / "results/dynamics_validation"
    root.mkdir(parents=True, exist_ok=True)
    run = root / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + f"_{os.getpid()}")
    run.mkdir()
    status = dict(status="running", run_directory=str(run), completed=[], errors=[])
    save(run / "status.json", status)
    plans = [
        ("acceleration", "acceleration", 1.e5, a.dt, a.particles, a.seed, 1000),
        ("relaxation", "relaxation", 0., a.dt, a.particles, a.seed, 2000),
        ("timing_pilot", "rf", 1.e7, a.dt, a.particles, a.seed, 10000),
        ("50kVcm", "rf", 5.e6, a.dt, a.particles, a.seed, None),
        ("100kVcm", "rf", 1.e7, a.dt, a.particles, a.seed, None),
        ("100kVcm_half_dt", "rf", 1.e7, a.dt/2, a.particles, a.seed, None),
        ("100kVcm_double_particles", "rf", 1.e7, a.dt, 2*a.particles, a.seed, None),
        ("100kVcm_seed2", "rf", 1.e7, a.dt, a.particles, a.seed+1, None)]
    planned_steps = sum(s if s else round(a.periods/a.frequency/dt) for _, _, _, dt, _, _, s in plans)
    if planned_steps > a.step_budget:
        raise ValueError(f"Planned {planned_steps} exceeds total budget {a.step_budget}")
    estimate = dict(total_steps=planned_steps, step_budget=a.step_budget, wall_budget_s=a.timeout,
        initial_macros=a.particles, refined_macros=2*a.particles, cells=a.cells,
        particle_ceiling=20000, memory_estimate_mb=200, threads=1, processes_simultaneous=1,
        previous_reference="existing 6399-macro/64-cell case: 1394 s / 358333 steps; pilot measures this driver")
    save(run / "resource_estimate.json", estimate)
    started = time.perf_counter()
    for label, mode, field, dt, count, seed, steps in plans:
        cmd = [sys.executable, str(Path(__file__).resolve()), "case", "--out", str(run/label),
            "--label", label, "--mode", mode, "--field", str(field), "--dt", str(dt),
            "--particles", str(count), "--seed", str(seed), "--pressure", str(a.pressure),
            "--frequency", str(a.frequency), "--periods", str(a.periods), "--frames", str(a.frames),
            "--cells", str(a.cells), "--gap", str(a.gap), "--timeout", str(min(2400, a.timeout)),
            "--reduced-period", str(a.reduced_period)]
        if steps:
            cmd += ["--steps", str(steps)]
        print("EXECUTE", subprocess.list2cmdline(cmd), flush=True)
        try:
            remaining = a.timeout - (time.perf_counter()-started)
            if remaining <= 0:
                raise RuntimeError("Suite wall budget exhausted")
            subprocess.run(cmd, check=True, timeout=min(2400, remaining))
            status["completed"].append(label)
            if label == "timing_pilot":
                pilot = json.loads((run/label/"meta.json").read_text())
                estimate["pilot_seconds_per_step"] = pilot["seconds_per_step"]
                estimate["projected_suite_s"] = pilot["seconds_per_step"] * planned_steps * 1.6
                save(run / "resource_estimate.json", estimate)
                if estimate["projected_suite_s"] > a.timeout:
                    raise RuntimeError("Measured timing exceeds bounded wall budget; reduce workload explicitly")
            if mode != "rf" and not json.loads((run/label/"meta.json").read_text())["benchmark_pass"]:
                raise RuntimeError("Analytic benchmark failed; targets not launched")
        except Exception as exc:
            status["errors"].append(dict(label=label, error=repr(exc)))
            break
        finally:
            save(run / "status.json", status)
    status["status"] = "bounded_runs_complete" if not status["errors"] else "blocked"
    status["wall_s"] = time.perf_counter()-started
    save(run / "status.json", status)
    # Use already installed fluid environment, read-only imports; do not install numba in wx.
    cmd = [str(HERE.parent/"rf_plasma/.venv/Scripts/python.exe"),
           str(HERE/"pic_dynamics_report.py"), str(run)]
    try:
        subprocess.run(cmd, check=True, timeout=120)
    except Exception as exc:
        status["errors"].append(dict(label="report", error=repr(exc)))
        status["status"] = "report_blocked"
        save(run / "status.json", status)
        raise
    print("OUTPUT", run, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["suite", "case"])
    p.add_argument("--out")
    p.add_argument("--label", default="local")
    p.add_argument("--mode", choices=["rf", "acceleration", "relaxation"], default="rf")
    p.add_argument("--field", type=float, default=1.e7)
    p.add_argument("--pressure", type=float, default=PRESSURE)
    p.add_argument("--frequency", type=float, default=2.95e9)
    p.add_argument("--periods", type=float, default=2.)
    p.add_argument("--dt", type=float, default=2.e-15)
    p.add_argument("--particles", type=int, default=2000)
    p.add_argument("--seed", type=int, default=271828)
    p.add_argument("--cells", type=int, default=32)
    p.add_argument("--gap", type=float, default=1.e-3)
    p.add_argument("--frames", type=int, default=65)
    p.add_argument("--steps", type=int)
    p.add_argument("--reduced-period", type=int, default=20)
    p.add_argument("--step-budget", type=int, default=2100000)
    p.add_argument("--timeout", type=float, default=7200)
    a = p.parse_args()
    if a.dt <= 0 or a.pressure <= 0 or a.frequency <= 0 or not 40 <= a.frames <= 80 or a.particles < 100:
        p.error("Invalid physical/diagnostic settings")
    (suite if a.action == "suite" else run_case)(a)


if __name__ == "__main__":
    main()
