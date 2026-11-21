"""Read-only BOLOS comparison, contract checks, quantitative report and plots."""
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "rf_plasma/src"))
from rf_plasma.gas_model import load_gas_model
from rf_plasma.plasma3d import P_IP
from pic_dynamics import QE, ME, MHE, KB, EPS0, C, save


WINDOWS_PER_CYCLE = 32


def load_reduced(case_dir):
    """High-cadence WarpX species totals (electrons); None for legacy runs without them."""
    d = Path(case_dir)/"reduced"
    if not (d/"ParticleMomentum.txt").exists():
        return None
    mom, ene, num = (np.loadtxt(d/f"{k}.txt") for k in
                     ("ParticleMomentum", "ParticleEnergy", "ParticleNumber"))
    k = min(len(mom), len(ene), len(num))
    assert np.array_equal(mom[:k, 0], ene[:k, 0]) and np.array_equal(mom[:k, 0], num[:k, 0])
    w = num[:k, 6]
    return dict(t=mom[:k, 1], weight=w, pz=mom[:k, 7], kinetic_ev=ene[:k, 3]/QE,
                vz=mom[:k, 7]/(ME*w), energy_ev=ene[:k, 3]/(QE*w))


def window_average(t, y, freq, t_end):
    """Mean of y in fixed RF-phase windows (1/32 period); returns centres, means, SE."""
    width = 1/(freq*WINDOWS_PER_CYCLE)
    edges = np.arange(0., t_end+.5*width, width)
    index = np.digitize(t, edges)-1
    centres, means, errors = [], [], []
    for i in range(len(edges)-1):
        s = y[index == i]
        if len(s) >= 4:
            centres.append(.5*(edges[i]+edges[i+1]))
            means.append(s.mean())
            errors.append(s.std(ddof=1)/np.sqrt(len(s)))
    return np.array(centres), np.array(means), np.array(errors)


def nrmse(x, reference):
    return float(np.sqrt(np.mean((x-reference)**2))/max(np.sqrt(np.mean(reference**2)), 1.e-300))


def validate_contract(doc):
    assert set(doc) == {"provenance", "limits", "cases"}
    for c in doc["cases"]:
        assert all(k in c for k in ("label", "pressure_pa", "e_peak_v_m", "f_hz", "dt_s",
            "seed", "time_s", "mean_energy_ev", "mean_vz_m_s", "relative_electron_count",
            "energy_bins_ev", "eedf_probability", "particles", "validation"))
        t = np.array(c["time_s"])
        assert 40 <= len(t) <= 80 and t[0] == 0 and np.all(np.diff(t) > 0)
        for k in ("mean_energy_ev", "mean_vz_m_s", "relative_electron_count"):
            assert len(c[k]) == len(t) and np.all(np.isfinite(c[k]))
        h = np.array(c["eedf_probability"])
        assert h.shape == (len(t), len(c["energy_bins_ev"])-1) and np.all(h >= 0)
        v = c["validation"]
        assert np.max(abs(h.sum(axis=1)+np.array(v["histogram_overflow_probability"])+
            np.array(v["histogram_underflow_probability"])-1)) < 1.e-10
        assert len(c["particles"]) == len(t)
        for i, p in enumerate(c["particles"]):
            assert p["time_s"] == t[i]
            assert len(p["z_m"]) == len(p["energy_ev"]) == len(p["weight"]) <= 200
            assert len(p["particle_id"]) == len(p["z_m"])
            assert np.all(np.array(p["weight"]) > 0)
    json.dumps(doc, allow_nan=False)


def main():
    run = Path(sys.argv[1])
    gm = load_gas_model("He", 800.e-9)
    eg = gm.eps_grid
    ip = float(gm.prm[P_IP])
    def rates(e):
        return np.array([np.interp(np.log(np.maximum(e, eg[0])), np.log(eg), row)
                         for row in gm.tabs_bulk])
    cases = []
    for path in run.glob("*/case.json"):
        if path.parent.name in ("acceleration", "relaxation", "timing_pilot"):
            continue
        c = json.loads(path.read_text())
        cfg = json.loads((path.parent/"config.json").read_text())
        r = np.load(path.parent/"raw.npz")["records"]
        raw = np.load(path.parent/"raw.npz")
        assert np.array_equal(raw["records"][:, 0], c["time_s"])
        assert np.array_equal(raw["records"][:, 1], c["mean_energy_ev"])
        assert np.array_equal(raw["records"][:, 2], c["mean_vz_m_s"])
        assert np.array_equal(raw["eedf_probability"], c["eedf_probability"])
        initial_ids = set(c["particles"][0]["particle_id"])
        for i, frame in enumerate(c["particles"]):
            assert frame["particle_id"] == sorted(set(frame["particle_id"]))
            assert set(frame["particle_id"]).issubset(initial_ids)
            for field, key in (("z_m", "particle_z_m"), ("energy_ev", "particle_energy_ev"),
                               ("weight", "particle_weight"), ("particle_id", "particle_id")):
                assert np.array_equal(raw[key][i], frame[field])
        t, n, e0, freq = r[:, 0], cfg["neutral_density_m3"], c["e_peak_v_m"], c["f_hz"]
        def rhs(tt, y):
            v, e, logcount = y
            e = max(e, 1.e-3)
            km, kiz, kexc, kel = rates(e)
            field = e0*np.sin(2*np.pi*freq*tt)
            return [-QE*field/ME-n*km*v,
                -field*v-n*(kexc+kel*(e-1.5*KB*300/QE)+kiz*(ip+e)), n*kiz]
        sol = solve_ivp(rhs, (0, t[-1]), [r[0, 2], r[0, 1], 0.], t_eval=t, dense_output=True,
            method="LSODA", rtol=1.e-8, atol=[1.e-3, 1.e-7, 1.e-10], max_step=2.e-13)
        if not sol.success:
            raise RuntimeError(sol.message)
        fluid = dict(time_s=t.tolist(), mean_vz_m_s=sol.y[0].tolist(),
            mean_energy_ev=sol.y[1].tolist(), relative_electron_count=np.exp(sol.y[2]).tolist(),
            momentum_frequency_s=(n*rates(sol.y[1])[0]).tolist(),
            ionization_frequency_s=(n*rates(sol.y[1])[1]).tolist(),
            excitation_loss_ev_s=(n*rates(sol.y[1])[2]).tolist())
        m = t >= t[-1]-1/freq
        avg = lambda y: float(np.trapezoid(y[m], t[m])/(t[m][-1]-t[m][0]))
        field = e0*np.sin(2*np.pi*freq*t)
        loss = r[:, 17]+r[:, 18]+ip*r[:, 16]
        expected_loss = float(np.trapezoid(r[:, 3]*loss, t))
        rd = load_reduced(path.parent)
        snapshot_drift_error = nrmse(r[m, 2], sol.y[0, m])
        if rd is None:
            # Legacy run: only noisy snapshot frames exist; keep the old diagnostic.
            energy_error = float(np.sqrt(np.mean((r[m, 1]-sol.y[1, m])**2))/avg(r[:, 1]))
            drift_error = snapshot_drift_error
            work = float(np.trapezoid(r[:, 3]*(-field*r[:, 2]), t))
            gain = r[-1, 3]*r[-1, 1]-r[0, 3]*r[0, 1]+r[-1, 11]
            window = None
        else:
            # Phase-window means of high-cadence totals: ~260 samples per window, each
            # separated by >1 momentum-relaxation time, so sampling noise drops ~16x.
            fl_t = sol.sol(np.clip(rd["t"], 0, t[-1]))
            tw, vw, vse = window_average(rd["t"], rd["vz"], freq, t[-1])
            _, ew, ese = window_average(rd["t"], rd["energy_ev"], freq, t[-1])
            _, fvw, _ = window_average(rd["t"], fl_t[0], freq, t[-1])
            _, few, _ = window_average(rd["t"], fl_t[1], freq, t[-1])
            mw = tw >= t[-1]-1/freq
            drift_error = nrmse(vw[mw], fvw[mw])
            energy_error = nrmse(ew[mw], few[mw])
            # Exact field work on electrons from step-resolved total momentum.
            work = float(np.trapezoid(-e0*np.sin(2*np.pi*freq*rd["t"])*rd["pz"]/ME, rd["t"]))
            gain = float(rd["kinetic_ev"][-1]-rd["kinetic_ev"][0]+r[-1, 11])
            window = dict(time_s=tw.tolist(), mean_vz_m_s=vw.tolist(), mean_vz_se_m_s=vse.tolist(),
                mean_energy_ev=ew.tolist(), mean_energy_se_ev=ese.tolist(),
                fluid_mean_vz_m_s=fvw.tolist(), fluid_mean_energy_ev=few.tolist(),
                reduced_samples=int(len(rd["t"])), windows_per_cycle=WINDOWS_PER_CYCLE,
                method="WarpX ParticleMomentum/ParticleEnergy/ParticleNumber totals every reduced_period steps, averaged in fixed 1/32-period RF phase windows")
        residual = (gain+expected_loss-work)/work
        v = c["validation"]
        v.update(fluid=fluid, fluid_energy_cycle_mean_ev=avg(sol.y[1]),
            pic_energy_cycle_mean_ev=avg(r[:, 1]), fluid_energy_nrmse_last_cycle=energy_error,
            fluid_drift_nrmse_last_cycle=drift_error,
            snapshot_drift_nrmse_last_cycle=snapshot_drift_error,
            snapshot_drift_note="single 2000-macro snapshots: per-frame SE ~ v_th/sqrt(N) is ~30% of drift amplitude; noise-dominated, not used for acceptance",
            phase_window=window,
            fluid_acceptance_energy_nrmse=.10, fluid_acceptance_drift_nrmse=.15,
            fluid_energy_accept=energy_error < .10, fluid_drift_accept=drift_error < .15,
            pic_eedf_momentum_frequency_s=r[:, 15].tolist(),
            pic_eedf_ionization_frequency_s=r[:, 16].tolist(),
            pic_eedf_excitation_loss_ev_s=r[:, 17].tolist(),
            momentum_cycle_ratio=avg(r[:, 15])/avg(n*rates(sol.y[1])[0]),
            excitation_loss_cycle_ratio=avg(r[:, 17])/max(avg(n*rates(sol.y[1])[2]), 1.e-30),
            eedf_ionization_cycle_ratio=avg(r[:, 16])/max(avg(n*rates(sol.y[1])[1]), 1.e-30),
            energy_work_field_ev_per_m2=work, energy_work_gain_ev_per_m2=gain,
            energy_work_expected_collision_loss_ev_per_m2=expected_loss,
            energy_work_residual_relative=float(residual),
            energy_work_note="(Delta K + escaped + EEDF-expected collision sinks - field work)/field work; work and K exact from step-resolved WarpX totals, sinks from frame EEDFs x tabulated sigma (inelastic thresholds, ionization potential, elastic recoil)",
            fluid_rate_comparison_applicability="bulk BOLOS local-mean-energy table, Maxwellian fallback below solved EEDF range; same data is not independent experiment",
            grid_and_domain_convergence="not refined in this step budget; collective/sheath observables not validated")
        events = v["ionization_events"]
        exposure = v["exposure_macro_s"]
        lower = float(.5*chi2.ppf(.05, 2*events)/exposure) if events else 0.
        fluid_expected_events = float(np.trapezoid(r[:, 4]*n*rates(sol.y[1])[1], t))
        v.update(ionization_poisson_95_lower_s=lower,
            multiplication_poisson_standard_error_absolute=float(np.sqrt(events)/v["initial_macro_count"]),
            fluid_expected_ionization_events=fluid_expected_events,
            fluid_ionization_event_ratio=events/fluid_expected_events if fluid_expected_events > 0 else None,
            fluid_ionization_comparison="not validated: fewer than 100 events" if events < 100 else
                "count interval test only; two-cycle transient, no low-field extrapolation",
            fluid_ionization_rate_inside_count_interval=lower <= fluid_expected_events/exposure <=
                v["ionization_poisson_95_upper_s"] if events >= 100 else None)
        elastic_xs = np.loadtxt(HERE/"xs_he/elastic.dat")
        threshold_iz = next(p["threshold_ev"] for p in cfg["cross_sections"] if p["name"] == "ionization")
        total_sigma = np.zeros(10001)
        test_energy = np.linspace(0., max(v["maximum_energy_ev"], 1.), len(total_sigma))
        for p in cfg["cross_sections"]:
            xs = np.loadtxt(HERE/"xs_he"/(p["name"]+".dat"))
            total_sigma += np.interp(test_energy, xs[:, 0], xs[:, 1])
        observed_frequency_bound = float(np.max(n*np.sqrt(2*QE*test_energy/ME)*total_sigma))
        v.update(estimated_max_collision_frequency_observed_range_s=observed_frequency_bound,
            estimated_max_nu_dt=observed_frequency_bound*c["dt_s"],
            estimated_max_collision_probability_observed_range=1-np.exp(-observed_frequency_bound*c["dt_s"]),
            cross_section_range_status="all observed laboratory electron energies below tabulated maxima" if all(
                v["maximum_energy_ev"] < p["range_ev"][1] for p in cfg["cross_sections"]) else
                "WARNING: some sampled laboratory energies exceed a table upper bound",
            initial_macro_weight_per_m2=cfg["initial_macro_weight_per_m2"],
            actual_column_density_per_m2=r[:, 3].tolist(),
            ion_column_density_per_m2=r[:, 6].tolist(),
            global_domain_average_ne_m3=(r[:, 3]/cfg["gap"]).tolist(),
            density_note="column count / 1mm = domain-averaged density; initial central seed density is 1e12 m^-3, initial domain mean 1e11",
            raw_json_consistency_pass=True,
            fixed_id_sample_consistency_pass=True,
            estimated_initial_thermal_mfp_m=float(1/(n*np.interp(.75, elastic_xs[:, 0], elastic_xs[:, 1]))),
            energy_standard_error_iid_ev=(np.sqrt(r[:, 20]/r[:, 4])).tolist(),
            energy_standard_error_note="equal-weight macro snapshot IID estimate; correlated ionization families and serial frames mean this is not an independent time-series confidence interval",
            probability_above_ionization_threshold=[float(np.sum(np.array(row)[np.array(c["energy_bins_ev"][:-1]) >=
                threshold_iz])) for row in c["eedf_probability"]])
        # A separate diagnostic of RF phase response, NOT a replacement for real histories.
        # Single-frame drift is a small difference of broad thermal velocities.
        phase = 2*np.pi*freq*t[m]
        design = np.column_stack([np.ones(len(phase)), np.sin(phase), np.cos(phase),
                                  np.sin(3*phase), np.cos(3*phase)])
        coef, _, _, _ = np.linalg.lstsq(design, r[m, 2], rcond=None)
        fcoef, _, _, _ = np.linalg.lstsq(design, sol.y[0, m], rcond=None)
        noise = r[m, 2] - design@coef
        variance = float(np.sum(noise**2)/(len(noise)-len(coef)))
        cov = variance*np.linalg.inv(design.T@design)
        response_error = float(np.sqrt(np.mean((design@(coef-fcoef))**2))/
            max(np.sqrt(np.mean((design@fcoef)**2)), 1))
        difference = coef[1:3]-fcoef[1:3]
        wald = float(difference@np.linalg.inv(cov[1:3, 1:3])@difference)
        v.update(drift_rf_harmonic_fit=dict(method="last-cycle OLS: constant, fundamental and third harmonic",
            pic_coefficients_m_s=coef.tolist(), fluid_coefficients_m_s=fcoef.tolist(),
            coefficient_standard_error_m_s=np.sqrt(np.diag(cov)).tolist(),
            residual_rms_m_s=float(np.sqrt(np.mean(noise**2))),
            response_nrmse=response_error, fundamental_wald_chi2=wald,
            fundamental_95_compatibility=wald < 5.991464547107979,
            uncertainty_note="approximate OLS residual noise, 2 coefficients; independent-frame assumption justified only by frame separation >> collision time; not a robust multi-seed CI"),
            validation_warnings=[message for condition, message in (
                (not v["fluid_drift_accept"], "phase-window drift vs local-mean-energy fluid fails 15% criterion (closure mismatch, not sampling noise)"),
                (abs(v["energy_work_residual_relative"]) > .10, "electron energy budget (field work vs Delta K + expected collision sinks) fails 10%"),
                (events < 100, "rare ionization; do not validate a low-field rate from ratio"))
                if condition])
        v["neutral_depletion_fraction_upper_using_ions"] = float((r[-1, 6]+r[-1, 13])/
            (n*cfg["gap"]))
        v["fixed_gas_heating_fraction_scale"] = float((r[:, 3]/cfg["gap"]).max()*50*QE/
            (1.5*n*KB*300))
        v["fixed_gas_heating_scale_note"] = "illustrative 50 eV deposited per live electron divided by neutral thermal energy; not solved gas temperature"
        cases.append(c)
    by = {c["label"]: c for c in cases}
    comparisons = []
    if "100kVcm" in by:
        base = by["100kVcm"]
        t = np.array(base["time_s"])
        m = t >= t[-1]-1/base["f_hz"]
        for name in ("100kVcm_half_dt", "100kVcm_double_particles", "100kVcm_seed2"):
            if name not in by:
                continue
            other = by[name]
            metrics = dict(label=name)
            bw, ow = base["validation"]["phase_window"], other["validation"]["phase_window"]
            for field, threshold in (("mean_energy_ev", .10), ("mean_vz_m_s", .15)):
                if bw and ow:
                    tb = np.array(bw["time_s"])
                    mb = tb >= t[-1]-1/base["f_hz"]
                    b = np.array(bw[field])[mb]
                    x = np.interp(tb[mb], ow["time_s"], ow[field])
                    noise = np.sqrt(np.mean(np.array(bw[field.replace("_m_s", "_se_m_s").replace("_ev", "_se_ev")])[mb]**2))
                    metrics[field+"_expected_noise_nrmse"] = float(np.sqrt(2)*noise/np.sqrt(np.mean(b**2)))
                    metrics[field+"_basis"] = "1/32-period phase windows of WarpX reduced diagnostics, last cycle"
                else:
                    b = np.array(base[field])[m]
                    x = np.interp(t, other["time_s"], other[field])[m]
                value = nrmse(x, b)
                metrics[field+"_nrmse"] = value
                metrics[field+"_threshold"] = threshold
                metrics[field+"_pass"] = value < threshold
            metrics["multiplication_absolute_difference"] = abs(other["relative_electron_count"][-1]-
                base["relative_electron_count"][-1])
            # Rare-event uncertainty only; do not confuse it with energy SE.
            metrics["ionization_counts"] = [base["validation"]["ionization_events"],
                other["validation"]["ionization_events"]]
            comparisons.append(metrics)
        base["validation"]["refinement_comparisons"] = comparisons
        if "100kVcm_seed2" in by:
            other = by["100kVcm_seed2"]
            x = np.interp(t, other["time_s"], other["mean_energy_ev"])
            b = np.array(base["mean_energy_ev"])
            base["validation"]["two_seed_cycle_energy_standard_error_ev"] = float(
                abs(np.trapezoid(x[m]-b[m], t[m]))/(2*(t[m][-1]-t[m][0])))
            base["validation"]["seed_uncertainty_note"] = "two seeds only, one degree of freedom; descriptive SE, not a robust confidence interval"
    status = json.loads((run/"status.json").read_text())
    cfgpath = next(iter(run.glob("*/config.json")), None)
    cfg = json.loads(cfgpath.read_text()) if cfgpath else {}
    ne, te_ev = 1.e12, .5
    omega = np.sqrt(ne*QE**2/(EPS0*ME))
    limits = dict(model="1D3V electrostatic PIC/MCC, dilute prescribed-field swarm; 1 mm representative local region",
        rf="E_z=E_peak sin(2*pi*f*t); uncalibrated 50/100 kV/cm peak; no fitted 3.9 multiplier",
        gas="pure He fixed thermal neutral reservoir, 300 K; no depletion or self-consistent gas heating",
        seed="illustrative Maxwellian 0.75 eV mean, ne=1e12 m^-3 in central 0.1 mm; laser seed only",
        boundaries="absorbing particles, Dirichlet plates; no emission; ions move with real He mass",
        omitted=["Coulomb collisions", "attachment", "recombination", "metastable populations",
            "superelastic/stepwise ionization", "molecular/dimer/three-body chemistry",
            "optical ionization", "RF screening/feedback", "cavity S-parameters", "discharge/sheath validation"],
        macro_weights="electron number per m^2 transverse area per macro; no physical total without specified transverse area",
        eedf_normalization="probability per energy bin, weight-normalized by ALL live electrons; density per eV = probability/bin width; overflow reported separately",
        particles="deterministic fixed initial particle IDs, up to 200; real samples only, births appear in full EEDF not tracked initial sample",
        plasma_frequency_rad_s=omega, plasma_frequency_hz=omega/(2*np.pi),
        debye_length_m=np.sqrt(EPS0*te_ev/(ne*QE)), collisionless_skin_depth_m=C/omega,
        dx_m=cfg.get("gap", .001)/cfg.get("cells", 32), omega_p_dt=omega*cfg.get("dt", 2.e-15),
        transit_time_1mm_at_1e6ms_s=1.e-9,
        solver_resolution="ES: no EM CFL constraint; dx resolves imposed constant field and dilute charge, not dense sheaths. dt resolves electron-neutral collisions; half-dt tested.",
        dense_plasma_cost="No density, c or mass scaling used. At ne=1e20 m^-3 and Te=1 eV Debye length ~0.74 um: full atmospheric cavity requires enormous cells/particles and ps plasma resolution; not attempted.",
        stochastic="rare ionization: exact one-sided 95% Poisson count upper limit / measured macro-time exposure; effective constant hazard only",
        validation_scope="analytic acceleration/relaxation plus numerical cross-code local observables; no experimental validation; grid/domain convergence remains untested")
    provenance = dict(status=status, run_directory=str(run), software=cfg.get("software"),
        python_warpx=cfg.get("python"), python_fluid=sys.executable,
        source="LoKI He / LXCat, existing local converted two-column data (eV,m^2); no raw download redistributed",
        cross_sections=cfg.get("cross_sections"), cross_section_source=cfg.get("xs_source"),
        scattering="isotropic elastic using momentum-transfer sigma; excitation removes threshold then scatters; WarpX background MCC impact ionization creates electron and He+; no custom redistribution",
        starting_example="https://warpx.readthedocs.io/en/latest/usage/examples/capacitive_discharge/README.html",
        algorithm_reference="https://warpx.readthedocs.io/en/latest/theory/multiphysics/collisions.html",
        fluid_provenance=gm.provenance if hasattr(gm, "provenance") else "existing GasModel.tabs_bulk",
        driver_sha256=hashlib.sha256((HERE/"pic_dynamics.py").read_bytes()).hexdigest())
    provenance["native_banner"] = "WarpX (Unknown), PICSAR 26.10, AMReX 26.11; distribution metadata pywarpx 26.11"
    provenance["observed_build"] = "native Windows CPU, double fields/particles, 1D Z, Boris, linear shape, energy-conserving gather, multigrid ES; OMP_NUM_THREADS=1; exact used inputs saved per case"
    provenance["report_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    provenance["resource_estimate"] = json.loads((run/"resource_estimate.json").read_text())
    commands = []
    for path in run.glob("*/config.json"):
        conf = json.loads(path.read_text())
        commands.append(dict(label=conf["label"], command=[conf["python"], str(HERE/"pic_dynamics.py"),
            "case", "--out", conf["out"], "--label", conf["label"], "--mode", conf["mode"],
            "--field", str(conf["field"]), "--pressure", str(conf["pressure"]), "--frequency",
            str(conf["frequency"]), "--dt", str(conf["dt"]), "--particles", str(conf["particles"]),
            "--seed", str(conf["seed"]), "--steps", str(conf["steps_actual"]), "--cells",
            str(conf["cells"]), "--gap", str(conf["gap"]), "--frames", str(conf["frames"]),
            "--timeout", str(conf["timeout"]), "--reduced-period", str(conf.get("reduced_period", 20))]))
    save(run/"commands.json", commands)
    doc = dict(provenance=provenance, limits=limits, cases=cases)
    validate_contract(doc)
    save(run/"pic_dynamics.json", doc)
    save(run.parent/"pic_dynamics.json", doc)
    summary = dict(status=status["status"], schema_valid=True, cases=[dict(label=c["label"],
        **{k: c["validation"][k] for k in ("ionization_events", "ionization_poisson_95_upper_s",
        "electron_weight_gain", "escaped_electron_macros", "histogram_normalization_error_max",
        "charge_accounting_relative_max", "pic_energy_cycle_mean_ev", "fluid_energy_cycle_mean_ev",
        "fluid_energy_nrmse_last_cycle", "fluid_drift_nrmse_last_cycle", "snapshot_drift_nrmse_last_cycle",
        "momentum_cycle_ratio", "energy_work_residual_relative", "runtime_s")}) for c in cases],
        convergence=comparisons, benchmarks={p.parent.name: json.loads(p.read_text())
            for p in run.glob("*/meta.json") if p.parent.name in ("acceleration", "relaxation")})
    if "relaxation" in summary["benchmarks"]:
        br = np.load(run/"relaxation/raw.npz")["records"]
        bc = json.loads((run/"relaxation/config.json").read_text())
        nu = bc["neutral_density_m3"]*1.e-21*1.e6
        thermal = 1.5*KB*300/QE
        expected = thermal+(br[0, 1]-thermal)*np.exp(-2*ME/MHE*nu*br[-1, 0])
        summary["benchmarks"]["relaxation"].update(mean_energy_initial_ev=float(br[0, 1]),
            mean_energy_final_ev=float(br[-1, 1]), approximate_elastic_cooling_expected_ev=float(expected),
            elastic_cooling_relative_error=float(abs(br[-1, 1]/expected-1)),
            elastic_cooling_note="near-constant-speed heavy-target exponential cooling, finite thermal reservoir; approximate analytic test")
    save(run/"summary.json", summary)
    for item, c in zip(summary["cases"], cases):
        item.update({key: c["validation"][key] for key in ("fluid_expected_ionization_events",
            "fluid_ionization_event_ratio", "fluid_ionization_comparison",
            "estimated_max_nu_dt", "estimated_max_collision_probability_observed_range",
            "cross_section_range_status", "initial_macro_weight_per_m2")})
    save(run/"summary.json", summary)
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    for label in ("50kVcm", "100kVcm"):
        if label not in by:
            continue
        c = by[label]; t = np.array(c["time_s"])*1.e9
        color = "C0" if label == "50kVcm" else "C1"
        fl = c["validation"]["fluid"]
        for column, key, ylabel in ((0, "mean_energy_ev", "mean energy (eV)"),
            (1, "mean_vz_m_s", "mean v_z (m/s)"), (2, "relative_electron_count", "electron weight / initial")):
            pw = c["validation"].get("phase_window")
            if key == "mean_vz_m_s" and pw:
                ax[0, 1].plot(t, c[key], ".", color=color, ms=2, alpha=.35, label=label+" snapshots")
                ax[0, 1].plot(np.array(pw["time_s"])*1.e9, pw[key], color=color, label=label+" PIC phase-window")
                ax[0, 1].plot(np.array(pw["time_s"])*1.e9, pw["fluid_"+key], "--", color=color, label=label+" fluid")
                ax[0, 1].set(xlabel="time (ns)", ylabel=ylabel)
                continue
            ax[0, column].plot(t, c[key], color=color, label=label+" PIC")
            ax[0, column].plot(t, fl[key], "--", color=color, label=label+" fluid")
            ax[0, column].set(xlabel="time (ns)", ylabel=ylabel)
        edges = np.array(c["energy_bins_ev"])
        h = np.array(c["eedf_probability"])
        ax[1, 0].plot((edges[1:]+edges[:-1])/2, h[-1]/np.diff(edges), label=label)
        pp = c["particles"]
        for j in range(min(12, len(pp[0]["z_m"]))):
            ax[1, 1].plot(t, [p["z_m"][j]*1.e3 for p in pp], color=color, alpha=.4, lw=.7)
    ax[1, 0].set(xlabel="energy (eV)", ylabel="probability density (1/eV)", title="Actual final EEDF", xlim=(0, 35))
    ax[1, 1].set(xlabel="time (ns)", ylabel="z (mm)", title="Actual fixed-ID electron samples")
    if comparisons:
        xs = np.arange(len(comparisons))
        ax[1, 2].bar(xs-.2, [x["mean_energy_ev_nrmse"] for x in comparisons], .4, label="energy (limit 0.10)")
        ax[1, 2].bar(xs+.2, [x["mean_vz_m_s_nrmse"] for x in comparisons], .4, label="drift (limit 0.15)")
        ax[1, 2].set_xticks(xs, [x["label"].replace("100kVcm_", "") for x in comparisons])
        ax[1, 2].axhline(.10, color="C0", ls="--")
        ax[1, 2].axhline(.15, color="C1", ls="--")
        ax[1, 2].legend(fontsize=8)
    ax[1, 2].set(ylabel="last-cycle NRMSE vs base", title="Numerical sensitivity (phase windows)")
    for axis in ax[0]:
        axis.legend(fontsize=8)
    ax[1, 0].legend()
    fig.suptitle("Real local WarpX PIC/MCC — He 200 psig, 300 K, 2.95 GHz; imposed RF, no cavity feedback")
    fig.tight_layout()
    fig.savefig(run/"diagnostics.png", dpi=140)
    targets = [by[label] for label in ("50kVcm", "100kVcm") if label in by]
    if targets:
        efig, eax = plt.subplots(len(targets), 1, figsize=(10, 3.5*len(targets)), squeeze=False)
        for axis, c in zip(eax[:, 0], targets):
            edges = np.array(c["energy_bins_ev"])
            density = np.array(c["eedf_probability"])/np.diff(edges)
            image = axis.pcolormesh(np.array(c["time_s"])*1.e9,
                (edges[1:]+edges[:-1])/2, density.T, shading="nearest", cmap="magma")
            axis.set(xlabel="time (ns)", ylabel="energy (eV)", ylim=(0, 35),
                title=c["label"]+" actual weight-normalized EEDF frames")
            efig.colorbar(image, ax=axis, label="probability density (1/eV)")
        efig.tight_layout()
        efig.savefig(run/"eedf_frames.png", dpi=140)
    report = "# Bounded local WarpX dynamics validation\n\n" + json.dumps(summary, indent=2) + "\n\n" + json.dumps(limits, indent=2)
    (run/"REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
