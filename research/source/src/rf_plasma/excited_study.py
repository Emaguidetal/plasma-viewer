"""Mechanism study for the excited-state avalanche (rf_plasma.excited_avalanche): breakdown thresholds by gas,
pressure and mechanism variant, laser-seeded (cylinder, heuristic 376-electron baseline) and RF-only (sphere, single electron).

CLI:
    python -m rf_plasma.excited_study thresholds [--workers 12] [--out results/excited_avalanche]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "excited_avalanche"
GASES = ("He", "Ar", "N2")
PSIG = (50.0, 100.0, 150.0, 200.0, 300.0)
FIELDS_KV_CM = (30.0, 50.0, 100.0, 150.0, 200.0, 300.0)

# variant name -> Case overrides (kinetics set handled separately via 'kset')
VARIANTS = {
    "direct": dict(excited=False),                       # direct ionization + Ar2+/He2+/N4+ recombination
    "direct_norec": dict(excited=False, kset="norec"),   # exact gas_compare limit (no recombination)
    "linear": dict(nonlinear=False),
    "nominal": dict(),
    "nominal_low": dict(kset="low"),
    "nominal_high": dict(kset="high"),
    "hm1": dict(hm_branch=1.0),
    "eta1e-3": dict(eta_gs=1e-3),
    "eta1e-2": dict(eta_gs=1e-2),
    "eta1e-1": dict(eta_gs=1e-1),
    "eta1e-1_hm1_high": dict(eta_gs=1e-1, hm_branch=1.0, kset="high"),
    "eta1": dict(eta_gs=1.0),                            # physical upper bound: f_geo * P_survival = 1
    "eta3e-1": dict(eta_gs=0.3),
    # breakdown declared when the peak density reaches RF screening (model validity limit)
    "direct_screen": dict(excited=False, criterion="screen"),
    "nominal_screen": dict(criterion="screen"),
    "high_screen": dict(kset="high", criterion="screen"),
    # Hypothetical prescribed dense seed (peak ionization degree x0, radius 5 um); no laser-formation model.
    "spark1e-4": dict(spark_x0=1e-4),
}
MODES = {
    # Baseline laser-line seed: heuristic He-law estimate (376 electrons), not derived from the user's laser data.
    "laser": dict(geometry="cylinder", n0=376.0, sigma0_m=4.1e-6, length_m=2.03e-3, t_flat_s=100e-9),
    # RF-only: one electron present at the start of a 5 us flat pulse (most favourable arrival time)
    "rf": dict(geometry="sphere", n0=1.0, sigma0_m=1e-6, t_flat_s=5e-6),
}


def make_case(gas, psig, mode, variant, **extra):
    from . import excited_avalanche as xa
    kw = dict(VARIANTS[variant])
    kset = kw.pop("kset", "nominal")
    x0 = kw.pop("spark_x0", None)
    kin = xa.scaled_kinetics(xa.load_kinetics()[gas], "nominal" if kset == "norec" else kset)
    if kset == "norec":
        kin["ion"]["alpha_dr"] = None
    case = xa.Case(gas=gas, psig=psig, kinetics=kin, **MODES[mode], **kw)
    if x0 is not None:
        n0 = x0 * case.n_gas * 2 * math.pi * 5e-6 ** 2 * case.length_m
        case = replace(case, n0=n0, sigma0_m=5e-6, ee="dense", n_br=max(1e13, 10 * n0))
    return replace(case, **extra) if extra else case


def _threshold_task(args):
    gas, psig, mode, variant = args
    from . import excited_avalanche as xa
    t0 = time.time()
    sw = xa.Swarm(gas)
    case = make_case(gas, psig, mode, variant)
    e = xa.threshold(case, sw)
    out = dict(gas=gas, psig=psig, mode=mode, variant=variant, thr_kV_cm=e / 1e5 if math.isfinite(e) else None,
               wall_s=time.time() - t0)
    if math.isfinite(e):
        r = xa.run_case(replace(case, e0_v_m=e * 1.003), sw)
        out.update(t_bd_ns=r["t_bd_s"] * 1e9 if r["t_bd_s"] else None, en_eff_Td=r["en_eff_Td"], te_eV=r["te_eV"],
                   counters=r["counters"], wall_fraction=r["wall_fraction"], ne_peak_m3=r["ne_peak_m3"],
                   rms_radius_um=r["rms_radius_m"] * 1e6, electron_balance_rel=r["electron_balance_rel"],
                   t_screen_ns=r["t_screen_s"] * 1e9 if r.get("t_screen_s") else None,
                   n_at_screen=r.get("n_at_screen"), n_screen_m3=r.get("n_screen_m3"))
    return out


def thresholds(out=OUT, workers=12, modes=("laser", "rf"), variants=tuple(VARIANTS), gases=GASES, psig=PSIG):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    tasks = [(g, p, m, v) for m in modes for v in variants for g in gases for p in psig]
    path = out / "thresholds.json"
    done = json.loads(path.read_text()) if path.exists() else []
    have = {(d["gas"], d["psig"], d["mode"], d["variant"]) for d in done}
    tasks = [t for t in tasks if t not in have]
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_threshold_task, t): t for t in tasks}
        for f in as_completed(futs):
            try:
                d = f.result()
            except Exception as exc:  # keep the sweep alive; record the failure
                g, p, m, v = futs[f]
                d = dict(gas=g, psig=p, mode=m, variant=v, error=repr(exc))
            done.append(d)
            path.write_text(json.dumps(done, indent=1, default=float))
            print(d["gas"], d["psig"], d["mode"], d["variant"], d.get("thr_kV_cm"), d.get("error", ""), flush=True)
    meta = dict(python=sys.version.split()[0], platform=platform.platform(), variants=VARIANTS, modes=MODES,
                command="PYTHONPATH=src python -m rf_plasma.excited_study thresholds")
    (out / "thresholds_meta.json").write_text(json.dumps(meta, indent=1))
    return done


def _delay_grid_task(args):
    gas, psig, mode, variant, e_kv, fgrid, wgrid = args
    from . import excited_avalanche as xa
    import numpy as np
    sw = xa.Swarm(gas)
    base = make_case(gas, psig, mode, variant)
    tb = np.full((len(fgrid), len(wgrid)), np.nan)
    # no breakdown even at the highest field and seed factor -> skip the grid (p_breakdown = 0)
    if not xa.run_case(replace(base, e0_v_m=e_kv * 1e5 * fgrid[-1], n0=base.n0 * wgrid[-1]), sw)["breakdown"]:
        return dict(gas=gas, psig=psig, mode=mode, variant=variant, E_kV_cm=e_kv, fgrid=list(fgrid),
                    wgrid=list(wgrid), t_bd_s=tb.tolist())
    for i, f in enumerate(fgrid):
        for j, w in enumerate(wgrid):
            r = xa.run_case(replace(base, e0_v_m=e_kv * 1e5 * f, n0=base.n0 * w), sw)
            if r["breakdown"]:
                tb[i, j] = r["t_bd_s"]
    return dict(gas=gas, psig=psig, mode=mode, variant=variant, E_kV_cm=e_kv, fgrid=list(fgrid), wgrid=list(wgrid),
                t_bd_s=tb.tolist())


def sample_delays(g: dict, n_shots=2000, seed=1, sig_p=0.02, sig_el=0.03, k_order=6.0) -> dict:
    """Monte Carlo over RF power (field ~ sqrt(P)), the 376-electron laser baseline (N0 ~ E^k) and Yule factor
    W ~ Gamma(N0)/N0
    (RF-only: N0 = 1, W ~ Exp(1)), by log-linear interpolation of the tabulated breakdown times."""
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator
    rng = np.random.default_rng(seed)
    lf, lw = np.log(g["fgrid"]), np.log(g["wgrid"])
    tb = np.asarray(g["t_bd_s"], float)
    ok = np.isfinite(tb)
    if not ok.any():
        return dict(p_breakdown=0.0, delay_median_ns=None, delay_p10_ns=None, delay_p90_ns=None,
                    jitter_p90_p10_ns=None, delay_std_ns=None, n_shots=n_shots)
    # breakdown decided by the interpolated 0/1 indicator; delays interpolated in log t from finite nodes only
    # (non-breakdown nodes filled with the largest finite delay so they never create fake long delays)
    ind = RegularGridInterpolator((lf, lw), ok.astype(float), bounds_error=False, fill_value=None)
    ltab = np.log(np.where(ok, tb, np.nanmax(tb)))
    itp = RegularGridInterpolator((lf, lw), ltab, bounds_error=False, fill_value=None)
    fs = np.sqrt(np.maximum(1 + sig_p * rng.standard_normal(n_shots), 1e-6))
    if g["mode"] == "laser":
        n0 = 376.0
        el = np.maximum(1 + sig_el * rng.standard_normal(n_shots), 1e-6)
        n_nom = n0 * el ** k_order
        k = rng.poisson(n_nom)
        w = np.where(k > 0, rng.gamma(np.maximum(k, 1), 1.0) / n0, 0.0)
    else:
        w = rng.exponential(1.0, n_shots)
    ok_w = w > 0
    pts = np.c_[np.clip(np.log(fs), lf[0], lf[-1]), np.clip(np.log(np.where(ok_w, w, 1.0)), lw[0], lw[-1])]
    t = np.exp(itp(pts))
    bd = ok_w & (ind(pts) >= 0.5)
    d = t[bd] * 1e9
    q = (lambda p: float(np.percentile(d, p))) if d.size else (lambda p: None)
    return dict(p_breakdown=float(bd.mean()), delay_median_ns=q(50), delay_p10_ns=q(10), delay_p90_ns=q(90),
                jitter_p90_p10_ns=(q(90) - q(10)) if d.size else None,
                delay_std_ns=float(d.std()) if d.size > 1 else None, n_shots=n_shots)


def delays(out=OUT, workers=14, modes=("laser", "rf"), variants=("nominal", "hm1"), gases=GASES, psig=PSIG,
           fields=FIELDS_KV_CM):
    import numpy as np
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fgrid = tuple(np.exp(np.linspace(math.log(0.97), math.log(1.03), 5)))
    grids = {"laser": tuple(np.exp(np.linspace(math.log(0.5), math.log(2.5), 5))),
             "rf": tuple(np.exp(np.linspace(math.log(0.02), math.log(5.0), 6)))}
    path = out / "delays.json"
    done = json.loads(path.read_text()) if path.exists() else []
    have = {(d["gas"], d["psig"], d["mode"], d["variant"], d["E_kV_cm"]) for d in done}
    tasks = [(g, p, m, v, e, fgrid, grids[m]) for m in modes for v in variants for g in gases for p in psig
             for e in fields if (g, p, m, v, e) not in have
             and not (v.startswith("spark") and (m == "rf" or g == "N2"))]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_delay_grid_task, t): t for t in tasks}
        for f in as_completed(futs):
            try:
                d = f.result()
                d["mc"] = sample_delays(d)
            except Exception as exc:
                g, p, m, v, e = futs[f][:5]
                d = dict(gas=g, psig=p, mode=m, variant=v, E_kV_cm=e, error=repr(exc))
            done.append(d)
            path.write_text(json.dumps(done, indent=1, default=float))
            mc = d.get("mc", {})
            print(d["gas"], d["psig"], d["mode"], d["variant"], d["E_kV_cm"], mc.get("p_breakdown"),
                  mc.get("delay_median_ns"), mc.get("jitter_p90_p10_ns"), d.get("error", ""), flush=True)
    return done


def _spark_task(args):
    gas, psig, x0, kset, drp = args
    from . import excited_avalanche as xa
    sw = xa.Swarm(gas)
    kin = xa.scaled_kinetics(xa.load_kinetics()[gas], kset)
    if drp != "nominal":
        kin["ion"]["dr_product"] = None if drp == "ground" else drp
    case0 = xa.Case(gas=gas, psig=psig, kinetics=kin, **MODES["laser"])
    sig = 5e-6
    n0 = max(x0 * case0.n_gas * 2 * math.pi * sig ** 2 * case0.length_m, 376.0)
    case = replace(case0, n0=n0, sigma0_m=sig, ee="dense", n_br=max(1e13, 10.0 * n0))
    e = xa.threshold(case, sw)
    out = dict(gas=gas, psig=psig, x0=x0, kset=kset, dr_product=drp, n0=n0,
               thr_kV_cm=e / 1e5 if math.isfinite(e) else None)
    if math.isfinite(e):
        r = xa.run_case(replace(case, e0_v_m=e * 1.003), sw)
        out.update(t_grow_ns=r["t_bd_s"] * 1e9 if r["t_bd_s"] else None, counters=r["counters"],
                   x_peak_at_end=r["x_peak"], dense_max_x=r["dense_max_x"], dense_used_filled=r["dense_used_filled"],
                   rms_radius_um=r["rms_radius_m"] * 1e6)
    return out


DR_ALT = {"Ar": ("Ar2t", "ground"), "He": ("He2", "ground")}


def spark(out=OUT, workers=14, gases=("He", "Ar"), psig=PSIG, x0s=(0.0, 1e-4)):
    """Counterfactual seed-density scan: threshold for a plasma of peak ionization degree x0 (radius 5 um, cylinder)
    to grow to N_br = max(1e13, 10 N0) electrons, with k_iz and Te from the dense-regime e-e/e-ion tables, for the
    nominal / low / high kinetics sets and alternative recombination products. The applied field is assumed to act
    inside the plasma (no depolarization: valid only for a thin line along E). x0 = 0 is the heuristic 376-electron
    baseline, not a measured laser seed."""
    out = Path(out)
    tasks = [(g, p, x, k, "nominal") for g in gases for p in psig for x in x0s for k in ("nominal", "low", "high")]
    tasks += [(g, p, x, "nominal", d) for g in gases for p in psig for x in x0s if x > 0 for d in DR_ALT[g]]
    res = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for f in as_completed([ex.submit(_spark_task, t) for t in tasks]):
            d = f.result()
            res.append(d)
            print(d["gas"], d["psig"], d["x0"], d["kset"], d["dr_product"], d["thr_kV_cm"], d.get("dense_max_x"),
                  d.get("dense_used_filled"), flush=True)
    (out / "spark.json").write_text(json.dumps(res, indent=1, default=float))
    return res


def resample(out=OUT):
    """Recompute the Monte Carlo statistics of every saved delay grid (no new solver runs)."""
    path = Path(out) / "delays.json"
    done = json.loads(path.read_text())
    for d in done:
        if "t_bd_s" in d:
            d["mc"] = sample_delays(d)
    path.write_text(json.dumps(done, indent=1, default=float))
    return done


def laser_seed_sensitivity(out=OUT, psig=200.0, seeds=(0.214, 1.0, 376.0)):
    """Recompute nominal laser-line thresholds for density-inferred, conditional-one-electron and baseline seeds."""
    from . import excited_avalanche as xa

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for seed in seeds:
        for gas in ("He", "Ar"):
            case = replace(make_case(gas, psig, "laser", "nominal"), n0=float(seed))
            e = xa.threshold(case, xa.Swarm(gas))
            rows.append(dict(gas=gas, psig=psig, seed_electrons=float(seed),
                             threshold_kV_cm=e / 1e5 if math.isfinite(e) else None))
    result = {
        "description": (
            "Nominal laser-line threshold sensitivity. 0.214 is the mean count implied by 1e6 cm^-3 "
            "over the model effective Gaussian volume, not a physical fractional electron."
        ),
        "rows": rows,
    }
    (out / "laser_seed_sensitivity.json").write_text(json.dumps(result, indent=2))
    return result


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["thresholds", "delays", "spark", "resample", "laser-seed"])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--modes", nargs="+", default=["laser", "rf"])
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS))
    ap.add_argument("--gases", nargs="+", default=list(GASES))
    ap.add_argument("--psig", nargs="+", type=float, default=list(PSIG))
    a = ap.parse_args(argv)
    if a.what == "thresholds":
        thresholds(a.out, a.workers, a.modes, a.variants, a.gases, a.psig)
    elif a.what == "spark":
        spark(a.out, a.workers)
    elif a.what == "resample":
        resample(a.out)
    elif a.what == "laser-seed":
        laser_seed_sensitivity(a.out)
    else:
        delays(a.out, a.workers, a.modes, ("nominal", "spark1e-4") if a.variants == list(VARIANTS) else a.variants,
               a.gases, a.psig)


if __name__ == "__main__":
    main()
