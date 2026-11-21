"""General scan runner for the coupled 3D solver (pressure, field, laser, seed, mesh, closure).

usage (from simulations/rf_plasma, PYTHONPATH=src):
  python results/scan_runner.py --root scan_pE --geom guide --p_psig 100 --E_kVcm 80 \
      [--laser --E_mJ 1 --fwhm_ns 1 --waist_um 20] [--physics fullee] [--seed 1e9] [--dx_mm 1] \
      [--t_end_ns 24] [--threads 2] [--gas He]
Pressure is gauge (psig): p_abs = psig * 6894.757 + 101325 Pa.
physics: base | full | fullee (full + solved e-e EEDF closure) | <...>_ff (filament free->ambipolar front) | nocr ...
"""
import argparse
import json

from rf_plasma.coupled3d import Case, run_case, RESULTS

PSI = 6894.757
ap = argparse.ArgumentParser()
ap.add_argument("--root", required=True)
ap.add_argument("--geom", default="guide")
ap.add_argument("--gas", default="He")
ap.add_argument("--p_psig", type=float, default=130.4)
ap.add_argument("--p_MPa", type=float, default=None, help="absolute pressure; overrides --p_psig")
ap.add_argument("--E_kVcm", type=float, default=80.0)
ap.add_argument("--laser", action="store_true")
ap.add_argument("--E_mJ", type=float, default=1.0)
ap.add_argument("--fwhm_ns", type=float, default=1.0)
ap.add_argument("--waist_um", type=float, default=20.0)
ap.add_argument("--wavelength_nm", type=float, default=532.0)
ap.add_argument("--t_laser_ns", type=float, default=5.0)
ap.add_argument("--t_end_ns", type=float, default=24.0)
ap.add_argument("--seed", type=float, default=1e9)
ap.add_argument("--laser_seed", type=float, default=0.0, help="peak laser-localized seed density m^-3")
ap.add_argument("--seed_k", type=float, default=6.0, help="effective order of the laser seed (0: fixed profile)")
ap.add_argument("--imp_frac", type=float, default=0.0, help="impurity (N2) mole fraction photoionized by PPT")
ap.add_argument("--dx_mm", type=float, default=1.0)
ap.add_argument("--physics", default="fullee")
ap.add_argument("--threads", type=int, default=2)
ap.add_argument("--frames", type=int, default=3)
ap.add_argument("--tag", default="")
a = ap.parse_args()

flags = dict(cr=True, transport=True, gas_heating=True, propagation=True, ee_tables=False, fil_free=False)
ph = a.physics
if ph == "base":
    flags = {k: False for k in flags}
else:
    if ph.startswith("fullee"):
        flags["ee_tables"] = True
        ph = ph[len("fullee"):].lstrip("_")
    elif ph.startswith("full"):
        ph = ph[len("full"):].lstrip("_")
    for part in filter(None, ph.split("_")):
        if part == "ff":
            flags["fil_free"] = True
            continue
        off = {"nocr": "cr", "notransport": "transport", "noheat": "gas_heating", "norefr": "propagation"}[part]
        flags[off] = False

p_abs = a.p_MPa * 1e6 if a.p_MPa is not None else a.p_psig * PSI + 101325.0
ptag = f"p{a.p_MPa:g}MPa" if a.p_MPa is not None else f"p{a.p_psig:g}psig"
tag = (f"{ptag}_E{a.E_kVcm:g}kVcm_{a.physics}"
       + (f"_L{a.E_mJ:g}mJ{a.fwhm_ns:g}ns_w{a.waist_um:g}" if a.laser else "_rfonly")
       + (f"_seed{a.seed:.0e}" if a.seed != 1e9 else "") + (f"_ls{a.laser_seed:.0e}" if a.laser_seed else "")
       + (f"_imp{a.imp_frac:.0e}" if a.imp_frac else "") + (f"_dx{a.dx_mm:g}" if a.dx_mm != 1.0 else ""))
if a.tag:  # append, never replace: runs differing only in field/pressure must not share a directory
    tag = a.tag if a.tag.startswith(tag) else f"{tag}_{a.tag}"
t_end = a.t_end_ns * 1e-9
snaps = [t * 1e-9 for t in (4, 5, 8, 12, 16, 20, 24, 30, 40) if t <= a.t_end_ns]
c = Case(geometry=a.geom, gas=a.gas, laser=a.laser, rf=True, fwhm_s=a.fwhm_ns * 1e-9, energy_j=a.E_mJ * 1e-3,
         waist_m=a.waist_um * 1e-6, wavelength_m=a.wavelength_nm * 1e-9, t_laser_s=a.t_laser_ns * 1e-9,
         t_end_s=t_end, pressure_pa=p_abs, field_target_v_m=a.E_kVcm * 1e5, seed_m3=a.seed, dx_m=a.dx_mm * 1e-3,
         laser_seed_m3=a.laser_seed, laser_seed_k=a.seed_k, impurity_frac=a.imp_frac,
         frame_every_blocks=a.frames, tag=tag, snapshot_s=snaps, laser_rmax_m=max(0.8e-3, 5 * a.waist_um * 1e-6),
         **flags)
s = run_case(c, threads=a.threads, out_root=RESULTS / a.root)
print("RESULT", c.name, s["status"], json.dumps({"gamma_end": s["rf"].get("gamma_end"),
      "fil": s["ledger_filament_J"].get("closure_residual_rel"),
      "bulk": s["ledger_bulk_J"].get("closure_residual_rel")}))
