import sys
from rf_plasma.coupled3d import Case, run_case, RESULTS
# usage: geom kind field_MVm threads [physics] [t_end_ns]
#   physics: base (original kernels) | full (CR + transport + gas heating + propagation)
#            | nocr / notransport / noheat / norefr (full minus one stage, ablations)
geom, kind = sys.argv[1], sys.argv[2]
pulses = {"rfonly": (False, 500e-12, 1e-3), "L1mJ500ps": (True, 500e-12, 1e-3), "L10mJ3ns": (True, 3e-9, 10e-3)}
laser, fw, e = pulses[kind]
field = float(sys.argv[3]) if len(sys.argv) > 3 else 8.0
threads = int(sys.argv[4]) if len(sys.argv) > 4 else 4
physics = sys.argv[5] if len(sys.argv) > 5 else "base"
t_end = float(sys.argv[6]) * 1e-9 if len(sys.argv) > 6 else 24e-9
flags = dict(cr=True, transport=True, gas_heating=True, propagation=True)
if physics == "base":
    flags = {k: False for k in flags}
elif physics != "full":
    off = {"nocr": "cr", "notransport": "transport", "noheat": "gas_heating", "norefr": "propagation"}[physics]
    flags[off] = False
tag = f"E{field:g}MVm" + ("" if physics == "base" else f"_{physics}")
root = "laser_trigger_He_1MPa" if physics == "base" else "full_physics_He_1MPa"
if t_end < 24e-9:
    root += "_smoke"
c = Case(geometry=geom, gas="He", laser=laser, rf=True, fwhm_s=fw, energy_j=e, t_laser_s=5e-9, t_end_s=t_end,
         pressure_pa=1e6, field_target_v_m=field * 1e6, frame_every_blocks=3, tag=tag,
         snapshot_s=[t for t in (4e-9, 5e-9, 8e-9, 12e-9, 16e-9, 20e-9, 24e-9) if t <= t_end], **flags)
s = run_case(c, threads=threads, out_root=RESULTS / root)
print("RESULT", c.name, s["status"], s["peaks"], s["laser"], s["rf"].get("gamma_end"),
      s["ledger_filament_J"].get("closure_residual_rel"), s["ledger_bulk_J"].get("closure_residual_rel"))
