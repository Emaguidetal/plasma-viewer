"""Independent matched-guide 3D slab checks for the standing-wave demonstration."""

import argparse
import json
import math
from pathlib import Path
import time

import numba
import numpy as np
from scipy.constants import elementary_charge as QE, electron_mass as ME, epsilon_0 as EPS0

from .cold_rf import available_memory_bytes, estimated_memory_bytes, memory_budget_bytes
from .drude import te10_slab_scattering
from .fdtd3d import FDTD, _plane_projection
from .geometry3d import TE10Port, raised_cosine_ramp, straight_guide


def validate(output: Path, threads=2):
    numba.set_num_threads(threads)
    meta = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    f0, nu, ne = meta["frequency_hz"], meta["collision_frequency_s"], meta["case"]["density_m3"]
    rows = []
    for dx in (1e-3, 0.5e-3):
        begin = time.perf_counter()
        lay = straight_guide(dx, .072, .034, .14)
        if estimated_memory_bytes(lay.node_inside.size) > memory_budget_bytes(available_memory_bytes()):
            raise MemoryError("estimated slab-check memory exceeds the shared-host budget")
        g = lay.grid
        fd = FDTD(g, lay.node_inside, lay.pml_lo, lay.pml_hi, f0)
        wp2 = np.zeros(lay.node_inside.shape)
        nun = np.zeros_like(wp2)
        nk = round(meta["sheet_thickness_m"] / dx)
        k0 = (lay.meas_k + lay.trans_k) // 2 - nk // 2
        wp2[:, :, k0:k0 + nk] = ne * QE**2 / (EPS0 * ME)
        nun[:, :, k0:k0 + nk] = nu
        fd.set_bulk_plasma(wp2, nun)
        port = TE10Port(lay, f0)
        port.amplitude = 1.0
        port.envelope = raised_cosine_ramp(1.5e-9)
        steps = math.ceil(4.5e-9 / g.dt)
        vin, iin, vout, iout = np.zeros((4, steps))
        for n in range(steps):
            fd.step(port.inject)
            vin[n], iin[n] = _plane_projection(fd.Ey, fd.Hx, lay.meas_k, port.profile, port.profile)
            vout[n], iout[n] = _plane_projection(fd.Ey, fd.Hx, lay.trans_k, port.profile, port.profile)
        ai, bi = port.modal_amplitudes(vin, iin, g.dt)
        ao, bo = port.modal_amplitudes(vout, iout, g.dt)
        period = round(1 / (f0 * g.dt))
        sl = slice(-2 * period, None)
        pf = port.power(ai[sl]).mean()
        rr, tt = float(port.power(bi[sl]).mean() / pf), float(port.power(ao[sl]).mean() / pf)
        ar, at = te10_slab_scattering(f0, lay.feed_width_m, nk * dx, ne, nu)
        row = dict(dx_m=dx, R_fdtd=rr, T_fdtd=tt, R_analytic=abs(ar)**2, T_analytic=abs(at)**2,
                   R_absolute_error=abs(rr - abs(ar)**2), T_absolute_error=abs(tt - abs(at)**2),
                   forward_power_tail_relative_std=float(port.power(ai[sl]).std() / pf),
                   elapsed_s=time.perf_counter() - begin)
        rows.append(row)
        print(json.dumps(row), flush=True)
    passed = all(r["R_absolute_error"] < .03 and r["T_absolute_error"] < .003
                 and r["forward_power_tail_relative_std"] < .02 for r in rows)
    result = dict(passed=passed, cases=rows, R_mesh_change=abs(rows[0]["R_fdtd"] - rows[1]["R_fdtd"]),
                  note="Matched-guide slab only. This is not a full cavity mesh-convergence study.")
    (output / "slab_validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    if not passed:
        raise ArithmeticError("3D slab comparison missed its declared R/T/settling tolerances")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    p.add_argument("--threads", type=int, default=2)
    args = p.parse_args()
    validate(args.output.resolve(), args.threads)
