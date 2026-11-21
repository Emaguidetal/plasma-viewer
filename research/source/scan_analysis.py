"""Collect scan runs into one table (CSV + JSON): breakdown time, reflection, channel size, closures.

usage: python results/scan_analysis.py <root> [<root> ...]   (roots under results/coupled3d)
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

KB = 1.380649e-23
BASE = Path(__file__).resolve().parent / "coupled3d"


def analyse(d: Path) -> dict:
    s = json.loads((d / "summary.json").read_text())
    c = s["case"]
    ts = np.load(d / "timeseries.npz")
    rf = np.load(d / "rf_series.npz")
    t = rf["t_phasor"]
    fwd = np.abs(rf["fwd"])
    g = np.abs(rf["bwd"]) / np.where(fwd > 0, fwd, np.inf)
    g_end = float(s["rf"]["gamma_end"])

    def t_cross(level):
        i = np.flatnonzero(g > level)
        return float(t[i[0]]) if i.size else None

    tt = ts["t"]
    n_tot = ts["fil_electrons"] + ts["bulk_electrons"]
    t03 = t_cross(0.3)
    n_br = float(np.interp(t03, tt, n_tot)) if t03 is not None else None
    n_gas = c["pressure_pa"] / (KB * c["temperature_k"])
    t_l = c["t_laser_s"] if c["laser"] else 0.0
    i_end = len(tt) - 1
    r50 = ts["fil_r50_focus"]
    ne_ax = ts["fil_ne_axis_focus"]
    i_l = int(np.searchsorted(tt, t_l + 3 * c["fwhm_s"])) if c["laser"] else 0
    i_l = min(i_l, i_end)
    pk = s["peaks"]
    i_nh = np.flatnonzero(n_tot >= 0.5 * n_tot.max())
    t_nhalf = (float(tt[i_nh[0]]) - t_l) * 1e9 if n_tot.max() > 0 and i_nh.size else None
    # focus |Ey| envelope (max over each RF period) vs its value at the laser time: cavity-safe delay
    ey, dte = np.abs(rf["ey_focus"]), float(rf["dt"])
    per = max(int(round(1.0 / (2.95e9 * dte))), 1)
    nper = len(ey) // per
    t_coll = None
    if nper > 2:
        env = ey[:nper * per].reshape(nper, per).max(axis=1)
        te_ = (np.arange(nper) + 1) * per * dte
        ref = float(np.interp(t_l, te_, env))
        j = np.flatnonzero((te_ > t_l) & (env < 0.5 * ref))
        t_coll = (float(te_[j[0]]) - t_l) * 1e9 if ref > 0 and j.size else None
    return dict(
        run=d.name, root=d.parent.name, status=s["status"], geometry=c["geometry"],
        p_MPa=c["pressure_pa"] / 1e6, p_psig=(c["pressure_pa"] - 101325.0) / 6894.757,
        E_kVcm=c["field_target_v_m"] / 1e5, E_N_Td=c["field_target_v_m"] / n_gas * 1e21,
        laser=c["laser"], E_mJ=c["energy_j"] * 1e3, fwhm_ns=c["fwhm_s"] * 1e9, waist_um=c["waist_m"] * 1e6,
        seed_m3=c["seed_m3"], dx_mm=c["dx_m"] * 1e3, ee=c.get("ee_tables", False), ff=c.get("fil_free", False),
        laser_seed_m3=c.get("laser_seed_m3", 0.0), laser_seed_k=c.get("laser_seed_k", 0.0),
        seed_electrons=s["laser"].get("seed_electrons"),
        t_ch21_ns=(float(tt[np.flatnonzero(ne_ax > 1e21)[0]]) - t_l) * 1e9 if np.any(ne_ax > 1e21) else None,
        cr=c["cr"], gamma_end=g_end,
        t_br03_ns=None if t03 is None else t03 * 1e9,
        delay03_ns=None if t03 is None else (t03 - t_l) * 1e9,
        t_half_ns=t_cross(0.5 * g_end) * 1e9 if g_end > 0.05 and t_cross(0.5 * g_end) else None,
        t_Nhalf_ns=t_nhalf, t_Ecollapse_ns=t_coll,
        n_electrons_at_Ecollapse=float(np.interp(t_l + t_coll * 1e-9, tt, n_tot)) if t_coll is not None else None,
        n_electrons_at_br=n_br,
        fil_electrons_after_laser=float(ts["fil_electrons"][i_l]),
        bulk_electrons_after_laser=float(ts["bulk_electrons"][i_l]),
        fil_ne_axis_max=pk["fil_ne_axis_max"], fil_te_axis_max=pk["fil_te_axis_max"],
        fil_tg_axis_max=pk["fil_tg_axis_max"], bulk_ne_max_end=pk["bulk_ne_max_end"],
        bulk_tg_max=pk.get("bulk_tg_max_max"),
        r50_um_end=float(r50[-1] * 1e6),
        r50_um_at_br=float(np.interp(t03, tt, r50) * 1e6) if t03 is not None else None,
        ne_axis_end=float(ne_ax[-1]),
        closure_fil=s["ledger_filament_J"].get("closure_residual_rel"),
        closure_bulk=s["ledger_bulk_J"].get("closure_residual_rel"),
        wall_s=s["wall_s"])


def main(roots):
    rows = []
    for r in roots:
        for d in sorted((BASE / r).glob("*")):
            if (d / "summary.json").exists():
                try:
                    rows.append(analyse(d))
                except Exception as exc:  # keep going; report the broken run
                    print("skip", d.name, exc)
    out = BASE / ("table_" + "+".join(roots))
    if rows:
        with open(out.with_suffix(".csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        out.with_suffix(".json").write_text(json.dumps(rows, indent=1))
    return rows


if __name__ == "__main__":
    for row in main(sys.argv[1:]):
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()
               if k in ("run", "E_N_Td", "gamma_end", "delay03_ns", "r50_um_end", "fil_tg_axis_max",
                        "n_electrons_at_br")})
