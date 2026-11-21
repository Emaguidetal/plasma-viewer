"""Jitter sensitivity of the laser-triggered delay (Technion He anchor, 0.2 MPa, 22.3 / 16.6 kV/cm).

Scans the RF-power rms sigma_P and the number of laser seed electrons N_L; reports mean delay and rms jitter
against the measured +-6 ns (high power) / +-9 ns (low power).
usage (simulations/rf_plasma, PYTHONPATH=src): python results/stochastic_jitter_sensitivity.py
"""
import json
from pathlib import Path

import numpy as np

from rf_plasma.stochastic_seed import Rates, Shot, simulate, stats

OUT = Path(__file__).resolve().parent / "stochastic_seed"
SIG_P = [0.02, 0.05, 0.1, 0.2]
N_L = [1e10, 1e8, 1e6, 1e4, 1e2, 10.0]
MEASURED = {"22.3": dict(jitter_ns=6.0), "16.6": dict(jitter_ns=9.0), "delay_difference_ns": [30.0, 15.0]}


def main(n_shots=800):
    r = Rates("He")
    res = {"sigma_P": SIG_P, "n_laser": N_L, "measured": MEASURED, "fields": {}}
    for ek in (22.3, 16.6):
        mean = np.full((len(SIG_P), len(N_L)), np.nan)
        std = np.full_like(mean, np.nan)
        pb = np.full_like(mean, np.nan)
        for i, sp in enumerate(SIG_P):
            for j, nl in enumerate(N_L):
                sh = Shot(e0_v_m=ek * 1e5, sig_p=sp, n_laser=nl)
                st = stats(simulate(sh, n_shots, seed=11, rates=r)[0], sh.t_laser_s)
                mean[i, j] = (st["delay_mean_s"] or np.nan) * 1e9
                std[i, j] = (st["delay_std_s"] or np.nan) * 1e9
                pb[i, j] = st["p_breakdown"]
        res["fields"][f"{ek}"] = dict(delay_mean_ns=mean.tolist(), jitter_ns=std.tolist(), p_breakdown=pb.tolist())
        print(ek, "jitter ns\n", np.round(std, 2), "\nmean ns\n", np.round(mean, 1), flush=True)
    m22 = np.array(res["fields"]["22.3"]["jitter_ns"])
    m16 = np.array(res["fields"]["16.6"]["jitter_ns"])
    d = np.array(res["fields"]["16.6"]["delay_mean_ns"]) - np.array(res["fields"]["22.3"]["delay_mean_ns"])
    cost = ((m22 - 6) / 2) ** 2 + ((m16 - 9) / 3) ** 2 + ((d - 30) / 15) ** 2
    i, j = np.unravel_index(np.nanargmin(cost), cost.shape)
    res["best"] = dict(sigma_P=SIG_P[i], n_laser=N_L[j], jitter_22=m22[i, j], jitter_16=m16[i, j],
                       delay_difference_ns=d[i, j], cost=cost[i, j])
    res["delay_difference_ns"] = d.tolist()
    OUT.mkdir(exist_ok=True)
    (OUT / "jitter_sensitivity.json").write_text(json.dumps(res, indent=1, default=float))
    print("best", res["best"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for a, (lab, z, tgt) in zip(ax, (("22.3 kV/cm jitter (ns), measured 6", m22, 6),
                                      ("16.6 kV/cm jitter (ns), measured 9", m16, 9),
                                      ("delay difference (ns), measured 30 +- 15", d, 30))):
        im = a.imshow(z, origin="lower", aspect="auto", cmap="viridis")
        a.set_xticks(range(len(N_L)), [f"{v:.0e}" for v in N_L])
        a.set_yticks(range(len(SIG_P)), [f"{v:g}" for v in SIG_P])
        a.set_xlabel("laser seed electrons N_L")
        a.set_ylabel("RF power rms sigma_P")
        a.set_title(lab, fontsize=9)
        for (ii, jj), v in np.ndenumerate(z):
            a.text(jj, ii, f"{v:.1f}", ha="center", va="center", fontsize=7,
                   color="w" if abs(v - tgt) > 0.3 * tgt else "r")
        fig.colorbar(im, ax=a)
    fig.suptitle("He 0.2 MPa laser-triggered breakdown: stochastic surrogate vs Technion measurements")
    fig.tight_layout()
    fig.savefig(OUT / "jitter_sensitivity.png", dpi=130)


if __name__ == "__main__":
    main()
