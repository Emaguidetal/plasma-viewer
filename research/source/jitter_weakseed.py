"""0D Monte-Carlo delay/jitter at the cavity anchors with the weak seed that the 3D runs need
(N0 ~ 1e5, N at field collapse ~ 1e12-4e12; see coupled3d/figures_scan/scan_summary.json)."""
import json
from pathlib import Path

from rf_plasma.stochastic_seed import Rates, Shot, simulate, stats

OUT = Path(__file__).resolve().parent / "stochastic_seed" / "anchor_weakseed.json"


def main():
    r = Rates("He")
    res = {}
    for n0, nbr in ((1e5, 2e12), (8.4e5, 6e13), (3.2e10, 2.5e13)):
        for sig_p in (0.02, 0.05, 0.10):
            row = {}
            for e_kv in (16.6, 22.3):
                sh = Shot(e0_v_m=e_kv * 1e5, n_laser=n0, n_br=nbr, sig_p=sig_p)
                st = stats(simulate(sh, 2000, seed=3, rates=r)[0], sh.t_laser_s)
                row[f"{e_kv}"] = {(k[:-2] + "_ns" if k.endswith("_s") else k): (v * 1e9 if k.endswith("_s") and v is not None else v)
                                  for k, v in st.items()}
            row["delay_difference_ns"] = row["16.6"]["delay_mean_ns"] - row["22.3"]["delay_mean_ns"]
            res[f"N0={n0:.1e}_Nbr={nbr:.1e}_sigP={sig_p}"] = row
    OUT.write_text(json.dumps(res, indent=1))
    for k, v in res.items():
        print(k, f"d16.6={v['16.6']['delay_mean_ns']:.1f}±{v['16.6']['delay_std_ns']:.1f}",
              f"d22.3={v['22.3']['delay_mean_ns']:.1f}±{v['22.3']['delay_std_ns']:.1f}",
              f"diff={v['delay_difference_ns']:.1f} ns")


if __name__ == "__main__":
    main()
