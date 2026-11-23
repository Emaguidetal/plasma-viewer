"""Dense-regime electron kinetics: ionization rate vs E/N and ionization degree x for He / Ar (/ N2), from the
two-term Boltzmann equation with e-e collisions (Hagelaar & Pitchford 2005) and optionally e-ion momentum
transfer (Coulomb), rf_plasma.swarm_ee. Non-converged points are reported as missing (no filling).

python -m rf_plasma.dense_regime [--workers 14]
"""
from __future__ import annotations

import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "excited_avalanche" / "dense_regime.json"
EN_TD = (50.0, 30.0, 20.0, 15.0, 10.0, 8.0, 6.0, 5.0, 4.0, 3.0)     # high -> low (continuation)
X = (0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
N_GAS = 2.414e26


def _task(args):
    gas, x, ei = args
    from bolos import solver as bsolver
    from .swarm import _row
    from .swarm_ee import EESolver
    s = EESolver(gas, N_GAS, ei=ei)
    f0 = s.b.maxwell(2.0)
    rows = []
    for en in EN_TD:
        try:
            f1 = s.solve(en, x, f0)
            s.b.EN = en * bsolver.TOWNSEND
            r = _row(s.b, f1, en)
            rows.append(dict(E_N_Td=en, mean_energy_eV=r[0], muN=r[2], kiz=r[4], converged=True))
            f0 = f1
        except Exception as exc:                      # ConvergenceError or numerical failure
            rows.append(dict(E_N_Td=en, converged=False, error=type(exc).__name__))
    return dict(gas=gas, x=x, ei=ei, rows=rows)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--gases", nargs="+", default=["He", "Ar"])
    a = ap.parse_args(argv)
    tasks = [(g, x, ei) for g in a.gases for x in X for ei in ((False,) if x == 0 else (False, True))]
    res = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(_task, t) for t in tasks]):
            d = f.result()
            res.append(d)
            print(d["gas"], d["x"], d["ei"], [(r["E_N_Td"], f"{r['kiz']:.2e}" if r["converged"] else "--")
                                             for r in d["rows"]], flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(dict(n_gas=N_GAS, en_td=EN_TD, x=X, results=res,
                                   method="two-term Boltzmann (BOLOS) + e-e (Hagelaar & Pitchford 2005 sec 2.6) "
                                          "+ optional Coulomb e-ion momentum transfer; LXCat IST-Lisbon sets"),
                              indent=1))


if __name__ == "__main__":
    main()
