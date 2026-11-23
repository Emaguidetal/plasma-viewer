"""Parallel builder for the e-e Boltzmann tables (rf_plasma.swarm_ee) for several gases at one gas density.

python -m rf_plasma.build_ee_tables --gases Ar N2 --n-gas 2.414e26 --workers 14
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed


def _task(args):
    gas, n_gas, x = args
    from .swarm import write_table
    from .swarm_ee import ee_table, table_path
    p = table_path(gas, n_gas, x)
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        write_table(p, ee_table(gas, n_gas, float(x)), gas, f"boltzmann + e-e, x = ne/N = {x:g}, N = {n_gas:.4e}")
    return str(p)


def main(argv=None):
    from .swarm_ee import X_LEVELS
    ap = argparse.ArgumentParser()
    ap.add_argument("--gases", nargs="+", default=["Ar", "N2"])
    ap.add_argument("--n-gas", type=float, default=2.414e26)
    ap.add_argument("--workers", type=int, default=14)
    a = ap.parse_args(argv)
    tasks = [(g, a.n_gas, float(x)) for g in a.gases for x in X_LEVELS]
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(_task, t) for t in tasks]):
            print(f.result(), flush=True)


if __name__ == "__main__":
    main()
