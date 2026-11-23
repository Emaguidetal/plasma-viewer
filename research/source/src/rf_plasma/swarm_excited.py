"""Excitation COUNT rate coefficients per excited-state group (He, Ar, N2) from the same two-term BOLOS solution
as swarm_gases (same cross sections, grid and E/N points).

The existing swarm tables store kexc = sum(threshold * rate), an energy-loss coefficient. A reservoir model of
excited species needs the number of excitations per unit time, grouped by the state into which the excitation
cascades on the ns time scale at 5-20 bar:
  Ar: 'm'   = Ar(3P2) + Ar(3P0)                (4s metastables)
      'r'   = Ar(3P1) + Ar(1P1)                (4s resonant, radiation trapped)
      'hi'  = 4p, 3d, 5s, 5p, ...              (radiative cascade mainly to 4s within ~10-100 ns)
  He: 't'   = triplets (2^3S, 2^3P, n>=3 triplets; cascade to 2^3S)
      's'   = singlets (2^1S, 2^1P, n>=3 singlets)
  N2: 'trip' = A, B, W, B', C, E               (C -> B -> A cascade; W, B' -> B)
      'sing' = a', a, w, a'' and higher singlets
Electronic excitation only; N2 vibrational/rotational excitations are energy losses, not reservoirs.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from bolos import solver as bsolver

from . import swarm_gases as sg
from .swarm import SWARM_DIR, read_table

GROUP_COLUMNS = {"Ar": ["k_m", "k_r", "k_hi"], "He": ["k_t", "k_s"], "N2": ["k_trip", "k_sing"]}
# Hornbeck-Molnar (associative ionization X** + X -> X2+ + e) is energetically open for levels above
# IP(X) - D0(X2+): Ar 15.760 - 1.33 = 14.43 eV, He 24.587 - 2.47 = 22.12 eV (n >= 3). 'k_ai' counts excitations
# into those levels (a subset already included in the groups above, so it must not be added to them).
AI_THRESHOLD_EV = {"Ar": 14.43, "He": 22.12, "N2": math.inf}


def group_of(gas: str, product: str) -> str | None:
    p = product.replace(" ", "")
    if gas == "Ar":
        if p in ("Ar(3P2)", "Ar(3P0)"):
            return "k_m"
        if p in ("Ar(3P1)", "Ar(1P1)"):
            return "k_r"
        return "k_hi"
    if gas == "He":
        return "k_t" if p.endswith("3)") else "k_s"
    if gas == "N2":
        if "v=0-v=" in p.replace("(", "").replace(")", "") or "(v=" in p or p.startswith("N2(J") or "J=" in p:
            return None
        if any(s in p for s in ("A3Su", "B3Pg", "W3Du", "B3Su", "C3Pu", "E3Sg")):
            return "k_trip"
        if any(s in p for s in ("a1Su", "a1Pg", "w1Du", "a1Sg", "sing")):
            return "k_sing"
        return None
    raise ValueError(gas)


def table(gas: str, en_td=sg.EN_GRID_TD) -> dict:
    b = sg.make_solver(gas)
    cols = GROUP_COLUMNS[gas]
    f0 = b.maxwell(0.5)
    out = {c: [] for c in ["E_N_Td", "mean_energy_eV", "kiz_m3_s"] + cols + ["k_ai"]}
    for en in en_td:
        b.EN = float(en) * bsolver.TOWNSEND
        b.init()
        f1 = b.converge(f0, maxn=400, rtol=1e-6)
        acc = dict.fromkeys(cols, 0.0)
        kai = 0.0
        kiz = 0.0
        for tgt, proc in b.iter_inelastic():
            if tgt.name != gas:
                continue
            r = float(b.rate(f1, proc)) * tgt.density
            if proc.kind == "IONIZATION":
                kiz += r
            elif proc.kind == "EXCITATION":
                g = group_of(gas, str(proc.product))
                if g is not None:
                    acc[g] += r
                    if float(proc.threshold) >= AI_THRESHOLD_EV[gas]:
                        kai += r
        out["k_ai"].append(max(kai, 0.0))
        out["E_N_Td"].append(float(en))
        out["mean_energy_eV"].append(float(b.mean_energy(f1)))
        out["kiz_m3_s"].append(max(kiz, 0.0))
        for c in cols:
            out[c].append(max(acc[c], 0.0))
        f0 = f1
    return {k: np.asarray(v) for k, v in out.items()}


def path(gas: str) -> Path:
    return SWARM_DIR / f"{gas}_excitation_groups.csv"


def write(gas: str, d: dict) -> Path:
    p = path(gas)
    keys = list(d)
    with p.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# {gas} excitation COUNT rate coefficients per cascade group (m^3/s), two-term BOLOS, DC E/N; "
                 "built by rf_plasma.swarm_excited (same cross sections as swarm_gases)\n")
        fh.write(",".join(keys) + "\n")
        for i in range(len(d[keys[0]])):
            fh.write(",".join(f"{d[k][i]:.6e}" for k in keys) + "\n")
    return p


def load(gas: str) -> dict:
    return read_table(path(gas))


def build(gases=("He", "Ar", "N2")):
    return {g: write(g, table(g)) for g in gases}


if __name__ == "__main__":
    for g, p in build().items():
        print(g, p)
