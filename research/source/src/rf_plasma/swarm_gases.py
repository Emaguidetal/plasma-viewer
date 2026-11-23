"""Gas-selectable electron swarm tables (He, Ar, N2) for streamer2d / rod3d / breakdown models.

Same two-term BOLOS solver, energy grid (quadratic, 0-120 eV, 400 cells), gas temperature (300 K) and
E/N grid (90 log points, 0.03-3000 Td) as the existing He table (cases/swarm/He_boltzmann.csv, built by
swarm.py). Output CSV keeps the He columns (so streamer2d.Swarm and stochastic_seed.Rates read it
unchanged) and appends:
  katt_m3_s   attachment rate coefficient (0: no ATTACHMENT process in the LoKI He/Ar/N2 sets; N2- is a
              short-lived shape resonance, dissociative attachment to N2 needs > 9 eV and is absent)
  w_m_s       drift velocity  = muN * E/N
  alphaN_m2   reduced Townsend ionization coefficient approximated as kiz / w (temporal-growth swarm)
  DN_1_ms     density-normalised (transverse-type) two-term diffusion coefficient N*D

Cross sections (IST-Lisbon database, www.lxcat.net, files distributed with LoKI-B):
  He: LoKI_He_LXCat.txt                       Ar: LoKI_Ar_LXCat.txt
  N2: LoKI_N2_LXCat.txt (EFFECTIVE + v=0->v' vibrational + electronic + 2 ionization channels)
      + LoKI_N2_rot_LXCat.txt (J -> J+2, Gerjuoy-Stein), J populations Boltzmann at Tg with
        rigid-rotor energies B J(J+1) and nuclear-spin weights 6 (even J) / 3 (odd J), as in the
        LoKI-B N2 swarm set-up (Code/Input/Nitrogen/N2_swarm_setup.in)
      + LoKI_N2_vib_LXCat.txt hot bands v -> v' (v >= 1), populations of a harmonic oscillator
        (omega_e = 0.2924 eV) at Tv (default Tv = Tg = 300 K -> f(v=1) = 1.3e-5, negligible).
  The N2 EFFECTIVE (total momentum-transfer) cross section is converted to ELASTIC by subtracting every
  inelastic cross section of N2, rotational and hot-band ones weighted by their populations (the LoKI
  convention). BOLOS 0.2 does NOT implement superelastic collisions: the '<->' de-excitations of the
  rotational and vibrational sets are omitted. This matters below ~1 Td (rotations) and whenever
  Tv >> Tg (vibrational heating); it is negligible at the 30-200 Td relevant for RF breakdown at Tv = Tg.
Growth model: BOLOS temporal growth, equal energy sharing at ionization (LoKI swarm runs use spatial
growth and SDCS sharing; differences are a few % below alpha/N ~ 1e-21 m^2).
"""
from __future__ import annotations

import logging
import math
import re
from pathlib import Path

import numpy as np
import scipy.constants as co
from scipy import integrate

from .swarm import CROSS_SECTION_DIR, SOURCE_NOTE, SWARM_DIR, COLUMNS

if not hasattr(integrate, "simps"):  # BOLOS uses the removed scipy alias
    integrate.simps = integrate.simpson

from bolos import grid as bgrid, parser as bparser, solver as bsolver  # noqa: E402

EN_GRID_TD = np.geomspace(0.03, 3000.0, 90)        # identical to swarm.boltzmann_table / He_boltzmann.csv
EXTRA_COLUMNS = ["katt_m3_s", "w_m_s", "alphaN_m2", "DN_1_ms"]
GAS_FILES = {
    "He": dict(base="LoKI_He_LXCat.txt"),
    "Ar": dict(base="LoKI_Ar_LXCat.txt"),
    "N2": dict(base="LoKI_N2_LXCat.txt", rot="LoKI_N2_rot_LXCat.txt", vib="LoKI_N2_vib_LXCat.txt"),
}
# N2 X state constants (NIST WebBook / Huber-Herzberg): B0 = 1.98957 cm^-1, omega_e = 2358.57 cm^-1
CM1_EV = 1.239841984e-4
B_N2_EV = 1.98957 * CM1_EV
WE_N2_EV = 2358.57 * CM1_EV
KB_EV = co.k / co.e


def swarm_path(gas: str, tag: str = "swarm") -> Path:
    return SWARM_DIR / f"{gas}_{tag}.csv"


# ----------------------------------------------------------------------------- new-format LXCat blocks
_RE_PROC = re.compile(r"PROCESS:\s*E \+ (\S+) <?->\s*E \+ (\S+) \(([vJ])=(\d+)-[vJ]=(\d+)\)")
_RE_PAR = re.compile(r"PARAM\.:\s*E = ([0-9.eE+-]+) eV(?:, g1/g0 = ([0-9.eE+-]+))?")


def parse_state_resolved(path: Path):
    """Blocks of the newer LXCat export (ROTATIONAL / VIBRATIONAL keyword, PROCESS/PARAM lines) ->
    list of dicts(kind, quantum ('J'/'v'), lower, upper, threshold, g_ratio, data[n,2])."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    out, i = [], 0
    while i < len(lines):
        key = lines[i].strip()
        if key in ("ROTATIONAL", "VIBRATIONAL"):
            proc = par = None
            j = i + 1
            while not lines[j].strip().startswith("-----"):
                proc = proc or _RE_PROC.search(lines[j])
                par = par or _RE_PAR.search(lines[j])
                j += 1
            k = j + 1
            data = []
            while not lines[k].strip().startswith("-----"):
                data.append([float(x) for x in lines[k].split()])
                k += 1
            if proc is None or par is None:
                raise ValueError(f"unparsed block at line {i + 1} of {path}")
            out.append(dict(kind=key, quantum=proc.group(3), lower=int(proc.group(4)), upper=int(proc.group(5)),
                            threshold=float(par.group(1)),
                            g_ratio=float(par.group(2)) if par.group(2) else 1.0, data=np.array(data)))
            i = k + 1
        else:
            i += 1
    return out


def n2_rot_populations(tg_k: float, jmax: int = 60) -> np.ndarray:
    j = np.arange(jmax + 1)
    g = (2 * j + 1) * np.where(j % 2 == 0, 6.0, 3.0)
    w = g * np.exp(-B_N2_EV * j * (j + 1) / (KB_EV * tg_k))
    return w / w.sum()


def n2_vib_populations(tv_k: float, vmax: int = 60) -> np.ndarray:
    v = np.arange(vmax + 1)
    x = math.exp(-WE_N2_EV / (KB_EV * tv_k))
    return (1 - x) * x ** v


# ----------------------------------------------------------------------------- solver
def make_solver(gas: str, e_max: float = 120.0, cells: int = 400, tg_k: float = 300.0, tv_k: float | None = None,
                rot: bool = True, vib: bool = True):
    spec = GAS_FILES[gas]
    with (CROSS_SECTION_DIR / spec["base"]).open(encoding="utf-8", errors="replace") as fh:
        procs = bparser.parse(fh)
    extra = []          # (target name, fraction, process dict)
    if gas == "N2":
        tv_k = tg_k if tv_k is None else tv_k
        blocks = []
        if rot and "rot" in spec:
            blocks += parse_state_resolved(CROSS_SECTION_DIR / spec["rot"])
        if vib and "vib" in spec:
            blocks += parse_state_resolved(CROSS_SECTION_DIR / spec["vib"])
        prot, pvib = n2_rot_populations(tg_k), n2_vib_populations(tv_k)
        for b in blocks:
            frac = (prot if b["quantum"] == "J" else pvib)[b["lower"]]
            name = f"N2({b['quantum']}={b['lower']})"
            extra.append((name, float(frac), dict(kind="EXCITATION", target=name,
                                                  product=f"N2({b['quantum']}={b['upper']})",
                                                  threshold=b["threshold"], data=b["data"].tolist(),
                                                  comment="IST-Lisbon state-resolved set")))
        # EFFECTIVE -> subtract population-weighted state-resolved inelastic cross sections (LoKI convention)
        eff = [p for p in procs if p["kind"] == "EFFECTIVE" and p["target"] == "N2"][0]
        d = np.array(eff["data"], float)
        e_all = np.unique(np.concatenate([d[:, 0]] + [np.array(p["data"])[:, 0] for _, _, p in extra]))
        sig = np.interp(e_all, d[:, 0], d[:, 1])
        for _, frac, p in extra:
            dd = np.array(p["data"])
            sig -= frac * np.interp(e_all, dd[:, 0], dd[:, 1], left=0.0, right=dd[-1, 1])
        eff["data"] = np.c_[e_all, np.maximum(sig, 0.0)].tolist()
    b = bsolver.BoltzmannSolver(bgrid.QuadraticGrid(0.0, e_max, cells))
    level = logging.getLogger().level
    logging.getLogger().setLevel(logging.ERROR)  # 'target has no ELASTIC' warnings for state targets
    try:
        b.load_collisions(procs + [p for _, _, p in extra])
    finally:
        logging.getLogger().setLevel(level)
    b.target[gas].density = 1.0
    for name, frac, _ in extra:
        b.target[name].density = frac
    b.kT = tg_k * co.k / co.eV
    b.EN = 1.0 * bsolver.TOWNSEND
    b.init()
    return b


def row(b, f0, en_td):
    """Same definitions as swarm._row, but every rate is weighted by its target density (state-resolved
    targets have fractional populations), plus the extra columns."""
    eps = float(b.mean_energy(f0))
    mun = float(b.mobility(f0))
    kiz = kizl = kexc = katt = 0.0
    for tgt, proc in b.iter_inelastic():
        r = float(b.rate(f0, proc)) * tgt.density
        if proc.kind == "IONIZATION":
            kiz += r
            kizl += r * float(proc.threshold)
        elif proc.kind == "EXCITATION":
            kexc += r * float(proc.threshold)
        elif proc.kind == "ATTACHMENT":
            katt += r
    fb = np.interp(b.benergy, b.cenergy, f0)
    dfb = np.gradient(fb, b.benergy)
    gam = math.sqrt(2 * co.e / co.m_e)
    p_el = gam * integrate.simpson(b.benergy ** 2 * b.sigma_eps * (fb + b.kT * dfb), x=b.benergy)
    eps_g = 1.5 * b.kT
    kel = p_el / (eps - eps_g) if eps > eps_g * 1.02 else np.nan
    p_in = (en_td * bsolver.TOWNSEND) ** 2 * mun
    resid = (p_in - p_el - kexc - kizl - kiz * eps) / p_in
    w = mun * en_td * bsolver.TOWNSEND
    dn = float(b.diffusion(f0))
    return [eps, en_td, mun, co.e / (co.m_e * mun), kiz, kexc, kel, resid, katt, w, kiz / w, dn]


def table(gas: str, en_td=EN_GRID_TD, **kw) -> np.ndarray:
    b = make_solver(gas, **kw)
    f0 = b.maxwell(0.5)
    rows = []
    for en in en_td:
        b.EN = float(en) * bsolver.TOWNSEND
        b.init()
        f1 = b.converge(f0, maxn=400, rtol=1e-6)
        rows.append(row(b, f1, float(en)))
        f0 = f1
    return np.array(rows)


def write(path: Path, data: np.ndarray, gas: str, note: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# {gas} electron swarm table (two-term Boltzmann EEDF, DC E/N), per unit gas density; {note}\n")
        fh.write(f"# source: {SOURCE_NOTE}; built by rf_plasma.swarm_gases\n")
        fh.write("# km = e/(m_e muN); kexc = sum threshold*rate (excitation incl. vib/rot, population weighted); "
                 "elastic loss = N ne kel (eps - 1.5 k Tg), Tg = 300 K; w = muN E/N; alphaN = kiz/w; "
                 "katt = 0 (no attachment in set); no superelastic collisions (BOLOS)\n")
        fh.write(",".join(COLUMNS + EXTRA_COLUMNS) + "\n")
        for r in data:
            fh.write(",".join(f"{v:.6e}" for v in r) + "\n")


NOTES = {"He": "LoKI He set (consistency copy; streamer2d/rod3d keep He_boltzmann.csv)",
         "Ar": "LoKI Ar set",
         "N2": "LoKI N2 set + rotational (Boltzmann J, Tg) + hot-band vibrational (Tv = Tg) cross sections"}


def build(gases=("Ar", "N2", "He")) -> dict:
    out = {}
    for g in gases:
        p = swarm_path(g)
        write(p, table(g), g, NOTES[g])
        out[g] = p
    return out


if __name__ == "__main__":
    for g, p in build().items():
        print(g, p)
