"""Electron swarm/rate tables vs mean energy from LXCat cross sections (BOLOS two-term solver).

Two tables per gas, both built from the same cross-section set:
  boltzmann: steady two-term EEDF for a DC reduced field E/N (valid when nu_m >> omega_RF,
             which holds at ~8 bar: nu_m ~ 1e12-1e13 s^-1 vs 2 pi 3 GHz).
  maxwell:   Maxwellian EEDF at Te = 2 eps / 3 (appropriate where electron-electron collisions
             thermalize the EEDF, e.g. the laser-filament core at ionization degree >~ 1e-4).

Columns (per unit gas density N): mean_energy_eV, E_N_Td (boltzmann only), muN_1_Vms,
km_m3_s = e / (m_e muN) [conductivity-equivalent momentum rate], kiz_m3_s,
kexc_eV_m3_s = sum(threshold * rate) over excitation + vibrational channels,
kel_m3_s = elastic power / (N n_e (eps - 1.5 kTg)) so the kernel loss is N n_e kel (eps - eps_g).
"""

from pathlib import Path
import math

import numpy as np
import scipy.constants as co
from scipy import integrate

if not hasattr(integrate, "simps"):  # BOLOS 0.1 uses the removed scipy alias
    integrate.simps = integrate.simpson

from bolos import grid as bgrid, parser as bparser, solver as bsolver  # noqa: E402

CROSS_SECTION_DIR = Path(__file__).resolve().parents[2] / "cases" / "cross_sections"
SWARM_DIR = Path(__file__).resolve().parents[2] / "cases" / "swarm"
CROSS_SECTION_FILES = {"Ar": "LoKI_Ar_LXCat.txt", "He": "LoKI_He_LXCat.txt", "N2": "LoKI_N2_LXCat.txt"}
SOURCE_NOTE = ("IST-Lisbon database, www.lxcat.net (files distributed with LoKI-B, "
               "github.com/IST-Lisbon/LoKI, retrieved 2022-09-23 by the LoKI authors); "
               "two-term Boltzmann solution by BOLOS (A. Luque, github.com/aluque/bolos)")
COLUMNS = ["mean_energy_eV", "E_N_Td", "muN_1_Vms", "km_m3_s", "kiz_m3_s", "kexc_eV_m3_s", "kel_m3_s",
           "balance_residual"]


def make_solver(gas: str, e_max: float = 120.0, cells: int = 400, tg_k: float = 300.0):
    path = CROSS_SECTION_DIR / CROSS_SECTION_FILES[gas]
    with path.open(encoding="utf-8", errors="replace") as fh:
        procs = bparser.parse(fh)
    b = bsolver.BoltzmannSolver(bgrid.QuadraticGrid(0.0, e_max, cells))
    b.load_collisions(procs)
    b.target[gas].density = 1.0
    b.kT = tg_k * co.k / co.eV
    b.EN = 1.0 * bsolver.TOWNSEND
    b.init()
    return b


def _row(b, f0, en_td):
    eps = float(b.mean_energy(f0))
    mun = float(b.mobility(f0))
    kiz = 0.0
    kizl = 0.0
    kexc = 0.0
    for _, proc in b.iter_inelastic():
        r = float(b.rate(f0, proc))
        if proc.kind == "IONIZATION":
            kiz += r
            kizl += r * float(proc.threshold)
        elif proc.kind == "EXCITATION":
            kexc += r * float(proc.threshold)
    fb = np.interp(b.benergy, b.cenergy, f0)
    dfb = np.gradient(fb, b.benergy)
    gam = math.sqrt(2 * co.e / co.m_e)
    p_el = gam * integrate.simpson(b.benergy**2 * b.sigma_eps * (fb + b.kT * dfb), x=b.benergy)
    eps_g = 1.5 * b.kT
    kel = p_el / (eps - eps_g) if eps > eps_g * 1.02 else np.nan
    p_in = (en_td * bsolver.TOWNSEND) ** 2 * mun if en_td is not None else np.nan
    resid = (p_in - p_el - kexc - kizl - kiz * eps) / p_in if en_td else np.nan
    return [eps, en_td if en_td is not None else np.nan, mun, co.e / (co.m_e * mun), kiz, kexc, kel, resid]


def boltzmann_table(gas: str, en_td=None):
    en_td = np.geomspace(0.03, 3000.0, 90) if en_td is None else en_td
    b = make_solver(gas)
    f0 = b.maxwell(0.5)
    rows = []
    for en in en_td:
        b.EN = en * bsolver.TOWNSEND
        b.init()
        f1 = b.converge(f0, maxn=400, rtol=1e-6)
        rows.append(_row(b, f1, float(en)))
        f0 = f1
    return np.array(rows)


def maxwell_table(gas: str, te_ev=None):
    te_ev = np.geomspace(0.03, 40.0, 120) if te_ev is None else te_ev
    b = make_solver(gas, e_max=600.0, cells=1200)
    rows = []
    for te in te_ev:
        f0 = b.maxwell(float(te))
        f0 = f0 / integrate.simpson(f0 * np.sqrt(b.cenergy), x=b.cenergy)
        rows.append(_row(b, f0, None))
    return np.array(rows)


def write_table(path: Path, data: np.ndarray, gas: str, kind: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# {gas} electron rate table ({kind} EEDF), per unit gas density\n")
        fh.write(f"# source: {SOURCE_NOTE}\n")
        fh.write("# km = e/(m_e muN); kexc = sum threshold*rate (excitation+vibrational); "
                 "elastic loss = N ne kel (eps - 1.5 k Tg), Tg = 300 K\n")
        fh.write(",".join(COLUMNS) + "\n")
        for r in data:
            fh.write(",".join(f"{v:.6e}" for v in r) + "\n")


def read_table(path: Path) -> dict:
    lines = [ln for ln in Path(path).read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    head = lines[0].split(",")
    arr = np.array([[float(x) for x in ln.split(",")] for ln in lines[1:]])
    return {h: arr[:, i] for i, h in enumerate(head)}


def build_all(gases=("Ar", "He", "N2")) -> dict:
    out = {}
    for gas in gases:
        for kind, fn in (("boltzmann", boltzmann_table), ("maxwell", maxwell_table)):
            path = SWARM_DIR / f"{gas}_{kind}.csv"
            write_table(path, fn(gas), gas, kind)
            out[(gas, kind)] = path
    return out


if __name__ == "__main__":
    for key, p in build_all().items():
        print(key, p)
