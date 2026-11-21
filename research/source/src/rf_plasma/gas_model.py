"""Gas model loader: LXCat-derived swarm tables + cited heavy-particle chemistry -> kernel inputs.

Electron rates (km, kiz, kexc, kel) are tabulated on a log-uniform mean-energy grid:
  bulk      two-term Boltzmann table (rows with energy-balance residual > 2% or non-monotonic
            mean energy are dropped); below its lowest valid energy the Maxwellian table is used.
  maxwell   Maxwellian table (EEDF thermalized by e-e collisions at high ionization degree).
`GasModel.tabs` stacks both (8 rows); the kernels blend them linearly, k = (1-w) kB + w kM, with the
local weight w = nu_ee / (nu_ee + nu_eps). No transition ionization degree is assumed: nu_ee is the
NRL Plasma Formulary electron collision rate and nu_eps = N (kexc/eps + kel + kiz Ip/eps) is the
energy-relaxation frequency from the same LXCat-derived Boltzmann table, so the Maxwellization point
follows from the textbook criterion "e-e energy exchange faster than energy loss to neutrals"
(Raizer 1991; Lieberman & Lichtenberg 2005). The ionization degree ne/N is the control parameter for
Coulomb collisions in two-term Boltzmann solvers (BOLSIG+ manual, Hagelaar & Pitchford 2005, PSST 14
722); a BOLSIG+ table computed at the local ionization degree is the planned swap-in. ee_scale
multiplies nu_ee (default 1) and is varied 0.1-10 only to bound the O(1) uncertainty of the criterion.
A user-supplied BOLSIG+ export can replace either table if it uses the CSV columns in swarm.py.

Heavy-particle chemistry comes from cases/gas_kinetics_model.json. Entries with verified=false
are rejected unless allow_unverified=True; all used entries are returned in `provenance`.
"""

from dataclasses import dataclass, field
import json
import math
from pathlib import Path

import numpy as np
from scipy.constants import atomic_mass, c as C0, elementary_charge as QE, h as PLANCK

from . import plasma3d as k3
from .ppt import photons_required, ppt_table
from .swarm import SWARM_DIR, read_table

MODEL_FILE = Path(__file__).resolve().parents[2] / "cases" / "gas_kinetics_model.json"
TABLE_NAMES = ("km_m3_s", "kiz_m3_s", "kexc_eV_m3_s", "kel_m3_s")


@dataclass
class GasModel:
    species: str
    prm: np.ndarray
    tabs_bulk: np.ndarray      # (4, n) on the log-eps grid
    tabs_filament: np.ndarray
    ln_e0: float
    h: float
    ln_i0: float
    hi: float
    ppt_tab: np.ndarray
    en_td_of_eps: np.ndarray   # Boltzmann E/N (Td) on the grid; nan outside the solved range
    provenance: list = field(default_factory=list)
    tabs_cr: np.ndarray | None = None   # (2, n) He* superelastic / stepwise rows (CR models only)
    optics: dict = field(default_factory=dict)
    tabs_ee: np.ndarray | None = None   # (4 * nlev, n) e-e Boltzmann blocks per ionization degree

    @property
    def tabs(self) -> np.ndarray:
        """(8, n): Boltzmann rows 0-3, Maxwellian rows 4-7 (blended in the kernels); rows 8-9 CR."""
        rows = [self.tabs_bulk, self.tabs_filament] + ([self.tabs_cr] if self.tabs_cr is not None else [])
        rows += [self.tabs_ee] if self.tabs_ee is not None else []
        return np.ascontiguousarray(np.vstack(rows))

    @property
    def eps_grid(self) -> np.ndarray:
        return np.exp(self.ln_e0 + self.h * np.arange(self.tabs_bulk.shape[1]))


def _clean(d: dict, boltzmann: bool) -> dict:
    eps = d["mean_energy_eV"]
    keep = np.isfinite(eps) & np.isfinite(d["kel_m3_s"])
    if boltzmann:
        keep &= np.abs(d["balance_residual"]) < 0.02
    idx = np.where(keep)[0]
    mono = [idx[0]]
    for i in idx[1:]:
        if eps[i] > eps[mono[-1]] * 1.0005:
            mono.append(i)
    return {k: v[mono] for k, v in d.items()}


def _resample(d: dict, eps_grid: np.ndarray) -> np.ndarray:
    x = np.log(d["mean_energy_eV"])
    lg = np.log(eps_grid)
    out = np.empty((4, len(eps_grid)))
    for r, name in enumerate(TABLE_NAMES):
        y = np.maximum(d[name], 0.0)
        if name in ("km_m3_s", "kel_m3_s"):
            out[r] = np.exp(np.interp(lg, x, np.log(np.maximum(y, 1e-30))))
        else:  # threshold rates: log-interpolate with a tiny floor, then zero the floor
            ly = np.interp(lg, x, np.log(np.maximum(y, 1e-300)))
            out[r] = np.where(ly < math.log(1e-290), 0.0, np.exp(ly))
    return out


def _entry(spec, key, allow_unverified, prov, species, required=True):
    if key not in spec:
        if required:
            raise KeyError(f"{species}: missing required entry '{key}'")
        return None
    ent = spec[key]
    ok = bool(ent.get("verified", False))
    prov.append({"species": species, "entry": key, "verified": ok, "source": ent.get("source"),
                 "note": ent.get("note")})
    if not ok and not allow_unverified:
        raise ValueError(f"{species}: '{key}' is not verified ({ent.get('source')}); "
                         "pass allow_unverified=True to run with it")
    return ent


def load_gas_model(species: str, wavelength_m: float, allow_unverified: bool = False,
                   model_file: Path = MODEL_FILE, bulk_table: Path | None = None,
                   filament_table: Path | None = None, eps_min: float = 0.03, eps_max: float = 40.0,
                   points: int = 500, ee_scale: float = 1.0, extended: bool = False, cr: bool = False,
                   transport: bool = False, ee_n_gas: float | None = None) -> GasModel:
    spec = json.loads(Path(model_file).read_text(encoding="utf-8"))["species"][species]
    prov: list = []

    def get(key, required=True):
        return _entry(spec, key, allow_unverified, prov, species, required)

    bulk_path = Path(bulk_table) if bulk_table else SWARM_DIR / f"{species}_boltzmann.csv"
    fil_path = Path(filament_table) if filament_table else SWARM_DIR / f"{species}_maxwell.csv"
    boltz = _clean(read_table(bulk_path), boltzmann=True)
    maxw = _clean(read_table(fil_path), boltzmann=False)
    prov.append({"species": species, "entry": "electron_rates_bulk", "verified": True, "source": str(bulk_path)})
    prov.append({"species": species, "entry": "electron_rates_filament", "verified": True, "source": str(fil_path)})

    ln_e0 = math.log(eps_min)
    hh = (math.log(eps_max) - ln_e0) / (points - 1)
    eps = np.exp(ln_e0 + hh * np.arange(points))
    tb = _resample(boltz, eps)
    tm = _resample(maxw, eps)
    below = eps < boltz["mean_energy_eV"][0]
    tb[:, below] = tm[:, below]
    en = np.interp(np.log(eps), np.log(boltz["mean_energy_eV"]), boltz["E_N_Td"], left=np.nan, right=np.nan)

    prm = np.zeros(k3.N_PARAMS)
    prm[k3.P_MASS] = float(get("mass_amu")["value"]) * atomic_mass
    prm[k3.P_IP] = float(get("ip_eV")["value"])
    prm[k3.P_DOF] = float(get("dof")["value"])
    kc = get("kconv_m6_s")
    prm[k3.P_KCONV0], prm[k3.P_KCONV_EXP] = kc["k0"], kc["tg_exp"]
    b1 = get("beta_atomic_m3_s")
    prm[k3.P_B1_0], prm[k3.P_B1_EXP] = b1["b0"], b1["te_exp"]
    b1n = get("beta_atomic_neutral_m6_s")
    prm[k3.P_B1N_0], prm[k3.P_B1N_EXP] = b1n["b0"], b1n["te_exp"]
    b2 = get("beta_dimer_m3_s")
    prm[k3.P_B2_0], prm[k3.P_B2_EXP] = b2["b0"], b2["te_exp"]
    prm[k3.P_B2_THETA], prm[k3.P_B2_C] = b2.get("theta_K", 0.0), b2.get("c", 0.0)
    b2n = get("beta_dimer_neutral_m6_s")
    prm[k3.P_B2N_0], prm[k3.P_B2N_EXP] = b2n["b0"], b2n["te_exp"]
    k3b = get("k3b_m6_s")
    prm[k3.P_K3B0], prm[k3.P_K3B_EXP] = k3b["k0"], k3b["te_exp"]
    mob = get("mobility_K0_m2_Vs")
    prm[k3.P_MU_XP], prm[k3.P_MU_X2P] = mob["Xp"], mob["X2p"]
    prm[k3.P_ETA] = float(get("eta_fast")["value"])
    d = get("dimer_dissociation", False)
    if d:
        prm[k3.P_KDISS0], prm[k3.P_EDISS] = d["k0"], d["E_eV"]
        prm[k3.P_DISS_A], prm[k3.P_DISS_T1] = d.get("a", 0.0), d.get("T1_K", 0.0)
    d = get("dimer_dissociative_excitation", False)
    if d:
        prm[k3.P_KDEX0], prm[k3.P_EDEX] = d["k0"], d["E_eV"]
    prm[k3.P_EE_SCALE] = ee_scale
    prov.append({"species": species, "entry": "eedf_closure", "verified": True,
                 "source": "w = nu_ee/(nu_ee+nu_eps): nu_ee NRL Plasma Formulary (2.91e-12 ne lnL Te^-1.5, SI); "
                           "nu_eps from the LXCat-derived Boltzmann table; criterion Raizer 1991, "
                           "Lieberman & Lichtenberg 2005; ne/N as Coulomb control parameter: BOLSIG+ manual",
                 "note": f"ee_scale = {ee_scale:g}; interpolation form k = (1-w) kB + w kM is a modeling choice"})
    prm[k3.P_PHOTON_EV] = PLANCK * C0 / wavelength_m / QE
    prm[k3.P_OMEGA_L] = 2 * math.pi * C0 / wavelength_m
    prm[k3.P_NPHOT] = photons_required(1e12, wavelength_m, prm[k3.P_IP])
    p = get("ppt")
    ln_i, ln_rate = ppt_table(wavelength_m, prm[k3.P_IP], float(p["z_eff"]), int(p["l"]))
    gm = GasModel(species, prm, np.ascontiguousarray(tb), np.ascontiguousarray(tm), ln_e0, hh,
                  float(ln_i[0]), float(ln_i[1] - ln_i[0]), np.ascontiguousarray(ln_rate), en, prov)
    if extended:
        _extend(gm, spec, get, eps, cr, transport)
        if ee_n_gas:
            _extend_ee(gm, maxw, eps, ee_n_gas, prov)
    return gm


def _extend_ee(gm: GasModel, maxw: dict, eps: np.ndarray, n_gas: float, prov: list):
    """Replace the w-blend closure by solved e-e Boltzmann tables (swarm_ee), one block per ne/N level."""
    from . import plasma3d_cr as kc
    from .swarm_ee import X_LEVELS, build, nearest_complete
    tm = _resample(maxw, eps)
    blocks = []
    n_req, n_gas = n_gas, nearest_complete(gm.species, n_gas)
    for path in build(gm.species, n_gas):
        d = _clean(read_table(path), boltzmann=True)
        tb = _resample(d, eps)
        below = eps < d["mean_energy_eV"][0]
        tb[:, below] = tm[:, below]
        blocks.append(tb)
    gm.tabs_ee = np.ascontiguousarray(np.vstack(blocks))
    gm.prm[kc.P_EE_ROW0] = 8 + (2 if gm.tabs_cr is not None else 0)
    gm.prm[kc.P_EE_NLEV] = len(X_LEVELS)
    prov.append({"species": gm.species, "entry": "eedf_closure_ee", "verified": True,
                 "source": "two-term Boltzmann + e-e collisions (Hagelaar & Pitchford 2005 PSST 14 722, sec. 2.6) "
                           f"solved with BOLOS; levels ne/N = {list(X_LEVELS)}; table N = {n_gas:.4e} m^-3 "
                           f"(run N = {n_req:.4e}; N enters only via ln Lambda)",
                 "note": "replaces the w = nu_ee/(nu_ee+nu_eps) interpolation; CR superelastic/stepwise rows remain Maxwellian"})


def _extend(gm: GasModel, spec: dict, get, eps: np.ndarray, cr: bool, transport: bool):
    """Extended parameter vector for plasma3d_cr; CR (helium only) and transport mobilities."""
    from . import plasma3d_cr as kc
    prm = kc.extend_params(gm.prm)
    if transport:
        mob = get("transport_mobility_K0_m2_Vs")
        prm[k3.P_MU_XP], prm[k3.P_MU_X2P] = mob["Xp"], mob["X2p"]
    o = get("optics", False)
    if o:
        gm.optics = dict(refractivity_stp=o["refractivity_stp"], n2_m2_W_per_atm=o["n2_m2_W_per_atm"])
    if cr:
        if gm.species != "He":
            raise ValueError("the lumped CR model is implemented for helium only")
        from .helium_cr import E_MANIFOLD, G_MANIFOLD, manifold_tables
        c = get("cr_model")
        k_sup, k_step = manifold_tables(eps)
        gm.tabs_cr = np.ascontiguousarray(np.vstack((k_sup, k_step)))
        ip = prm[k3.P_IP]
        prm[kc.P_CR_ON] = 1.0
        prm[kc.P_ESTAR] = E_MANIFOLD
        prm[kc.P_DION] = c["dion_eV"]["value"]
        prm[kc.P_EX] = ip - c["dion_eV"]["value"] - c["kxi"]["E_eV"]
        prm[kc.P_I2] = c["i2_eV"]["value"]
        prm[kc.P_KC3] = c["kc3_m6_s"]["value"]
        prm[kc.P_KMM_A], prm[kc.P_KMM_M] = c["kmm_atomic_m3_s"]["value"], c["kmm_dimer_m3_s"]["value"]
        prm[kc.P_KMX_A], prm[kc.P_KMX_M] = c["kmx_atomic_m3_s"]["value"], c["kmx_dimer_m3_s"]["value"]
        prm[kc.P_KXX] = c["kxx_m3_s"]["value"]
        prm[kc.P_KXQ] = c["kxq_m3_s"]["value"]
        prm[kc.P_KXI0], prm[kc.P_KXIE] = c["kxi"]["k0"], c["kxi"]["E_eV"]
        prm[kc.P_AX] = c["excimer_rad_s"]["value"]
        prm[kc.P_KAREC0], prm[kc.P_KAREC_EXP] = c["alpha_rec_m6_s"]["k0"], c["alpha_rec_m6_s"]["te_exp"]
        # resonance broadening (NIST AtSpec eq. 41) -> density-independent Lorentz line-centre k0
        lam, f_r, a_r = c["resonance_lambda_m"], c["resonance_f"], c["resonance_A_s"]
        lam_a = lam * 1e10
        dlam_per_n = 8.6e-30 * math.sqrt(1 / 3) * lam_a ** 3 * f_r * 1e-10 / 1e6  # m per m^-3
        dnu_per_n = C0 * dlam_per_n / lam ** 2                                      # Hz per m^-3
        k0 = lam ** 2 / (8 * math.pi) * 3.0 * a_r * 2.0 / (math.pi * dnu_per_n)    # 1/m, N cancels
        g = min(1.0, 1.115 / math.sqrt(math.pi * k0 * c["trap_radius_m"]))
        prm[kc.P_ARES] = a_r * g
        prm[kc.P_FRES] = 3.0 / G_MANIFOLD
        gm.optics["resonance_k0_per_m"] = k0
        gm.optics["resonance_escape_factor"] = g
    gm.prm = prm
