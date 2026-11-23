"""He / Ar / N2 comparison of the laser-seeded CW-RF breakdown model at 50-200 psig, 300 K, 2.95 GHz.

Everything uses the SAME criterion and machinery as breakdown_probability.py (He model):
  * breakdown when N0 exp(int nu_iz dt) >= N_br = 1e12 electrons, nu_loss = 0 (no attachment in any of
    the three cross-section sets; diffusion only spreads the avalanche, it is not a loss of electrons
    from the cm-scale uniform-field region);
  * nu_iz from the effective field of the RF amplitude (stochastic_seed.Rates.nu_iz):
    E_rms = E0/sqrt2, E_eff = E_rms nu_m / sqrt(nu_m^2 + omega^2), nu_m = N k_m(E_rms/N);
  * field history = breakdown_probability._time_axis: laser at 4.9 us, 100 ns flat RF, then cavity
    ring-down (tau = 35.7 ns) for 10 tau;
  * common seed baseline: the He k=6 seed law (assumed, fitted) at 0.75 mJ -> N0 = 376 electrons. The same
    N0 is used for all gases so that differences are those of the avalanche only. Per-gas multiphoton
    (PPT) seeds are reported separately and used as a sensitivity case.
Thresholds are physical peak fields with e_cal = 1. (The He fit e_cal = 3.912 says the He model needed
~3.9x the nominal measured field; gas RATIOS are the robust output, not absolute values.)

Reference swarm data embedded below are transcribed from A. V. Phelps' JILA transport tables
("eletrans.txt", revised 5/7/98, retrieved via web.archive.org of jila.colorado.edu/~avp/collision_data):
  N2: Phelps & Pitchford, Phys. Rev. A 31, 2932 (1985) / JILA Report 26 -- calculated, stated to agree
      with experiment to +-10 % in drift velocity (to 2000 Td) and ionization (90-1000 Td).
  Ar: Phelps (cross sections revised 10/15/97), calculated; alpha/N stated to fit Kruithof & Penning
      data to ~10 % up to 1000 Td.
  Secondary Ar check: Norman et al., Eur. Phys. J. C 82, 52 (2022) (arXiv:2107.07521), Ar DC breakdown
      projected at 52.6 kV for a 1 cm gap at 10 bar (PyBoltz-augmented Townsend, gamma = 0.01).
"""
from __future__ import annotations

import json
import math
import platform
import sys
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import scipy

from .breakdown_probability import Setup, _time_axis, seed_per_length, simulate, stats
from .ppt import photons_required, ppt_rate_s
from .seed_sources import mpi_per_length
from .stochastic_seed import Rates
from .streamer2d import SWARM_FILES
from .swarm import read_table

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "gas_compare"
PSI = 6894.757
ATM = 101325.0
KB = 1.380649e-23
TD = 1e-21
GASES = ("He", "Ar", "N2")
PSIG = (50.0, 100.0, 150.0, 200.0)
FIELDS_KV_CM = (20.0, 50.0, 100.0, 200.0)
IP_EV = {"He": 24.587, "Ar": 15.760, "N2": 15.581}     # NIST ionization energies
E_CAL_HE = 3.912                                        # breakdown_probability He fit (illustrative)

PHELPS = {
    "N2": dict(
        source="Phelps & Pitchford PRA 31, 2932 (1985); JILA eletrans.txt (rev. 5/7/98), calculated, "
               "+-10% vs experiment (w to 2000 Td, alpha 90-1000 Td)",
        w_m_s={50: 62728, 70: 81653, 100: 107300, 150: 148390, 200: 187380, 300: 260960, 500: 388000,
               1000: 643000},
        mean_energy_eV={100: 2.24, 150: 3.81, 200: 5.03, 300: 6.73, 500: 9.41},
        alphaN_m2={50: 4.29e-28, 70: 9.53e-26, 100: 3.69e-24, 150: 6.93e-23, 200: 2.86e-22, 300: 1.23e-21,
                   500: 4.70e-21, 1000: 1.58e-20},
        kiz_m3_s={100: 3.96e-19, 150: 1.03e-17, 200: 5.36e-17, 300: 3.22e-16, 500: 1.82e-15, 1000: 1.02e-14}),
    "Ar": dict(
        source="Phelps JILA eletrans.txt (Ar cross sections rev. 10/15/97), calculated; alpha/N fits "
               "Kruithof & Penning to ~10% up to 1000 Td",
        w_m_s={10: 10540, 20: 20050, 30: 28300, 50: 43100, 70: 56600, 100: 75800, 200: 134800, 300: 190000,
               500: 294000, 1000: 529000},
        mean_energy_eV={10: 5.37, 20: 5.64, 50: 6.27, 100: 6.84, 200: 7.68},
        alphaN_m2={10: 1.84e-27, 20: 1.58e-24, 30: 1.79e-23, 50: 1.54e-22, 70: 4.46e-22, 100: 9.72e-22,
                   200: 3.30e-21, 300: 5.75e-21, 500: 1.05e-20, 1000: 2.08e-20}),
}
NORMAN_AR = dict(p_pa=10e5, t_k=293.15, gap_m=0.01, gamma=0.01, v_kv_pyboltz=52.6,
                 source="Norman et al., EPJC 82, 52 (2022), arXiv:2107.07521, Table 2 (10 bar, 1 cm)")


def p_abs(psig):
    return psig * PSI + ATM


def rates(gas: str) -> Rates:
    return Rates(gas, path=SWARM_FILES[gas])


def setup(psig: float, **kw) -> Setup:
    return Setup(p_pa=p_abs(psig), **kw)


def seed_n0(st: Setup, energy_j=0.75e-3) -> float:
    """He k=6 seed law (assumed, fitted) electrons at the focus, as breakdown_probability.main."""
    z = np.linspace(-0.03, 0.03, 4801)
    return float(seed_per_length(st, energy_j, z).sum() * (z[1] - z[0]))


# ----------------------------------------------------------------------------- validation
def _interp_loglog(en_tab, y_tab, en):
    m = np.asarray(y_tab) > 0
    return float(np.exp(np.interp(np.log(en), np.log(en_tab[m]), np.log(np.asarray(y_tab)[m]))))


def validation() -> dict:
    out = {}
    for gas, ref in PHELPS.items():
        d = read_table(SWARM_FILES[gas])
        en = np.asarray(d["E_N_Td"])
        cols = dict(w_m_s=np.asarray(d["w_m_s"]), mean_energy_eV=np.asarray(d["mean_energy_eV"]),
                    alphaN_m2=np.asarray(d["alphaN_m2"]), kiz_m3_s=np.asarray(d["kiz_m3_s"]))
        g = dict(source=ref["source"])
        for q in ("w_m_s", "mean_energy_eV", "alphaN_m2", "kiz_m3_s"):
            if q not in ref:
                continue
            rows = []
            for e, v in ref[q].items():
                mod = _interp_loglog(en, cols[q], e)
                rows.append(dict(E_N_Td=e, ref=v, model=mod, dev=mod / v - 1))
            g[q] = rows
        out[gas] = g
    # Ar DC breakdown, Townsend criterion alpha d = ln(1 + 1/gamma); alpha = (kiz/w) N (temporal approx.)
    d = read_table(SWARM_FILES["Ar"])
    en, an = np.asarray(d["E_N_Td"]), np.asarray(d["alphaN_m2"])
    nr = NORMAN_AR
    n = nr["p_pa"] / (KB * nr["t_k"])
    need = math.log(1 + 1 / nr["gamma"]) / (n * nr["gap_m"])
    ok = an > 0
    en_b = float(np.exp(np.interp(math.log(need), np.log(an[ok]), np.log(en[ok]))))
    v = en_b * TD * n * nr["gap_m"]
    out["Ar_DC_breakdown_10bar_1cm"] = dict(nr, E_N_Td=en_b, v_kv_model=v / 1e3,
                                            dev=v / 1e3 / nr["v_kv_pyboltz"] - 1)
    return out


# ----------------------------------------------------------------------------- avalanche / threshold
def gain(r: Rates, st: Setup, e0_v_m: float) -> float:
    """int nu_iz dt from the laser to the end of the ring-down (envelope of breakdown_probability)."""
    t, env = _time_axis(st)
    nu = np.maximum(r.nu_iz(e0_v_m * st.e_cal * env, st.n_gas, st.f_hz) - st.nu_loss_s, 0.0)
    return float(np.trapezoid(nu, t))


def _bisect_log(f, target, lo=1e4, hi=1e9, it=80):
    if f(hi) < target:
        return math.inf
    for _ in range(it):
        mid = math.sqrt(lo * hi)
        if f(mid) >= target:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


def threshold_envelope(r, st: Setup, g_req: float) -> float:
    return _bisect_log(lambda e: gain(r, st, e), g_req)


def threshold_flat(r, st: Setup, g_req: float, window_s: float) -> float:
    return _bisect_log(lambda e: float(r.nu_iz(e, st.n_gas, st.f_hz)) * window_s, g_req)


def eff_en_td(r: Rates, e0, n, f):
    e_rms = e0 / math.sqrt(2)
    nu_m = n * r.km_at(e_rms / n / TD)
    return float(e_rms * nu_m / math.sqrt(nu_m ** 2 + (2 * math.pi * f) ** 2) / n / TD), float(nu_m)


def slope(r, st, e):
    """s = d ln nu_iz / d ln E0 at e."""
    h = 0.01
    a, b = r.nu_iz(e * (1 + h), st.n_gas, st.f_hz), r.nu_iz(e / (1 + h), st.n_gas, st.f_hz)
    return float(math.log(a / b) / (2 * math.log(1 + h)))


def avalanche(n_shots=300, seed=1, overvolt=(1.2, 1.5)) -> dict:
    base = setup(200.0)
    n0 = seed_n0(base)
    g_req = math.log(base.n_br / n0)
    res = dict(n0_common=n0, g_req=g_req, gases={})
    for gas in GASES:
        r = rates(gas)
        rows = []
        for psig in PSIG:
            st = setup(psig, energy_j=0.75e-3)
            n = st.n_gas
            e_env = threshold_envelope(r, st, g_req)
            row = dict(psig=psig, p_kpa=st.p_pa / 1e3, n_m3=n,
                       thr_env_kV_cm=e_env / 1e5,
                       thr_flat50_kV_cm=threshold_flat(r, st, g_req, 50e-9) / 1e5,
                       thr_flat100_kV_cm=threshold_flat(r, st, g_req, 100e-9) / 1e5)
            en_eff, nu_m = eff_en_td(r, e_env, n, st.f_hz)
            row.update(thr_env_Eeff_N_Td=en_eff, thr_env_E0_N_Td=e_env / n / TD, nu_m_at_thr_s=nu_m,
                       nu_m_over_omega=nu_m / (2 * math.pi * st.f_hz), slope_at_thr=slope(r, st, e_env))
            grow = []
            for ekv in FIELDS_KV_CM:
                e0 = ekv * 1e5
                nu = float(r.nu_iz(e0, n, st.f_hz))
                en_eff, nu_m = eff_en_td(r, e0, n, st.f_hz)
                d, _ = simulate(replace(st, e_peak_v_m=e0), n_shots, seed, r)
                neg = nu * 1e-6 < 1e-3          # < 1e-3 e-folds in 1 us: below table noise, "no growth"
                grow.append(dict(E0_kV_cm=ekv, E0_N_Td=e0 / n / TD, Eeff_N_Td=en_eff, nu_m_s=nu_m,
                                 nu_iz_s=None if neg else nu, tau_ns=None if neg else 1e9 / nu,
                                 t_form_flat_ns=None if neg else 1e9 * g_req / nu,
                                 mc=_mc(d)))
            row["growth"] = grow
            ov = []
            for k in overvolt:
                if math.isfinite(e_env):
                    d, _ = simulate(replace(st, e_peak_v_m=k * e_env), n_shots, seed, r)
                    ov.append(dict(overvoltage=k, E0_kV_cm=k * e_env / 1e5, mc=_mc(d)))
            row["overvoltage"] = ov
            rows.append(row)
        res["gases"][gas] = rows
    res["ordering"] = ordering(res)
    return res


def _mc(d):
    s = stats(d)
    ok = np.isfinite(d)
    s["delay_std_ns"] = float(np.std(d[ok]) * 1e9) if ok.sum() > 1 else None
    s["jitter_p90_p10_ns"] = (s["delay_p90_ns"] - s["delay_p10_ns"]) if ok.any() else None
    return s


def ordering(res) -> dict:
    g = res["gases"]
    out = {}
    for i, psig in enumerate(PSIG):
        he, ar, n2 = (g[x][i] for x in GASES)
        o = dict(psig=psig,
                 thr_ratio_N2_He=n2["thr_env_kV_cm"] / he["thr_env_kV_cm"],
                 thr_ratio_Ar_He=ar["thr_env_kV_cm"] / he["thr_env_kV_cm"])
        for k, lab in ((0, "1.2x"), (1, "1.5x")):
            def m(x, key):
                v = x["overvoltage"][k]["mc"][key]
                return v if v is not None else math.nan
            for gas, x in (("Ar", ar), ("N2", n2)):
                o[f"delay_ratio_{gas}_He_at_{lab}_thr"] = m(x, "delay_median_ns") / m(he, "delay_median_ns")
                o[f"jitter_ratio_{gas}_He_at_{lab}_thr"] = m(x, "jitter_p90_p10_ns") / m(he, "jitter_p90_p10_ns")
        out[str(int(psig))] = o
    return out


# ----------------------------------------------------------------------------- laser seed (PPT)
def laser_seed(energies=(0.75e-3, 1e-3)) -> dict:
    """Intrinsic-gas multiphoton seed from PPT (Perelomov-Popov-Terent'ev, as seed_sources.mpi_per_length,
    m = 0, molecules as effective atoms). Order-of-magnitude only: no measured generalized K-photon
    cross section was available for comparison; PPT in the multiphoton regime (Keldysh gamma >> 1) is
    known to be uncertain by orders of magnitude for molecules and for near-resonant intermediate states."""
    st0 = Setup()
    ph_ev = 1239.841984e-9 / st0.wavelength_m
    out = dict(photon_eV=ph_ev, wavelength_nm=st0.wavelength_m * 1e9, w0_um=st0.w0_m * 1e6,
               fwhm_ns=st0.fwhm_s * 1e9, gases={})
    z = np.linspace(-0.03, 0.03, 1201)
    dz = z[1] - z[0]
    variants = {"He": [("Z=1,l=0", 1.0, 0)], "Ar": [("Z=1,l=0", 1.0, 0), ("Z=1,l=1 (3p)", 1.0, 1)],
                "N2": [("Z=1,l=0", 1.0, 0), ("Zeff=0.9,l=0 (Talebpour, unverified)", 0.9, 0)]}
    for gas in GASES:
        rows = []
        for psig in (50.0, 200.0):
            st = setup(psig)
            for en in energies:
                i_pk = st.i_peak(en)
                k = photons_required(i_pk, st.wavelength_m, IP_EV[gas])
                for lab, zeff, l in variants[gas]:
                    if (zeff, l) == (1.0, 0):
                        n_tot = float(mpi_per_length(st, IP_EV[gas], st.n_gas, z, en).sum() * dz)
                    else:
                        n_tot = float(_mpi_total(st, IP_EV[gas], zeff, l, en, z))
                    tau = st.fwhm_s / (2 * math.sqrt(math.log(2)))
                    w_pk = ppt_rate_s(i_pk, st.wavelength_m, IP_EV[gas], zeff, l)
                    p_pk = -math.expm1(-w_pk * tau * math.sqrt(math.pi))
                    rows.append(dict(psig=psig, energy_mJ=en * 1e3, I_peak_W_cm2=i_pk * 1e-4, photons=k,
                                     variant=lab, rate_peak_s=w_pk, n_e_focus_m3=st.n_gas * p_pk,
                                     N0_total=n_tot, he_law_N0=seed_n0(st, en)))
        out["gases"][gas] = rows
    return out


def _mpi_total(st, ip, zeff, l, en, z, nr=48, nt=81):
    from .ppt import ppt_table
    ln_i, ln_w = ppt_table(st.wavelength_m, ip, zeff, l, i_min_w_m2=1e10, i_max_w_m2=1e18, points=161)
    wz = st.w0_m * np.sqrt(1 + ((z - st.focus_m) / st.z_r) ** 2)
    i0 = st.i_peak(en) * (st.w0_m / wz) ** 2
    tau = st.fwhm_s / (2 * math.sqrt(math.log(2)))
    t = np.linspace(-3, 3, nt) * tau
    u = np.linspace(0, 2.0, nr)
    inten = i0[:, None, None] * np.exp(-2 * u[None, :, None] ** 2) * np.exp(-(t[None, None, :] / tau) ** 2)
    lw = np.interp(np.log(np.maximum(inten, 1e-30)), ln_i, ln_w, left=-700.0)
    w = np.exp(np.where(inten < math.exp(ln_i[0]), -700.0, lw))
    prob = -np.expm1(-np.trapezoid(w, t, axis=2))
    per_len = st.n_gas * np.trapezoid(2 * math.pi * u[None, :] * prob, u, axis=1) * wz ** 2
    return per_len.sum() * (z[1] - z[0])


def seed_sensitivity(seed_res, av) -> dict:
    """Envelope threshold with the per-gas PPT seed (0.75 mJ, Z=1, l=0) instead of the common He-law N0."""
    out = {}
    for gas in GASES:
        r = rates(gas)
        rows = []
        for psig in (50.0, 200.0):
            st = setup(psig)
            n0 = [x for x in seed_res["gases"][gas] if x["psig"] == psig and x["energy_mJ"] == 0.75
                  and x["variant"] == "Z=1,l=0"][0]["N0_total"]
            if n0 >= st.n_br:
                rows.append(dict(psig=psig, N0=n0, thr_env_kV_cm=0.0, note="seed alone exceeds N_br"))
                continue
            if n0 < 1:
                rows.append(dict(psig=psig, N0=n0, thr_env_kV_cm=None, note="PPT seed < 1 electron"))
                continue
            g_req = math.log(st.n_br / n0)
            rows.append(dict(psig=psig, N0=n0, g_req=g_req, thr_env_kV_cm=threshold_envelope(r, st, g_req) / 1e5))
        out[gas] = rows
    return out


class ScaledRates:
    """nu_iz scaled by a constant factor (swarm-data uncertainty sensitivity)."""

    def __init__(self, r: Rates, factor: float):
        self.r, self.f = r, factor

    def nu_iz(self, e0, n, f):
        return self.f * self.r.nu_iz(e0, n, f)


def rate_sensitivity(av, factors=(0.5, 2.0)) -> dict:
    """Envelope threshold change when k_iz is scaled (covers the Phelps deviations found in validation)."""
    out = {}
    for gas in GASES:
        r = rates(gas)
        rows = []
        for psig in (50.0, 200.0):
            st = setup(psig)
            base = threshold_envelope(r, st, av["g_req"])
            rows.append(dict(psig=psig, base_kV_cm=base / 1e5, **{
                f"kiz_x{k}_rel": threshold_envelope(ScaledRates(r, k), st, av["g_req"]) / base - 1
                for k in factors}))
        out[gas] = rows
    return out


def rate_factor_to_match(gas="Ar", ref="He", psig_list=PSIG, g_req=None) -> list:
    """Constant multiplier on nu_iz(gas) that brings its envelope threshold down to that of `ref`."""
    if g_req is None:
        g_req = math.log(Setup().n_br / seed_n0(setup(200.0)))
    out = []
    for psig in psig_list:
        st = setup(psig)
        target = threshold_envelope(rates(ref), st, g_req)
        r = rates(gas)
        lo, hi = 1.0, 1e12
        for _ in range(80):
            mid = math.sqrt(lo * hi)
            if threshold_envelope(ScaledRates(r, mid), st, g_req) <= target:
                hi = mid
            else:
                lo = mid
        out.append(dict(psig=psig, factor=math.sqrt(lo * hi), ref_thr_kV_cm=target / 1e5))
    return out


# ----------------------------------------------------------------------------- driver
MISSING_PHYSICS = dict(
    included=["N2 rotational excitation J->J+2 (Boltzmann J populations at Tg = 300 K)",
              "N2 vibrational hot bands v>=1 at Tv = Tg (negligible at 300 K)",
              "N2 ionization to N2+(X) and N2+(B) channels; Ar single ionization; He ionization",
              "effective-field RF heating, nu_m(E/N) from each gas' own two-term momentum-transfer rate"],
    missing=["three-body processes (e.g. Ar+ + 2Ar -> Ar2+, N2+ + 2N2 -> N4+, He+ + 2He): change the ion and "
             "the recombination in the afterglow / channel, not the swarm ionization rate",
             "molecular-ion recombination (Ar2+, N4+, He2+) -- nu_loss = 0 in the breakdown criterion",
             "Ar metastable step-wise ionization and Ar*/N2(A) pooling / associative ionization (lower Ar "
             "threshold further, raise N2 rate late in the avalanche)",
             "N2 vibrational heating (Tv > Tg) and all superelastic collisions (BOLOS 0.2)",
             "electron-electron collisions; spatial (vs temporal) growth; non-local / non-stationary EEDF "
             "during the ns laser pulse",
             "attachment: none in He/Ar/N2 sets (correct for pure gases; O2/H2O impurities would add it)",
             "laser cascade (inverse-bremsstrahlung) breakdown in the focus during the 1 ns pulse",
             "He table energy grid truncated at 120 eV (inaccurate above ~300 Td, outside the He thresholds)"])


def main(out=OUT, n_shots=300):
    out.mkdir(parents=True, exist_ok=True)
    val = validation()
    av = avalanche(n_shots=n_shots)
    sd = laser_seed()
    res = dict(
        description="He/Ar/N2 comparison of laser-seeded CW RF breakdown (rf_plasma.gas_compare)",
        time_origin="t = 0 at the laser pulse peak (4.9 us after RF start); RF flat to 100 ns, then "
                    "ring-down tau = 35.7 ns",
        units=dict(field="kV/cm peak (E0 of E0 cos wt)", E_N="Td (1e-21 V m^2)", time="ns", density="m^-3"),
        criterion="N0 exp(int nu_iz dt) >= 1e12, nu_loss = 0, effective field (breakdown_probability)",
        setup=asdict(Setup()), e_cal=1.0, e_cal_he_fit=E_CAL_HE,
        swarm_files={g: str(SWARM_FILES[g].relative_to(ROOT)) for g in GASES},
        validation=val, avalanche=av, rate_sensitivity=rate_sensitivity(av),
        ar_kiz_factor_to_reach_he_threshold=rate_factor_to_match("Ar", "He", g_req=av["g_req"]),
        laser_seed=sd,
        seed_sensitivity=seed_sensitivity(sd, av),
        physics=MISSING_PHYSICS,
        environment=dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                         platform=platform.platform(), n_shots=n_shots, rng_seed=1,
                         command="PYTHONPATH=src python -m rf_plasma.gas_compare"))
    (out / "gas_compare.json").write_text(json.dumps(res, indent=1, default=_json_default))
    figures(res, out)
    return res


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return str(o)


def figures(res, out=OUT):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"He": "tab:blue", "Ar": "tab:green", "N2": "tab:red"}
    av = res["avalanche"]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    a = ax[0, 0]
    for g in GASES:
        rows = av["gases"][g]
        p = [r["psig"] for r in rows]
        a.plot(p, [r["thr_env_kV_cm"] for r in rows], "o-", c=col[g], label=f"{g} (100 ns + ring-down)")
        a.fill_between(p, [r["thr_flat100_kV_cm"] for r in rows], [r["thr_flat50_kV_cm"] for r in rows],
                       color=col[g], alpha=0.2)
    a.set(xlabel="pressure (psig)", ylabel="threshold peak field (kV/cm)", yscale="log",
          title=f"CW-RF threshold, N0 = {av['n0_common']:.0f}, ln(Nbr/N0) = {av['g_req']:.1f}\n"
                "band: flat 100-50 ns window")
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")
    a = ax[0, 1]
    for g in GASES:
        for r, ls in zip(av["gases"][g], (":", "-.", "--", "-")):
            e = [x["E0_kV_cm"] for x in r["growth"]]
            tau = [np.nan if x["tau_ns"] is None else x["tau_ns"] for x in r["growth"]]
            a.plot(e, tau, ls, c=col[g], marker="o", ms=3,
                   label=f"{g}" if r["psig"] == 200 else None)
    a.set(xlabel="peak field E0 (kV/cm)", ylabel="avalanche e-folding time 1/nu_iz (ns)", xscale="log",
          yscale="log", ylim=(1e-3, 1e5), title="growth time; line style 50 (:) ... 200 psig (-)")
    a.axhspan(100 / av["g_req"], 1e5, color="0.9")
    a.text(21, 100 / av["g_req"] * 1.5, "no breakdown in 100 ns above here", fontsize=7)
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")
    for j, key, lab in ((0, "delay_median_ns", "median delay after laser (ns)"),
                        (1, "jitter_p90_p10_ns", "jitter p90-p10 (ns)")):
        a = ax[1, j]
        for g in GASES:
            for k, mk in ((0, "o-"), (1, "s--")):
                rows = av["gases"][g]
                v = [r["overvoltage"][k]["mc"][key] if r["overvoltage"] else None for r in rows]
                a.plot([r["psig"] for r in rows], [np.nan if x is None else x for x in v], mk, c=col[g],
                       label=f"{g} {rows[0]['overvoltage'][k]['overvoltage']}x thr" if rows[0]["overvoltage"] else None)
        a.set(xlabel="pressure (psig)", ylabel=lab, title=f"Monte Carlo (same N0, sig_p 2%, sig_E 3%)")
        a.legend(fontsize=7)
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "gas_compare.png", dpi=140)
    plt.close(fig)

    val = res["validation"]
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    for g, c in (("Ar", col["Ar"]), ("N2", col["N2"])):
        d = read_table(SWARM_FILES[g])
        en = np.asarray(d["E_N_Td"])
        for a, q, lab in zip(ax, ("w_m_s", "alphaN_m2", "mean_energy_eV"),
                             ("drift velocity (m/s)", "alpha/N (m^2) [model: kiz/w]", "mean energy (eV)")):
            a.plot(en, np.asarray(d[q]), c=c, label=f"{g} BOLOS (this work)")
            rr = val[g].get(q, [])
            a.plot([x["E_N_Td"] for x in rr], [x["ref"] for x in rr], "o", mfc="none", c=c, label=f"{g} Phelps")
            a.set(xscale="log", yscale="log", xlabel="E/N (Td)", ylabel=lab, xlim=(5, 2000))
            a.grid(alpha=0.3, which="both")
    ax[1].set_ylim(1e-28, 1e-19)
    ax[2].set_ylim(0.3, 30)
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "swarm_validation.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    r = main()
    for g in GASES:
        for x in r["avalanche"]["gases"][g]:
            print(g, x["psig"], round(x["thr_env_kV_cm"], 2), [y["tau_ns"] for y in x["growth"]])
    print(json.dumps(r["avalanche"]["ordering"], indent=1))
