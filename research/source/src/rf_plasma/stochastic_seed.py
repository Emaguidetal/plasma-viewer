"""Stochastic seed-electron / breakdown-delay surrogate (statistical lag + formative time + jitter).

Breakdown in a pre-filled RF cavity needs (i) a first effective electron and (ii) an avalanche from the
seed population to a detectable/reflecting electron number N_br. Both are stochastic:

* RF-only: background electrons appear as a Poisson process with rate R = S_bg * V_hf (S_bg: natural
  ion-pair production by cosmic rays + terrestrial radioactivity, ~10 cm^-3 s^-1 in STP air
  [Raizer, Gas Discharge Physics, 1991, sec. 7.1], scaled to helium by mass stopping power and W-value:
  rho_He/rho_air * (dE/dx)_He/(dE/dx)_air * W_air/W_He = (4/28.96)(1.94/1.82)(34/41) = 0.12 at equal N;
  PDG minimum-ionizing dE/dx and ICRU-31 W-values). Each electron starts a lineage that survives with
  probability 1 - nu_loss/nu_iz.
* Laser-triggered: the pulse leaves N_L electrons (N_L ~ (E_L/E_L0)^k, k = effective multiphoton/cascade
  order fitted to the 3D filament runs) at t_L.
* Avalanche: linear birth-death process (birth nu_iz, death nu_loss). N(t) exp(-int nu_net dt) -> W with
  W = sum of K i.i.d. Exp(nu_iz/nu_net) variables, K ~ Binomial(n0, 1 - nu_loss/nu_iz) surviving lineages
  (Kendall 1948, Ann. Math. Stat. 19, 1). Breakdown when int_0^t nu_net dt' = ln(N_br / W).

Rates come from the two-term Boltzmann table (BOLOS + LXCat) at the effective field E_eff = E0/sqrt(2)
* nu_m/sqrt(nu_m^2 + omega^2): in helium the energy relaxation frequency (2m/M) nu_m ~ 1e8-1e9 s^-1 is
far below omega, so the EEDF responds to the RMS field (effective-field concept, Raizer sec. 7.2).
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from .swarm import read_table

KB = 1.380649e-23
SWARM = Path(__file__).resolve().parents[2] / "cases" / "swarm"
N_ATM = 101325.0 / (KB * 293.15)
S_AIR_STP = 10e6           # ion pairs m^-3 s^-1 (10 cm^-3 s^-1)
HE_OVER_AIR = (4.0 / 28.96) * (1.94 / 1.82) * (34.0 / 41.0)


class Rates:
    def __init__(self, gas="He", closure="boltzmann"):
        d = read_table(SWARM / f"{gas}_{closure}.csv")
        en = np.asarray(d["E_N_Td"])
        o = np.argsort(en)
        self.ln_en = np.log(en[o])
        self.ln_kiz = np.log(np.maximum(np.asarray(d["kiz_m3_s"])[o], 1e-60))
        self.km = np.asarray(d["km_m3_s"])[o]

    def kiz(self, en_td):
        x = np.log(np.clip(en_td, math.exp(self.ln_en[0]), math.exp(self.ln_en[-1])))
        return np.exp(np.interp(x, self.ln_en, self.ln_kiz))

    def km_at(self, en_td):
        x = np.log(np.clip(en_td, math.exp(self.ln_en[0]), math.exp(self.ln_en[-1])))
        return np.interp(x, self.ln_en, self.km)

    def nu_iz(self, e0_v_m, n_gas, f_hz):
        """Ionization frequency at the effective field of an RF amplitude e0."""
        e_rms = np.asarray(e0_v_m) / math.sqrt(2.0)
        nu_m = n_gas * self.km_at(e_rms / n_gas * 1e21)
        w = 2 * math.pi * f_hz
        e_eff = e_rms * nu_m / np.sqrt(nu_m ** 2 + w ** 2)
        return n_gas * self.kiz(e_eff / n_gas * 1e21)


@dataclass
class Shot:
    gas: str = "He"
    p_pa: float = 0.2e6
    t_k: float = 300.0
    e0_v_m: float = 22.3e5            # RF amplitude at the breakdown point when the cavity is filled
    f_hz: float = 2.767e9
    t_fill_s: float = 0.0              # field envelope 1 - exp(-t/t_fill) from t = 0 (0: already filled)
    t_rf_s: float = 2.4e-6             # RF pulse length (window for breakdown)
    v_hf_m3: float = 1e-6              # volume where the field is within ~10 % of the peak
    s_bg_scale: float = 1.0            # multiplier on the natural background ion-pair rate
    nu_loss_s: float = 0.0             # attachment to impurities / diffusion out of V_hf
    n_br: float = 1e13                 # electrons in the avalanche at breakdown (calibrated on 3D runs)
    laser: bool = True
    t_laser_s: float = 100e-9
    n_laser: float = 1e8               # electrons left by the nominal pulse
    k_order: float = 6.0               # d ln N_L / d ln E_L
    sig_el: float = 0.03               # relative laser-energy rms (Nd:YAG Surelite: ~3 %)
    sig_p: float = 0.02                # relative magnetron-power rms; field ~ sqrt(P)
    sig_tl_s: float = 1e-9             # laser trigger jitter relative to the RF envelope
    dt_s: float = 0.05e-9

    @property
    def n_gas(self):
        return self.p_pa / (KB * self.t_k)

    @property
    def s_bg(self):
        scale = HE_OVER_AIR if self.gas == "He" else 1.0
        return S_AIR_STP * scale * self.n_gas / N_ATM * self.s_bg_scale


def _env(t, t_fill):
    return np.ones_like(t) if t_fill <= 0 else 1.0 - np.exp(-np.maximum(t, 0.0) / t_fill)


def simulate(shot: Shot, n_shots=2000, seed=1, rates: Rates | None = None, sources=None):
    """Monte Carlo of breakdown times (s, from RF start); NaN = no breakdown inside the RF pulse.
    sources: subset of {'el','p','tl','yule','bg'} to keep (default all) for the jitter decomposition."""
    rates = rates or Rates(shot.gas)
    rng = np.random.default_rng(seed)
    src = set(sources) if sources is not None else {"el", "p", "tl", "yule", "bg"}
    t = np.arange(0.0, shot.t_rf_s + shot.dt_s, shot.dt_s)
    env = _env(t, shot.t_fill_s)
    out = np.full(n_shots, np.nan)
    first = np.full(n_shots, np.nan)
    for s in range(n_shots):
        e_scale = math.sqrt(max(1.0 + (shot.sig_p * rng.standard_normal() if "p" in src else 0.0), 0.0))
        nu = rates.nu_iz(shot.e0_v_m * e_scale * env, shot.n_gas, shot.f_hz)
        nu_net = np.maximum(nu - shot.nu_loss_s, 0.0)
        surv = np.where(nu > 0, np.clip(1.0 - shot.nu_loss_s / np.maximum(nu, 1e-300), 0.0, 1.0), 0.0)
        g = np.concatenate([[0.0], np.cumsum(0.5 * (nu_net[1:] + nu_net[:-1]) * shot.dt_s)])
        cands = []
        if shot.laser:
            tl = shot.t_laser_s + (shot.sig_tl_s * rng.standard_normal() if "tl" in src else 0.0)
            el = 1.0 + (shot.sig_el * rng.standard_normal() if "el" in src else 0.0)
            n0 = shot.n_laser * max(el, 0.0) ** shot.k_order
            cands.append((tl, n0, False))
        if "bg" in src:
            # thinning of the Poisson process with rate S_bg V_hf surv(t)
            lam = shot.s_bg * shot.v_hf_m3 * surv
            cum = np.concatenate([[0.0], np.cumsum(0.5 * (lam[1:] + lam[:-1]) * shot.dt_s)])
            u = rng.exponential()
            if u < cum[-1]:
                tb = float(np.interp(u, cum, t))
                cands.append((tb, 1.0, True))
                first[s] = tb
        best = np.inf
        for t0, n0, conditioned in cands:
            i0 = min(max(int(t0 / shot.dt_s), 0), len(t) - 1)
            p_s = 1.0 if conditioned else float(surv[i0])  # background arrivals were thinned by surv
            if "yule" in src:
                k = rng.binomial(int(min(n0, 1e15)), p_s) if n0 < 1e6 else n0 * p_s
                if k <= 0:
                    continue
                mean = nu[i0] / max(nu_net[i0], 1e-300)
                w = rng.gamma(k, mean) if k < 1e6 else k * mean
            else:
                w = n0 * p_s
            target = g[i0] + math.log(max(shot.n_br / w, 1.0))
            if target <= g[-1]:
                tb = float(np.interp(target, g, t))
                best = min(best, tb)
        if np.isfinite(best):
            out[s] = best
    return out, first


def stats(tb, t_ref=0.0):
    ok = np.isfinite(tb)
    d = tb[ok] - t_ref
    return dict(p_breakdown=float(ok.mean()), delay_mean_s=float(d.mean()) if ok.any() else None,
                delay_std_s=float(d.std()) if ok.sum() > 1 else None,
                delay_p10_s=float(np.percentile(d, 10)) if ok.any() else None,
                delay_p90_s=float(np.percentile(d, 90)) if ok.any() else None)


def decompose(shot: Shot, n_shots=2000, rates=None):
    """Jitter contributed by each source alone (others frozen at nominal)."""
    rates = rates or Rates(shot.gas)
    res = {}
    for name, srcs in (("all", None), ("laser_energy", ("el",)), ("rf_power", ("p",)),
                       ("trigger", ("tl",)), ("avalanche_statistics", ("yule",)), ("none", ())):
        tb, _ = simulate(shot, n_shots if name != "none" else 4, seed=7, rates=rates,
                         sources=srcs if srcs is not None else None)
        res[name] = stats(tb, shot.t_laser_s if shot.laser else 0.0)
    return res


def main(out=Path(__file__).resolve().parents[2] / "results" / "stochastic_seed"):
    out.mkdir(parents=True, exist_ok=True)
    r = Rates("He")
    summary = {"assumptions": asdict(Shot()), "s_bg_note": "natural ionization, He scaled 0.12 x air at equal N"}
    # experimental anchor (Beilin 2015 EPL / Yatom 2016 PSST): He 0.2 MPa, 22.3 and 16.6 kV/cm, 2.767 GHz
    anchor = {}
    for e_kv in (22.3, 16.6):
        sh = Shot(e0_v_m=e_kv * 1e5)
        anchor[f"{e_kv}kVcm_laser"] = decompose(sh, rates=r)
        rf = Shot(e0_v_m=e_kv * 1e5, laser=False)
        tb, first = simulate(rf, 4000, rates=r)
        anchor[f"{e_kv}kVcm_rfonly"] = dict(stats(tb), mean_wait_first_electron_s=1.0 / (rf.s_bg * rf.v_hf_m3))
    summary["anchor"] = anchor
    # (p, E) map at constant 2.767 GHz, laser-triggered, nominal fluctuations
    psig = np.array([50, 100, 150, 200, 250.0])
    e_kv = np.array([40, 60, 80, 100, 200.0])
    dmean = np.full((len(psig), len(e_kv)), np.nan)
    dstd = np.full_like(dmean, np.nan)
    pbd_rf = np.full_like(dmean, np.nan)
    for i, ps in enumerate(psig):
        for j, ek in enumerate(e_kv):
            p = ps * 6894.757 + 101325.0
            st = stats(*simulate(Shot(p_pa=p, e0_v_m=ek * 1e5, t_laser_s=20e-9, t_rf_s=1e-6), 600, rates=r)[:1],
                       t_ref=20e-9)
            dmean[i, j] = st["delay_mean_s"] if st["delay_mean_s"] is not None else np.nan
            dstd[i, j] = st["delay_std_s"] if st["delay_std_s"] is not None else np.nan
            pbd_rf[i, j] = stats(simulate(Shot(p_pa=p, e0_v_m=ek * 1e5, laser=False, t_rf_s=1e-6), 600,
                                          rates=r)[0])["p_breakdown"]
    summary["map"] = dict(psig=psig.tolist(), E_kVcm=e_kv.tolist(), delay_mean_s=dmean.tolist(),
                          delay_std_s=dstd.tolist(), p_breakdown_rfonly_1us=pbd_rf.tolist())
    (out / "summary.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


if __name__ == "__main__":
    s = main()
    for k, v in s["anchor"].items():
        print(k, json.dumps(v.get("all", v), default=float))
