"""Seed-electron budget for the laser-triggered RF breakdown in the sealed 200 psig helium WR-284 cavity.

Question: which process supplies the first electrons that the RF field multiplies, and does it reproduce the
observations (Krasik / Leopold experiments, user report)?
  O1  laser alone (RF off): no visible plasma;
  O2  at high RF field the laser energy hardly matters;
  O3  at 20-50 kV/cm the laser energy changes the breakdown probability between ~70 % and 100 %;
  O4  moving the focus along the beam changes the probability from 0 to 100 %;
  O5  RF alone (no laser) breaks down in about 10 % of pulses (reported by the user, sealed He).

Candidates (all per shot, counted inside the high-field part of the beam, |x| < x_hf where the TE101 field
is within 10 % of its maximum):
  * He multiphoton ionization (11 photons at 532 nm);
  * multiphoton ionization of impurities O2, H2O, N2 and a generic organic (outgassing) at ppm level;
  * photodetachment of negative ions (O2-) left by background ionization and by the previous shot;
  * background (cosmic / terrestrial) free electrons already present, and their production during the RF pulse;
  * Fowler-Nordheim field emission at the broad walls (electrons at the wall, not at the laser);
  * Cu photoemission (laser hits the narrow walls at x = 0, a where the TE101 field vanishes -> irrelevant).

Rates: PPT (ppt.py, hydrogen-like l = 0 model applied to molecules -> order-of-magnitude, typically an
overestimate for O2 by 10-100); Langevin recombination (high-pressure limit, mobility-controlled);
three-body attachment e + O2 + He (rate constant order of magnitude, unverified, see K3_O2_HE).
Impurity levels are parametric (He 5.0 specification order: O2 <= 2 ppm, H2O <= 3 ppm, N2 <= 5 ppm,
hydrocarbons <= 0.1-0.5 ppm; a sealed cavity may outgas more H2O and organics).
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np

from .breakdown_probability import Setup
from .ppt import ppt_table
from .stochastic_seed import HE_OVER_AIR, N_ATM, S_AIR_STP, Rates

QE = 1.602176634e-19
EPS0 = 8.8541878128e-12
KB = 1.380649e-23
H_PL = 6.62607015e-34
C0 = 299792458.0
PSI = 6894.757
ATM = 101325.0
A_WG, B_WG, D_CAV = 0.072136, 0.034036, 0.06927
OUT = Path(__file__).resolve().parents[2] / "results" / "seed_sources"

# ionization energies (eV) [NIST Chemistry WebBook]
SPECIES = {"He": 24.587, "N2": 15.581, "H2O": 12.621, "O2": 12.070, "organic (8.8 eV)": 8.83}
K3_O2_HE = 3e-44            # m^6/s, e + O2 + He -> O2- + He (order 1e-32 cm^6/s; unverified, x3)
SIGMA_PD_O2M = 8e-23        # m^2, O2- photodetachment near 2.3 eV (order 1e-18 cm^2; unverified, x3)
MU_ION = 1.255e-4           # m^2/Vs, He ion mobility at 1.48 MPa (streamer2d table); used for all ions
P_RF_ONLY_MEASURED = 0.10   # RF-only breakdown probability per pulse, user report (O5)
FN_A, FN_B = 1.54e-6, 6.83e9  # Fowler-Nordheim constants (A eV / V^2, eV^-3/2 V/m)


@dataclass
class Budget:
    p_pa: float = 200 * PSI + ATM
    t_k: float = 300.0
    f_hz: float = 2.95e9
    energy_j: float = 0.75e-3
    focus_m: float = 0.0
    ppm: dict = field(default_factory=lambda: {"He": 1e6, "N2": 5.0, "H2O": 3.0, "O2": 2.0,
                                              "organic (8.8 eV)": 0.1})
    shot_interval_s: float = 1.0        # unknown; repetition of the experiment
    n_prev_m3: float = 1e23             # channel density left by the previous shot
    v_prev_m3: float = 1e-9             # volume of the previous channel (~ 30 mm x 0.2 mm diameter)
    t_rf_s: float = 5.0e-6
    n_br: float = 1e12
    beta_fn: float = 100.0              # field enhancement of wall micro-protrusions
    a_fn_m2: float = 1e-12              # total emitting area
    work_fn_ev: float = 4.65            # Cu
    x_hf_m: float = 0.0                 # 0: |x| where cos(pi x / a) = 0.9

    @property
    def n_gas(self):
        return self.p_pa / (KB * self.t_k)

    @property
    def x_hf(self):
        return self.x_hf_m or A_WG / math.pi * math.acos(0.9)

    def setup(self):
        return Setup(p_pa=self.p_pa, axis="x", focus_m=self.focus_m, energy_j=self.energy_j)


# ----------------------------------------------------------------------------- laser processes
_TABLES = {}


def _rate_table(ip):
    if ip not in _TABLES:
        _TABLES[ip] = ppt_table(532e-9, ip, 1.0, 0, i_min_w_m2=1e10, i_max_w_m2=1e18, points=161)
    return _TABLES[ip]


def mpi_per_length(st: Setup, ip_ev, n_species, x, energy_j, nr=48, nt=81):
    """Electrons per metre along the beam from multiphoton ionization of a species with density
    n_species; Gaussian beam in space and time (1 - exp(-int W dt) saturation kept)."""
    ln_i, ln_w = _rate_table(ip_ev)
    wz = st.w0_m * np.sqrt(1 + ((x - st.focus_m) / st.z_r) ** 2)
    i0 = st.i_peak(energy_j) * (st.w0_m / wz) ** 2
    tau = st.fwhm_s / (2 * math.sqrt(math.log(2)))      # I ~ exp(-(t/tau)^2)
    t = np.linspace(-3, 3, nt) * tau
    u = np.linspace(0, 2.0, nr)
    inten = i0[:, None, None] * np.exp(-2 * u[None, :, None] ** 2) * np.exp(-(t[None, None, :] / tau) ** 2)
    lw = np.interp(np.log(np.maximum(inten, 1e-30)), ln_i, ln_w, left=-700.0)
    w = np.exp(np.where(inten < math.exp(ln_i[0]), -700.0, lw))
    prob = -np.expm1(-np.trapezoid(w, t, axis=2))
    return n_species * np.trapezoid(2 * math.pi * u[None, :] * prob, u, axis=1) * wz ** 2


def photodetach_per_length(st: Setup, n_neg, x, energy_j, sigma=SIGMA_PD_O2M, nr=64):
    """Electrons per metre from one-photon detachment of a uniform negative-ion density."""
    wz = st.w0_m * np.sqrt(1 + ((x - st.focus_m) / st.z_r) ** 2)
    phot = energy_j / (H_PL * C0 / 532e-9)
    f0 = 2 * phot / (math.pi * wz ** 2)                  # peak photon fluence
    u = np.linspace(0, 3.0, nr)
    prob = -np.expm1(-sigma * f0[:, None] * np.exp(-2 * u[None, :] ** 2))
    return n_neg * np.trapezoid(2 * math.pi * u[None, :] * prob, u, axis=1) * wz ** 2


# ----------------------------------------------------------------------------- background population
def background(b: Budget):
    n = b.n_gas
    s = S_AIR_STP * HE_OVER_AIR * n / N_ATM                       # ion pairs m^-3 s^-1
    lam2 = 1.0 / ((math.pi / A_WG) ** 2 + (math.pi / B_WG) ** 2 + (math.pi / D_CAV) ** 2)
    mu_e = 2.163e25 / n                                          # thermal muN (BOLOS He table, low E/N)
    kt = KB * b.t_k / QE
    tau_de = lam2 / (mu_e * kt)
    tau_di = lam2 / (2 * MU_ION * kt)                            # ambipolar ~ 2 D_i
    nu_att = K3_O2_HE * n * n * b.ppm["O2"] * 1e-6
    beta_ii = QE * 2 * MU_ION / EPS0                              # Langevin ion-ion
    alpha_ei = QE * (mu_e + MU_ION) / EPS0                       # Langevin e-ion (upper bound)
    # free electrons: production s, loss attachment + diffusion + recombination with n_+
    nu_e = nu_att + 1 / tau_de
    # ions: s = n/tau_di + beta n^2  (negative ions dominate the negative charge when nu_att >> 1/tau_de)
    a_, b_, c_ = beta_ii, 1 / tau_di, -s
    n_ion = (-b_ + math.sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_)
    n_e = s / (nu_e + alpha_ei * n_ion)
    frac_neg = nu_att / nu_e
    # previous shot: channel ions spread by diffusion and recombine (Langevin ion-ion after attachment)
    d_i = MU_ION * kt
    r_spread = math.sqrt(4 * d_i * b.shot_interval_s)
    v_spread = max(b.v_prev_m3, 2 * math.pi * r_spread ** 2 * 0.03)  # 30 mm rod line diffused radially
    n_rec = 1 / (1 / b.n_prev_m3 + beta_ii * b.shot_interval_s)
    n_prev = min(n_rec, b.n_prev_m3 * b.v_prev_m3 / v_spread)
    return dict(s_pairs_m3_s=s, lambda_diff_m=math.sqrt(lam2), mu_e_thermal=mu_e, tau_diff_e_s=tau_de,
                tau_diff_ion_s=tau_di, nu_attach_s=nu_att, tau_attach_s=1 / nu_att if nu_att else math.inf,
                beta_ion_ion_m3_s=beta_ii, alpha_e_ion_langevin_m3_s=alpha_ei, n_ion_m3=n_ion, n_e_m3=n_e,
                frac_electrons_attached=frac_neg, n_neg_bg_m3=n_ion * frac_neg,
                prev_shot_r_spread_m=r_spread, n_neg_prev_shot_m3=n_prev * frac_neg)


def fowler_nordheim(b: Budget, e_wall):
    f = b.beta_fn * np.asarray(e_wall, float)
    phi = b.work_fn_ev
    j = FN_A * f ** 2 / phi * np.exp(-FN_B * phi ** 1.5 / np.maximum(f, 1.0))
    # RF-averaged over the half period with E > 0 (crude: factor 0.25 for the peak-dominated emission)
    return 0.25 * j * b.a_fn_m2 * b.t_rf_s / QE


# ----------------------------------------------------------------------------- RF-only check
def rf_only(b: Budget, e0_v_m, e_cal=1.0, n_free_m3=None, s_m3_s=None, nx=145, nz=139):
    """Expected number of successful RF-only lineages in a TE101 cavity during the RF pulse:
    pre-existing free electrons (n_free) plus production s during the pulse; a lineage needs
    nu_iz (t_rf - t_birth) >= ln N_br. Pure helium (attachment ~ 1e3-1e4 s^-1 << nu_iz)."""
    rates = Rates("He")
    x = (np.arange(nx) + 0.5) / nx * A_WG
    z = (np.arange(nz) + 0.5) / nz * D_CAV
    shape = np.sin(math.pi * x / A_WG)[:, None] * np.sin(math.pi * z / D_CAV)[None, :]
    dv = A_WG / nx * D_CAV / nz * B_WG
    out = []
    for e0 in np.atleast_1d(e0_v_m):
        nu = rates.nu_iz(e_cal * e0 * shape, b.n_gas, b.f_hz)
        t_need = math.log(b.n_br) / np.maximum(nu, 1e-30)
        pre = (t_need <= b.t_rf_s) * n_free_m3
        prod = s_m3_s * np.clip(b.t_rf_s - t_need, 0, None)
        out.append(float(((pre + prod) * dv).sum()))
    return np.array(out)


# ----------------------------------------------------------------------------- budget
def budget(b: Budget, nx=801):
    st = b.setup()
    xh = b.x_hf
    x = np.linspace(-xh, xh, nx)
    bg = background(b)
    rows = []
    for sp, ip in SPECIES.items():
        n_sp = b.n_gas * b.ppm[sp] * 1e-6
        if n_sp <= 0:
            continue
        rows.append(dict(source=f"MPI {sp}", ip_ev=ip, photons=int(math.ceil(ip / 2.3305)),
                         n_species_m3=n_sp,
                         electrons=float(np.trapezoid(mpi_per_length(st, ip, n_sp, x, b.energy_j), x))))
    for label, nn in (("photodetachment O2- (background)", bg["n_neg_bg_m3"]),
                      ("photodetachment O2- (previous shot)", bg["n_neg_prev_shot_m3"])):
        rows.append(dict(source=label, n_species_m3=nn,
                         electrons=float(np.trapezoid(photodetach_per_length(st, nn, x, b.energy_j), x))))
    beam_vol = math.pi * (3 * st.w0_m) ** 2 * 2 * xh
    rows.append(dict(source="free background electrons (laser-independent)", n_species_m3=bg["n_e_m3"],
                     electrons=bg["n_e_m3"] * beam_vol, note="counted in the beam volume only"))
    return dict(budget=asdict(b), background=bg, x_hf_m=xh, rows=rows)


def seed_line_density(b: Budget, x):
    """Laser seed per metre along the beam used by rod3d (seed_model='organic'): 4-photon MPI of the organic
    impurity plus photodetachment of the background O2- (the previous-shot reservoir is excluded: its
    saturated, focus-independent yield would contradict the measured focus dependence)."""
    st = b.setup()
    n_org = b.n_gas * b.ppm["organic (8.8 eV)"] * 1e-6
    return (mpi_per_length(st, SPECIES["organic (8.8 eV)"], n_org, x, b.energy_j)
            + photodetach_per_length(st, background(b)["n_neg_bg_m3"], x, b.energy_j))


def energy_slope(b: Budget, source_fn, e1=0.5e-3, e2=1.0e-3):
    n1 = source_fn(replace(b, energy_j=e1))
    n2 = source_fn(replace(b, energy_j=e2))
    return math.log(max(n2, 1e-300) / max(n1, 1e-300)) / math.log(e2 / e1)


def per_source(b: Budget):
    """Electrons in |x| < x_hf for each laser-dependent source (dict label -> callable(Budget))."""
    bg = background(b)
    fns = {}
    for sp, ip in SPECIES.items():
        def f(bb, ip=ip, sp=sp):
            st = bb.setup(); x = np.linspace(-bb.x_hf, bb.x_hf, 801)
            return float(np.trapezoid(mpi_per_length(st, ip, bb.n_gas * bb.ppm[sp] * 1e-6, x, bb.energy_j), x))
        fns[f"MPI {sp}"] = f
    for label, key in (("O2- photodetach. (bg)", "n_neg_bg_m3"), ("O2- photodetach. (prev. shot)",
                                                                    "n_neg_prev_shot_m3")):
        def g(bb, key=key):
            st = bb.setup(); x = np.linspace(-bb.x_hf, bb.x_hf, 801)
            return float(np.trapezoid(photodetach_per_length(st, background(bb)[key], x, bb.energy_j), x))
        fns[label] = g
    return fns


def main(out=OUT):
    out.mkdir(parents=True, exist_ok=True)
    b = Budget()
    res = budget(b)
    fns = per_source(b)
    energies = np.array([0.25, 0.5, 0.75, 1.0, 1.5]) * 1e-3
    focus = np.linspace(0, 30e-3, 13)
    res["energy_scan_j"] = energies.tolist()
    res["focus_scan_m"] = focus.tolist()
    res["energy_scan"] = {k: [f(replace(b, energy_j=e)) for e in energies] for k, f in fns.items()}
    res["focus_scan"] = {k: [f(replace(b, focus_m=s)) for s in focus] for k, f in fns.items()}
    res["energy_slope_0p5_1mJ"] = {k: energy_slope(b, f) for k, f in fns.items()}
    # RF-only
    bg = res["background"]
    e0 = np.array([20, 30, 50, 75, 100, 150, 200, 300]) * 1e5
    res["rf_only"] = dict(
        e0_kv_cm=(e0 / 1e5).tolist(),
        lineages_nominal=rf_only(b, e0, 1.0, bg["n_e_m3"], bg["s_pairs_m3_s"]).tolist(),
        lineages_ecal=rf_only(b, e0, 3.91, bg["n_e_m3"], bg["s_pairs_m3_s"]).tolist(),
        lineages_negative_ions_detached=rf_only(b, e0, 1.0, bg["n_ion_m3"], bg["s_pairs_m3_s"]).tolist(),
        fn_electrons_per_pulse=fowler_nordheim(b, e0).tolist())
    # O5: lineages scale linearly with the background ion-pair rate, so the measured probability fixes that rate
    lam_meas = -math.log(1 - P_RF_ONLY_MEASURED)
    res["rf_only"]["measured_probability"] = P_RF_ONLY_MEASURED
    res["rf_only"]["source_scale_for_measured_ecal"] = [lam_meas / max(v, 1e-300) for v in res["rf_only"]["lineages_ecal"]]
    (out / "seed_sources.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    figures(res, out)
    return res


def figures(res, out=OUT):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    rows = res["rows"]
    lab = [r["source"] for r in rows]
    val = [max(r["electrons"], 1e-30) for r in rows]
    a = ax[0, 0]
    a.barh(range(len(rows)), np.log10(val), color=["C0"] * 5 + ["C1", "C1", "C2"][: len(rows) - 5])
    a.set_yticks(range(len(rows))); a.set_yticklabels(lab, fontsize=8)
    a.axvline(0, color="k", lw=0.8)
    a.set_xlabel("log10 electrons per shot in |x| < x_hf (0.75 mJ, focus at field max)")
    a.set_title("Seed budget (impurity ppm: O2 2, H2O 3, N2 5, organic 0.1)")
    a.set_xlim(-30, 8)
    e = np.array(res["energy_scan_j"]) * 1e3
    a = ax[0, 1]
    for k, v in res["energy_scan"].items():
        v = np.maximum(v, 1e-30)
        if v.max() > 1e-12:
            a.loglog(e, v, "o-", label=f"{k} (k={res['energy_slope_0p5_1mJ'][k]:.1f})")
    a.set_xlabel("laser energy (mJ)"); a.set_ylabel("electrons in high-field region"); a.legend(fontsize=7)
    a.set_title("Energy dependence (O2/O3)")
    s = np.array(res["focus_scan_m"]) * 1e3
    a = ax[1, 0]
    for k, v in res["focus_scan"].items():
        v = np.asarray(v)
        if v.max() > 1e-12:
            a.semilogy(s, v / v[0], "o-", label=k)
    a.set_xlabel("focus offset from field maximum along the beam (mm)")
    a.set_ylabel("electrons in |x| < x_hf / value at 0"); a.legend(fontsize=7)
    a.set_title("Focus dependence (O4)")
    r = res["rf_only"]
    a = ax[1, 1]
    for k, l in (("lineages_nominal", "free e- only, BOLOS swarm"), ("lineages_ecal", "field x 3.91 (e_cal)"),
                 ("lineages_negative_ions_detached", "if all O2- were field-detached")):
        a.semilogy(r["e0_kv_cm"], 1 - np.exp(-np.maximum(r[k], 1e-30)), "o-", label=l)
    a.axhline(P_RF_ONLY_MEASURED, color="k", ls="--", label="measured, sealed He (~10 %)")
    a.axvspan(20, 50, color="0.85", zorder=0, label="experimental operating range")
    a.set_xlabel("cavity peak field E0 (kV/cm)"); a.set_ylabel("RF-only breakdown probability per 5 us pulse")
    a.set_ylim(1e-6, 2); a.legend(fontsize=7); a.set_title("RF-only check (O5)")
    fig.tight_layout()
    fig.savefig(out / "seed_sources.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    r = main()
    for row in r["rows"]:
        print(f"{row['source']:45s} {row['electrons']:.3e}")
    print(json.dumps(r["background"], indent=1, default=float))
    print(json.dumps(r["energy_slope_0p5_1mJ"], indent=1))
    print(json.dumps(r["rf_only"], indent=1))
