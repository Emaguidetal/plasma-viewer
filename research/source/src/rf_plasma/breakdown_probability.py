"""Breakdown probability and delay for a laser-seeded RF avalanche in a sealed high-pressure helium cavity.

Geometry (Technion WR-284 cavity experiment): the laser crosses the cavity perpendicular to the RF E
vector, focused by an f = 150 mm lens; the focus can be moved along the beam axis. The laser fires at
t_L = 4.9 us of a 5 us RF pulse, after which the cavity field rings down with amplitude time constant
tau = 2 Q_L / omega.

Model per shot
* Seed: n(r, z) = n_sat (1 - exp(-n_ref (I/I_ref)^k / n_sat)) (power law of effective order k, optional
  saturation density n_sat), I from a Gaussian beam w(z) = w0 sqrt(1 + ((z - s)/z_R)^2); integrated over
  the beam cross-section to dN0/dz.  n_ref = on-axis seed density at the focus for the reference pulse.
* Field along the beam: E(z) = E_peak * profile(z), profile = cos(2 pi z / lambda_g) (beam along the guide,
  standing wave) or cos(pi z / a) (beam across the broad dimension); z = 0 is the field maximum.
* Avalanche (pure He, no attachment unless nu_loss > 0): G(z, t) = int nu_iz(E(z) env(t)) dt from t_L,
  env = 1 until t_rf_end then exp(-(t - t_rf_end)/tau). Breakdown when sum_z dN0(z) exp(G) >= N_br.
* Fluctuations: magnetron power (field ~ sqrt(P)), laser energy, Poisson number of seed electrons and the
  Yule (Furry) avalanche factor W ~ Gamma(K, 1)/K.
e_cal multiplies the field seen by the electrons (calibration knob for the low-field discrepancy).
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np

from .stochastic_seed import KB, Rates

C0 = 299792458.0
PSI = 6894.757
ATM = 101325.0
OUT = Path(__file__).resolve().parents[2] / "results" / "breakdown_probability"


@dataclass
class Setup:
    p_pa: float = 200 * PSI + ATM       # sealed He, 200 psig
    t_k: float = 300.0
    f_hz: float = 2.95e9
    e_peak_v_m: float = 50e5            # peak cavity field at the field maximum on the beam line
    e_cal: float = 1.0
    axis: str = "z"                     # 'z': along the guide (standing wave); 'x': across the broad wall
    a_m: float = 0.07214                # WR-284 broad dimension
    focus_m: float = 0.0                # focus position along the beam relative to the field maximum
    energy_j: float = 1e-3
    fwhm_s: float = 1e-9
    wavelength_m: float = 532e-9
    w0_m: float = 20e-6
    k_order: float = 6.0
    n_ref_m3: float = 1e16              # seed density on axis at the focus for the reference pulse
    e_ref_j: float = 1e-3
    n_sat_m3: float = 0.0               # 0: no saturation
    t_laser_s: float = 4.9e-6
    t_rf_end_s: float = 5.0e-6
    tau_amp_s: float = 35.7e-9          # 2 Q_L / omega, Q_L = 331 (3D cavity calibration)
    n_br: float = 1e12
    nu_loss_s: float = 0.0
    sig_p: float = 0.02
    sig_el: float = 0.03
    z_half_m: float = 0.03
    nz: int = 961
    dt_s: float = 0.2e-9
    n_tau_after: float = 10.0

    @property
    def n_gas(self):
        return self.p_pa / (KB * self.t_k)

    @property
    def z_r(self):
        return math.pi * self.w0_m ** 2 / self.wavelength_m

    def i_peak(self, energy_j):
        p = energy_j * math.sqrt(4 * math.log(2) / math.pi) / self.fwhm_s
        return 2 * p / (math.pi * self.w0_m ** 2)


def field_profile(st: Setup, z):
    if st.axis == "x":
        return np.cos(math.pi * z / st.a_m)
    lam = C0 / st.f_hz
    lam_g = lam / math.sqrt(1 - (lam / (2 * st.a_m)) ** 2)
    return np.cos(2 * math.pi * z / lam_g)


def seed_per_length(st: Setup, energy_j, z, nu_r=64):
    """dN0/dz (m^-1) along the beam for a pulse of energy_j (time-integrated seed, fixed pulse length)."""
    wz = st.w0_m * np.sqrt(1 + ((z - st.focus_m) / st.z_r) ** 2)
    i0 = st.i_peak(energy_j) * (st.w0_m / wz) ** 2
    u = np.linspace(0.0, 3.0, nu_r)
    lin = st.n_ref_m3 * (i0[:, None] / st.i_peak(st.e_ref_j)) ** st.k_order \
        * np.exp(-2 * st.k_order * u[None, :] ** 2)
    n = st.n_sat_m3 * -np.expm1(-lin / st.n_sat_m3) if st.n_sat_m3 > 0 else lin
    return np.trapezoid(2 * math.pi * u[None, :] * n, u, axis=1) * wz ** 2


def _time_axis(st: Setup):
    t_end = st.t_rf_end_s - st.t_laser_s + st.n_tau_after * st.tau_amp_s
    t = np.arange(0.0, t_end + st.dt_s, st.dt_s)
    env = np.where(t + st.t_laser_s <= st.t_rf_end_s, 1.0,
                   np.exp(-(t + st.t_laser_s - st.t_rf_end_s) / max(st.tau_amp_s, 1e-30)))
    return t, env


def simulate(st: Setup, n_shots=1000, seed=1, rates: Rates | None = None):
    """Breakdown delay after the laser (s) per shot; NaN = no breakdown. Returns (delays, N0 mean)."""
    rates = rates or Rates("He")
    rng = np.random.default_rng(seed)
    z = np.linspace(-st.z_half_m, st.z_half_m, st.nz)
    dz = z[1] - z[0]
    t, env = _time_axis(st)
    dn_nom = seed_per_length(st, st.energy_j, z) * dz
    if not dn_nom.max() > 0:
        return np.full(n_shots, np.nan), 0.0
    # avalanche on coarse z bins (the cavity field varies on a cm scale)
    nb = 8
    zb = z[: len(z) // nb * nb].reshape(-1, nb).mean(1)
    prof = field_profile(st, zb)
    # cumulative growth table G(amplitude, t) for relative field amplitudes around the nominal value
    amp = np.exp(np.linspace(math.log(0.5), math.log(1.2), 281))
    nu_t = np.maximum(rates.nu_iz(st.e_peak_v_m * st.e_cal * amp[:, None] * env[None, :], st.n_gas, st.f_hz)
                      - st.nu_loss_s, 0.0)
    g_tab = np.concatenate([np.zeros((len(amp), 1)),
                            np.cumsum(0.5 * (nu_t[:, 1:] + nu_t[:, :-1]) * st.dt_s, 1)], 1)
    la = np.log(amp)
    out = np.full(n_shots, np.nan)
    for s in range(n_shots):
        fs = math.sqrt(max(1.0 + st.sig_p * rng.standard_normal(), 0.0))
        el = max(1.0 + st.sig_el * rng.standard_normal(), 0.0) * st.energy_j
        dnf = seed_per_length(st, el, z) * dz
        dn = dnf[: len(zb) * nb].reshape(-1, nb).sum(1)
        n0 = dn.sum()
        if n0 <= 0:
            continue
        k = rng.poisson(n0) if n0 < 1e7 else n0
        if k == 0:
            continue
        w = rng.gamma(k, 1.0) if k < 1e7 else float(k)
        sel = dn > 1e-12 * n0
        x = np.clip(np.log(fs * prof[sel]), la[0], la[-1])
        i = np.minimum(np.searchsorted(la, x) - 1, len(la) - 2).clip(0)
        f = ((x - la[i]) / (la[i + 1] - la[i]))[:, None]
        g = (1 - f) * g_tab[i] + f * g_tab[i + 1]
        lw = np.log(dn[sel] / n0)[:, None] + g
        m = lw.max(0)
        ln_n = math.log(w) + m + np.log(np.exp(lw - m).sum(0))
        hit = np.nonzero(ln_n >= math.log(st.n_br))[0]
        if hit.size:
            j = hit[0]
            if j == 0:
                out[s] = 0.0
            else:
                a, b = ln_n[j - 1], ln_n[j]
                out[s] = t[j - 1] + (math.log(st.n_br) - a) / (b - a) * st.dt_s
    return out, float(dn_nom.sum())


def stats(d):
    ok = np.isfinite(d)
    return dict(p=float(ok.mean()),
                delay_median_ns=float(np.median(d[ok]) * 1e9) if ok.any() else None,
                delay_p10_ns=float(np.percentile(d[ok], 10) * 1e9) if ok.any() else None,
                delay_p90_ns=float(np.percentile(d[ok], 90) * 1e9) if ok.any() else None)


def run(st: Setup, n_shots=600, rates=None, seed=1):
    d, n0 = simulate(st, n_shots, seed, rates)
    return dict(stats(d), n0_mean=n0)


def calibrate_e_cal(st: Setup, p_target, rates, lo=0.5, hi=20.0, n_shots=600, it=22):
    """Field multiplier giving breakdown probability p_target (bisection in log)."""
    for _ in range(it):
        mid = math.sqrt(lo * hi)
        if run(replace(st, e_cal=mid), n_shots, rates)["p"] >= p_target:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


def scan(base: Setup, rates, n_shots=600):
    e_kv = [10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 125, 150, 200]
    en_mj = [0.5, 0.75, 1.0]
    focus_mm = list(np.arange(-15.0, 15.01, 1.0))
    res = dict(e_kv=e_kv, energy_mJ=en_mj, focus_mm=focus_mm)
    res["vs_field"] = {f"{m}": [run(replace(base, e_peak_v_m=e * 1e5, energy_j=m * 1e-3), n_shots, rates)
                                for e in e_kv] for m in en_mj}
    res["vs_focus"] = {f"{e}": [run(replace(base, e_peak_v_m=e * 1e5, energy_j=0.75e-3, focus_m=f * 1e-3),
                                    n_shots, rates) for f in focus_mm] for e in (20, 30, 50)}
    return res


def main(out=OUT, n_shots=600):
    out.mkdir(parents=True, exist_ok=True)
    r = Rates("He")
    base = Setup()
    summary = dict(assumptions=asdict(base), z_r_mm=base.z_r * 1e3,
                   i_peak_w_cm2=base.i_peak(base.energy_j) * 1e-4,
                   n0_focus={f"{m}mJ": float(seed_per_length(base, m * 1e-3, np.linspace(-0.03, 0.03, 4801)).sum()
                                             * 0.06 / 4800) for m in (0.5, 0.75, 1.0)})
    summary["uncalibrated"] = scan(base, r, n_shots)
    # calibration: e_cal such that 0.5 mJ at the focus breaks down with P = 0.7 at 20 kV/cm
    e_cal = calibrate_e_cal(replace(base, e_peak_v_m=20e5, energy_j=0.5e-3), 0.7, r, n_shots=n_shots)
    cal = replace(base, e_cal=e_cal)
    summary["e_cal"] = e_cal
    summary["nu_iz_ratio_at_20kVcm"] = float(r.nu_iz(20e5 * e_cal, base.n_gas, base.f_hz)
                                             / max(r.nu_iz(20e5, base.n_gas, base.f_hz), 1e-300))
    summary["calibrated"] = scan(cal, r, n_shots)
    sens = {}
    for name, kw in (("tau_10ns", dict(tau_amp_s=10e-9)), ("tau_200ns", dict(tau_amp_s=200e-9)),
                     ("k4", dict(k_order=4.0)), ("k8", dict(k_order=8.0)),
                     ("n_ref_1e12", dict(n_ref_m3=1e12)), ("n_ref_1e20", dict(n_ref_m3=1e20)),
                     ("w0_10um", dict(w0_m=10e-6)), ("w0_30um", dict(w0_m=30e-6)),
                     ("n_br_1e13", dict(n_br=1e13)), ("sig_p_5pct", dict(sig_p=0.05))):
        s2 = replace(cal, **kw)
        sens[name] = dict(
            vs_energy_20kV={f"{m}": run(replace(s2, e_peak_v_m=20e5, energy_j=m * 1e-3), n_shots, r)
                            for m in (0.5, 0.75, 1.0)},
            focus_10mm_20kV=run(replace(s2, e_peak_v_m=20e5, energy_j=0.75e-3, focus_m=10e-3), n_shots, r))
    summary["sensitivity_calibrated"] = sens
    (out / "breakdown_probability.json").write_text(json.dumps(summary, indent=1))
    figures(summary, out)
    return summary


def figures(s, out=OUT):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    for row, key in enumerate(("uncalibrated", "calibrated")):
        d = s[key]
        lab = "rates as computed" if key == "uncalibrated" else f"field x {s['e_cal']:.2f} (fit: 0.5 mJ, 20 kV/cm, P=0.7)"
        for m in d["energy_mJ"]:
            rr = d["vs_field"][f"{m}"]
            ax[row, 0].plot(d["e_kv"], [x["p"] for x in rr], "o-", label=f"{m} mJ")
            ax[row, 1].semilogy(d["e_kv"], [x["delay_median_ns"] or np.nan for x in rr], "o-", label=f"{m} mJ")
        for e, rr in d["vs_focus"].items():
            ax[row, 2].plot(d["focus_mm"], [x["p"] for x in rr], "o-", label=f"{e} kV/cm")
        ax[row, 0].axvspan(20, 50, color="0.9", zorder=0)
        ax[row, 0].set(xscale="log", xlabel="peak cavity field (kV/cm)", ylabel="breakdown probability",
                       title=f"P vs field ({lab})")
        ax[row, 1].set(xscale="log", xlabel="peak cavity field (kV/cm)", ylabel="median delay after laser (ns)",
                       title="delay vs field")
        ax[row, 2].set(xlabel="focus offset along beam (mm)", ylabel="breakdown probability",
                       title="P vs focus position, 0.75 mJ")
        for a in ax[row]:
            a.grid(alpha=0.3)
            a.legend(fontsize=8)
    fig.suptitle("Laser-seeded RF breakdown, He 200 psig (1.48 MPa), laser perpendicular to E, "
                 "RF ends 100 ns after the laser")
    fig.tight_layout()
    fig.savefig(out / "breakdown_probability.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
