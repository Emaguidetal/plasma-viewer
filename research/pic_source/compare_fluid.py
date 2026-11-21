"""Compare WarpX PIC/MCC He swarms with the fluid (BOLOS-table, local-mean-energy) closure used by
the coupled solver in the swarm limit (ne -> 0, so the Maxwellian blend weight is zero).

Run with the rf_plasma venv after he_swarm_mcc.py has produced results/*/swarm.npz:
    ..\\rf_plasma\\.venv\\Scripts\\python.exe compare_fluid.py
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "rf_plasma" / "src"))
from rf_plasma.gas_model import load_gas_model  # noqa: E402
from rf_plasma.plasma3d import P_IP  # noqa: E402

QE, ME = 1.602176634e-19, 9.1093837015e-31
RES = HERE / "results"
gm = load_gas_model("He", 800e-9)
EPS = gm.eps_grid
KM, KIZ, KEXC, KEL = gm.tabs_bulk
IP = float(gm.prm[P_IP])
TWO_M_OVER_M = 2 * ME / (4.002602 * 1.66053906660e-27)


def tab(row, e):
    return np.interp(np.log(np.maximum(e, EPS[0])), np.log(EPS), row)


def fluid_dc(td, n):
    """Mean energy from the BOLOS E/N(eps) table (what the fluid's local-field tables encode)."""
    ok = np.isfinite(gm.en_td_of_eps)
    eps = float(np.interp(td, gm.en_td_of_eps[ok], EPS[ok]))
    e_field = td * 1e-21 * n
    return {"eps_eV": eps, "vd_m_s": QE * e_field / (ME * n * tab(KM, eps)),
            "nu_iz_s": n * tab(KIZ, eps)}


def fluid_rhs(t, y, e0, w, n, eps_g):
    """Swarm-limit fluid: Drude momentum + local mean energy with the kinetic_cell loss terms."""
    v, eps = y
    eps = max(eps, 1e-3)
    nu = n * tab(KM, eps)
    e_t = e0 * np.sin(w * t)
    dv = -QE * e_t / ME - nu * v
    gain = -e_t * v  # eV/s per electron (QE cancels)
    kiz = tab(KIZ, eps)
    loss = n * tab(KEXC, eps) + n * tab(KEL, eps) * (eps - eps_g) + n * kiz * (IP + eps)
    return [dv, gain - loss]


def ion_rate(t, w, t0):
    m = t >= t0
    return float(np.polyfit(t[m], np.log(w[m]), 1)[0]) if m.sum() > 3 else float("nan")


def main():
    out = {"fluid_tables": "BOLOS two-term (LoKI He), Boltzmann rows of GasModel.tabs_bulk",
           "pic": "WarpX 1D MCC, 42 LoKI excitations + elastic + ionization, swarm limit", "dc": [], "rf": []}
    dc_dirs = sorted((p for p in RES.glob("dc_*Td") if (p / "meta.json").exists()), key=lambda p: float(p.name[3:-2]))
    rf_dirs = sorted(p for p in RES.glob("rf_*") if (p / "meta.json").exists())
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    for d in dc_dirs:
        meta = json.loads((d / "meta.json").read_text())
        s = np.load(d / "swarm.npz")
        t, n = s["t"], meta["N_m3"]
        t0 = 0.4 * t[-1]
        m = t >= t0
        pic = {"eps_eV": float(s["eps"][m].mean()), "vd_m_s": float(-s["vz"][m].mean()),
               "nu_iz_s": ion_rate(t, s["w"], t0)}
        n_ev = (s["w"][-1] - s["w"][m][0]) / (s["w"][0] / s["n_macro"][0])
        pic["ionization_events"] = float(n_ev)
        pic["nu_iz_rel_stat_err"] = float(1 / np.sqrt(max(n_ev, 1.0)))
        fl = fluid_dc(meta["E0_over_N_Td"], n)
        f = s["hist_f"][s["hist_t"] >= t0].mean(axis=0)
        eb = s["ebins"]
        ec = 0.5 * (eb[1:] + eb[:-1])
        de = np.diff(eb)
        out["dc"].append({"E_N_Td": meta["E0_over_N_Td"], "pic": pic, "fluid": fl,
                          "ratio_eps": pic["eps_eV"] / fl["eps_eV"], "ratio_vd": pic["vd_m_s"] / fl["vd_m_s"],
                          "ratio_nu_iz": pic["nu_iz_s"] / fl["nu_iz_s"] if fl["nu_iz_s"] > 0 else None,
                          "frac_above_19_8eV": float(f[ec > 19.8].sum()), "wall_s": meta["wall_s"]})
        eedf = f / de / np.sqrt(np.maximum(ec, 1e-9))  # eV^-3/2
        te = 2 * pic["eps_eV"] / 3
        mx = 2 / np.sqrt(np.pi) * te**-1.5 * np.exp(-ec / te)
        l, = ax[0, 2].semilogy(ec, eedf, label=f"PIC {meta['E0_over_N_Td']:.1f} Td")
        ax[0, 2].semilogy(ec, mx, "--", color=l.get_color(), lw=0.8)
    if out["dc"]:
        td = [r["E_N_Td"] for r in out["dc"]]
        ax[0, 0].plot(td, [r["pic"]["eps_eV"] for r in out["dc"]], "o", label="WarpX PIC/MCC")
        ax[0, 0].plot(td, [r["fluid"]["eps_eV"] for r in out["dc"]], "s-", mfc="none", label="fluid table (BOLOS)")
        ax[0, 0].set(xlabel="E/N (Td)", ylabel="mean energy (eV)", title="DC swarm: mean energy")
        ax[0, 1].plot(td, [r["pic"]["vd_m_s"] for r in out["dc"]], "o", label="PIC drift")
        ax[0, 1].plot(td, [r["fluid"]["vd_m_s"] for r in out["dc"]], "s-", mfc="none", label="fluid eE/(m N k_m)")
        ax2 = ax[0, 1].twinx()
        ax2.semilogy(td, [max(r["pic"]["nu_iz_s"], 1e3) for r in out["dc"]], "^", color="C3", label="PIC nu_iz")
        ax2.semilogy(td, [max(r["fluid"]["nu_iz_s"], 1e3) for r in out["dc"]], "v-", mfc="none", color="C3",
                     label="fluid N k_iz")
        ax2.set_ylabel("ionization frequency (1/s)", color="C3")
        ax[0, 1].set(xlabel="E/N (Td)", ylabel="drift velocity (m/s)", title="DC swarm: drift and ionization")
        ax[0, 0].legend(); ax[0, 1].legend(loc="upper left"); ax2.legend(loc="lower right")
    ax[0, 2].set(xlabel="energy (eV)", ylabel="EEDF (eV$^{-3/2}$)", title="PIC EEDF (dashed: Maxwellian, same <eps>)",
                 ylim=(1e-7, 1))
    ax[0, 2].legend(fontsize=8)
    for k, d in enumerate(rf_dirs):
        meta = json.loads((d / "meta.json").read_text())
        s = np.load(d / "swarm.npz")
        t, n, e0, f = s["t"], meta["N_m3"], meta["E0_V_m"], meta["f_Hz"]
        w = 2 * np.pi * f
        sol = solve_ivp(fluid_rhs, (0, t[-1]), [0.0, 0.75], args=(e0, w, n, 1.5 * 1.380649e-23 * 300 / QE),
                        method="LSODA", t_eval=t, rtol=1e-7, atol=[1e-3, 1e-6], max_step=2e-13)
        eps_f = sol.y[1]
        per = 1 / f
        t0 = t[-1] - per
        m = t >= t0
        nu_f = n * tab(KIZ, eps_f)
        pic = {"eps_cycle_mean_eV": float(s["eps"][m].mean()), "eps_max_eV": float(s["eps"][m].max()),
               "eps_min_eV": float(s["eps"][m].min()), "nu_iz_cycle_s": ion_rate(t, s["w"], t0)}
        fl = {"eps_cycle_mean_eV": float(eps_f[m].mean()), "eps_max_eV": float(eps_f[m].max()),
              "eps_min_eV": float(eps_f[m].min()), "nu_iz_cycle_s": float(nu_f[m].mean())}
        out["rf"].append({"E0_V_m": e0, "f_Hz": f, "E0_over_N_Td": meta["E0_over_N_Td"], "pic": pic, "fluid": fl,
                          "ratio_eps_mean": pic["eps_cycle_mean_eV"] / fl["eps_cycle_mean_eV"],
                          "ratio_nu_iz": pic["nu_iz_cycle_s"] / fl["nu_iz_cycle_s"] if fl["nu_iz_cycle_s"] > 0
                          else None, "wall_s": meta["wall_s"]})
        c = f"C{k}"
        ax[1, 0].plot(t * 1e9, s["eps"], color=c, label=f"PIC {e0 / 1e6:g} MV/m")
        ax[1, 0].plot(t * 1e9, eps_f, "--", color=c, label=f"fluid {e0 / 1e6:g} MV/m")
        ax[1, 1].semilogy(t * 1e9, s["w"] / s["w"][0], color=c, label=f"PIC {e0 / 1e6:g} MV/m")
        nf = np.exp(np.concatenate([[0], np.cumsum(0.5 * (nu_f[1:] + nu_f[:-1]) * np.diff(t))]))
        ax[1, 1].semilogy(t * 1e9, nf, "--", color=c, label=f"fluid {e0 / 1e6:g} MV/m")
        ax[1, 2].plot(t * 1e9, -s["vz"], color=c, label=f"PIC {e0 / 1e6:g} MV/m")
        ax[1, 2].plot(t * 1e9, sol.y[0] * -1, "--", color=c, label=f"fluid {e0 / 1e6:g} MV/m")
    ax[1, 0].set(xlabel="t (ns)", ylabel="mean energy (eV)", title="3 GHz RF swarm, 1 MPa He: mean energy")
    ax[1, 1].set(xlabel="t (ns)", ylabel="n_e(t)/n_e(0)", title="RF swarm: electron multiplication")
    ax[1, 2].set(xlabel="t (ns)", ylabel="-<v_z> (m/s)", title="RF swarm: drift velocity")
    for a in ax[1]:
        a.legend(fontsize=8)
    fig.suptitle("PIC spot check: WarpX MCC vs fluid EEDF closure (He, 1 MPa, 300 K, swarm limit)")
    fig.tight_layout()
    fig.savefig(RES / "pic_vs_fluid.png", dpi=130)
    (RES / "comparison.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
