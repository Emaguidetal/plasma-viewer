"""Write the scan plans (one scan_runner argument line per 3D run) for scan_queue.py.

t_end from the surrogate delay ln(N_br/N_seed)/nu_iz (stochastic_seed.Rates, field = amplitude); laser cases whose
expected seed is < 1 electron are not run in 3D (a deterministic fluid cannot represent them; see stochastic_seed).
Seed model: laser_seed_m3 at the reference pulse (1 mJ, 1 ns, 20 um), effective order k = 6, normalized so the
reference pulse leaves N_L ~ 1e8 electrons, as in the stochastic surrogate.
"""
import json
import math
from pathlib import Path

from rf_plasma.stochastic_seed import Rates

HERE = Path(__file__).resolve().parent
KB, PSI = 1.380649e-23, 6894.757
N_BR = 4.6e15             # electrons at gamma > 0.3 in the 1 MPa / 80 kV/cm 3D run (fixes_ab)
K, N_REF = 6.0, 1e8
W_REF, T_REF, E_REF = 20e-6, 1e-9, 1e-3
Z_EFF = math.sqrt(math.pi) * math.gamma(K - 1.5) / math.gamma(K - 1)   # int (1+u^2)^-(k-1) du
ZR = math.pi * W_REF ** 2 / 532e-9
SEED_M3 = N_REF / (math.pi * W_REF ** 2 / (2 * K) * Z_EFF * ZR)
NU_3D_OVER_SURR = 1.4    # 3D (with stepwise ionization) grew 1.3e9 vs surrogate 0.93e9 at 1 MPa / 80 kV/cm
r = Rates()


def n_seed(e_mj, fwhm_ns, w_um, k=K):
    """Electrons left by the pulse (focus-integrated) relative to the reference normalization."""
    i_rat = (e_mj * 1e-3 / E_REF) / (fwhm_ns * 1e-9 / T_REF) / (w_um * 1e-6 / W_REF) ** 2
    zr = math.pi * (w_um * 1e-6) ** 2 / 532e-9
    zf = math.sqrt(math.pi) * math.gamma(k - 1.5) / math.gamma(k - 1)
    return SEED_M3 * i_rat ** k * (fwhm_ns * 1e-9 / T_REF) * math.pi * (w_um * 1e-6) ** 2 / (2 * k) * zf * zr


def t_end(p_pa, e_kv, ns, lo=14.0, hi=60.0):
    nu = NU_3D_OVER_SURR * float(r.nu_iz(e_kv * 1e5, p_pa / (KB * 293.0), 3e9))
    td = math.log(max(N_BR / max(ns, 1.0), 10.0)) / max(nu, 1e-30) * 1e9
    return nu, td, min(max(5.0 + 1.5 * td + 4.0, lo), hi)


def line(root, extra, te, geom="guide"):
    return (f"--root {root} --geom {geom} --physics fullee_ff --seed 0 --laser_seed {SEED_M3:.3e} "
            f"--t_end_ns {te:.0f} --threads 2 {extra}")


def main():
    plans, meta = {}, {"seed_m3_ref": SEED_M3, "n_br": N_BR, "k": K, "cases": []}
    # pressure x field, baseline laser (reference pulse)
    ns0 = n_seed(1, 1, 20)
    pe = []
    for ps in (50, 100, 150, 200, 250):
        p = ps * PSI + 101325.0
        for ek in (40, 60, 80, 100, 200):
            nu, td, te = t_end(p, ek, ns0)
            run = nu >= 2.4e8 or (ps, ek) in ((100, 40), (200, 60), (250, 80))
            if nu < 2.4e8:
                te = 20.0  # below-threshold check: no breakdown expected inside the window
            meta["cases"].append(dict(scan="pE", psig=ps, E=ek, nu=nu, td_ns=td, t_end=te, run3d=run, n_seed=ns0))
            if run:
                pe.append(line("scan_pE", f"--p_psig {ps} --E_kVcm {ek} --laser --E_mJ 1 --fwhm_ns 1 --waist_um 20",
                               te))
    plans["pE"] = pe
    # laser focus / pulse at 100 psig, 80 kV/cm (baseline 1 mJ / 1 ns / 20 um is in scan_pE)
    p = 100 * PSI + 101325.0
    las = []
    for e_mj, fw, wu, k in ([(1, f, w, K) for w in (20, 50, 100, 200) for f in (0.5, 1, 2)]
                            + [(0.5, 1, 20, K), (0.5, 1, 50, K), (1, 1, 100, 3.0)]):
        ns = n_seed(e_mj, fw, wu, k)
        nu, td, te = t_end(p, 80, ns)
        run = ns >= 1.0 and not (e_mj, fw, wu, k) == (1, 1, 20, K)
        meta["cases"].append(dict(scan="laser", E_mJ=e_mj, fwhm_ns=fw, waist_um=wu, k=k, n_seed=ns, td_ns=td,
                                  t_end=te, run3d=run))
        if run:
            extra = f"--p_psig 100 --E_kVcm 80 --laser --E_mJ {e_mj} --fwhm_ns {fw} --waist_um {wu}"
            if k != K:
                extra += f" --seed_k {k:g} --tag p100psig_E80kVcm_fullee_ff_L{e_mj}mJ{fw}ns_w{wu}_k{k:g}"
            las.append(line("scan_laser", extra, te))
    plans["laser"] = las
    # experimental anchor (Technion He 0.2 MPa, 22.3 / 16.6 kV/cm amplitude; 10 mJ, 5 ns; focus assumed 20 um)
    anc = []
    for ek in (22.3, 16.6):
        ns = n_seed(10, 5, 20)
        nu, td, te = t_end(0.2e6, ek, ns, hi=90.0)
        meta["cases"].append(dict(scan="anchor", E=ek, n_seed=ns, nu=nu, td_ns=td, t_end=te, run3d=True))
        anc.append(line("anchor", f"--p_MPa 0.2 --E_kVcm {ek} --laser --E_mJ 10 --fwhm_ns 5 --waist_um 20", te,
                        geom="cavity"))
    plans["anchor"] = anc
    # convergence: seed normalization x0.01 / x100, and mesh 0.5 mm (short)
    base = "--p_psig 100 --E_kVcm 80 --laser --E_mJ 1 --fwhm_ns 1 --waist_um 20"
    _, _, te = t_end(p, 80, ns0)
    conv = [line("conv", base + " --tag seedx0.01", te).replace(f"{SEED_M3:.3e}", f"{SEED_M3 * 0.01:.3e}"),
            line("conv", base + " --tag seedx100", te).replace(f"{SEED_M3:.3e}", f"{SEED_M3 * 100:.3e}"),
            line("conv", base + " --dx_mm 0.5 --tag dx0.5", 14).replace("--threads 2", "--threads 4")]
    plans["conv"] = conv
    for k_, v in plans.items():
        (HERE / f"plan_{k_}.txt").write_text("\n".join(v) + "\n")
    (HERE / "scan_plan_meta.json").write_text(json.dumps(meta, indent=1))
    print("seed_m3_ref", f"{SEED_M3:.3e}", {k_: len(v) for k_, v in plans.items()})
    for c in meta["cases"]:
        print({k_: (f"{v:.3g}" if isinstance(v, float) else v) for k_, v in c.items()})


if __name__ == "__main__":
    main()
