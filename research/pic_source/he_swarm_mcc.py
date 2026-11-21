"""WarpX 1D PIC/MCC electron swarm in He (spot check of the fluid EEDF closure).

A uniform, externally applied field is produced by the Dirichlet plate potentials of a 1 mm gap; the
electron density is tiny so space charge is negligible (swarm limit). Electrons start mid-gap and,
at 1 MPa, drift/diffuse only tens of micrometres, so wall losses are zero. All 42 LoKI excitation
levels, elastic (isotropic, momentum-transfer cross section) and ionization are included.

Run with the WarpX env:
    D:\\mm\\wx\\python.exe he_swarm_mcc.py dc  <E/N Td>  [t_ns]
    D:\\mm\\wx\\python.exe he_swarm_mcc.py rf  <E0 V/m> <f Hz> [t_ns]
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from pywarpx import callbacks, picmi

HERE = Path(__file__).resolve().parent
XS = HERE / "xs_he"
QE, ME, C = 1.602176634e-19, 9.1093837015e-31, 299792458.0
M_HE = 4.002602 * 1.66053906660e-27
P_PA, T_GAS = 1.0e6, 300.0
N_GAS = P_PA / (1.380649e-23 * T_GAS)
GAP, NZ = 1.0e-3, 64
DT = 2.0e-15
N_MACRO = 6000
SAMPLE = 250  # steps between diagnostics
EBINS = np.linspace(0.0, 40.0, 161)


def main():
    mode = sys.argv[1]
    if mode == "dc":
        td = float(sys.argv[2])
        e0, freq = td * 1e-21 * N_GAS, 0.0
        t_end = float(sys.argv[3]) * 1e-9 if len(sys.argv) > 3 else 0.3e-9
        volt = f"{-e0 * GAP:.9e}"
        tag = f"dc_{td:g}Td"
    else:
        e0, freq = float(sys.argv[2]), float(sys.argv[3])
        t_end = float(sys.argv[4]) * 1e-9 if len(sys.argv) > 4 else 2.0 / freq + 0.05e-9
        volt = f"{-e0 * GAP:.9e}*sin(2*pi*{freq:.9e}*t)"
        tag = f"rf_{e0 / 1e6:g}MVm_{freq / 1e9:g}GHz"
    out = HERE / "results" / tag
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((XS / "manifest.json").read_text())

    grid = picmi.Cartesian1DGrid(
        number_of_cells=[NZ], lower_bound=[0.0], upper_bound=[GAP], warpx_max_grid_size=NZ,
        lower_boundary_conditions=["dirichlet"], upper_boundary_conditions=["dirichlet"],
        lower_boundary_conditions_particles=["absorbing"], upper_boundary_conditions_particles=["absorbing"],
        warpx_potential_lo_z="0.0", warpx_potential_hi_z=volt)
    solver = picmi.ElectrostaticSolver(grid=grid)
    vth = np.sqrt(QE * 0.5 / ME)  # 0.75 eV initial mean energy
    electrons = picmi.Species(particle_type="electron", name="electrons",
                              initial_distribution=picmi.UniformDistribution(
                                  density=1.0e12, lower_bound=[None, None, 0.45 * GAP],
                                  upper_bound=[None, None, 0.55 * GAP], rms_velocity=[vth] * 3))
    ions = picmi.Species(particle_type="He", name="he_ions", charge="q_e", mass=M_HE)

    procs = {}
    for p in manifest["processes"]:
        entry = {"cross_section": str(XS / p["file"])}
        if p["energy_eV"] is not None:
            entry["energy"] = p["energy_eV"]
        if p["name"] == "ionization":
            entry["species"] = ions
        procs[p["name"]] = entry
    mcc = picmi.MCCCollisions(name="coll_elec", species=electrons, background_density=N_GAS,
                              background_temperature=T_GAS, background_mass=M_HE,
                              scattering_processes=procs)
    steps = int(round(t_end / DT))
    sim = picmi.Simulation(solver=solver, time_step_size=DT, max_steps=steps, warpx_collisions=[mcc],
                           warpx_collisions_split_momentum_push=0, verbose=0)
    sim.add_species(electrons, layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=N_MACRO // 6, grid=grid))
    sim.add_species(ions, layout=None)

    rec = {"t": [], "eps": [], "vz": [], "w": [], "n_macro": [], "zmin": [], "zmax": []}
    hist = {"t": [], "f": []}

    def sample():
        step = sim.extension.warpx.getistep(lev=0)
        if step % SAMPLE:
            return
        pc = sim.particles.get("electrons")
        ux, uy, uz, w, z = [], [], [], [], []
        for pti in pc.iterator(level=0):
            ux.append(np.asarray(pti["ux"])); uy.append(np.asarray(pti["uy"]))
            uz.append(np.asarray(pti["uz"])); w.append(np.asarray(pti["w"])); z.append(np.asarray(pti["z"]))
        ux, uy, uz, w, z = (np.concatenate(a) for a in (ux, uy, uz, w, z))
        u2 = ux**2 + uy**2 + uz**2
        gam = np.sqrt(1.0 + u2 / C**2)
        eps = (gam - 1.0) * ME * C**2 / QE
        vz = uz / gam
        ws = w.sum()
        rec["t"].append(step * DT); rec["eps"].append(float((w * eps).sum() / ws))
        rec["vz"].append(float((w * vz).sum() / ws)); rec["w"].append(float(ws))
        rec["n_macro"].append(int(w.size)); rec["zmin"].append(float(z.min())); rec["zmax"].append(float(z.max()))
        h, _ = np.histogram(eps, bins=EBINS, weights=w)
        hist["t"].append(step * DT); hist["f"].append(h / ws)

    callbacks.installafterstep(sample)
    t0 = time.time()
    sim.step(steps)
    wall = time.time() - t0
    np.savez(out / "swarm.npz", **{k: np.asarray(v) for k, v in rec.items()},
             hist_t=np.asarray(hist["t"]), hist_f=np.asarray(hist["f"]), ebins=EBINS)
    meta = {"mode": mode, "E0_V_m": e0, "f_Hz": freq, "N_m3": N_GAS, "E0_over_N_Td": e0 / N_GAS * 1e21,
            "dt_s": DT, "steps": steps, "t_end_s": t_end, "n_macro0": rec["n_macro"][0] if rec["n_macro"] else None, "wall_s": wall,
            "gap_m": GAP, "nz": NZ, "processes": len(procs),
            "cross_sections": manifest["source"], "code": "WarpX (pywarpx PICMI, 1D, MCC)"}
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
