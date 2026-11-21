# Bounded local PIC/MCC validation

New scripts do not change `he_swarm_mcc.py` or `compare_fluid.py`.

```powershell
Set-Location 'd:\CST projects\Plasma research\simulations\pic'
& D:\mm\wx\python.exe pic_dynamics.py suite
```

Each run gets a UTC/PID directory in `results/dynamics_validation/`. The suite
uses one CPU process/thread at a time, runs analytic acceleration and constant-
cross-section isotropic relaxation tests, a 10,000-step timing pilot, then two
RF periods at 50 and 100 kV/cm peak. The 100 kV/cm checks are half timestep,
double initial macros and an independent seed. Defaults: 2000 initial macros,
2 fs timestep, 32 cells across 1 mm, seed 271828, pure He at
`200*6894.757+101325 = 1480276.4 Pa`, 300 K, 2.95 GHz. The 3.9 fitted field factor
is not used. Initial electrons are an illustrative 0.75 eV Maxwellian seed,
not a simulated laser ionization calculation.

Hard bounds: 2,100,000 total planned steps (default actual 2,046,898), 7200 s
suite timeout, 2400 s per child, 700,000 steps per child, 8000 initial macros
and 20,000 live macros checked at frames. Initial memory estimate 200 MB
per solver process. Pilot timing gates launching longer cases. OMP, MKL and
OpenBLAS thread counts are one. No installation is performed. Per-case options
are shown by `pic_dynamics.py --help`; a suite invocation may override pressure,
seed, particles, timestep, frequency, periods, cells and gap, subject to bounds.
Use new output directories when replaying the saved `commands.json`.

## Software and data

The installed native Windows distribution is pywarpx 26.11 / AMReX 26.11;
the native WarpX banner says `Unknown`, so no unverified git revision is claimed.
Double precision, 1D3V electrostatic multigrid, Boris pusher, linear shape,
energy-conserving gathering, supported background MCC. These are patterned on
the current WarpX capacitive-discharge PICMI/MCC example and existing project
swarm scripts. Exact configurations and `warpx_used_inputs` are retained.
The existing rf_plasma Python environment is used *read-only* for BOLOS report
imports; WarpX Python lacks numba, so no dependency is installed or shared
environment modified.

Cross sections are existing `xs_he` converted LoKI He / LXCat tables: energy
eV, sigma m², elastic momentum-transfer treated isotropically, 42 excitation
thresholds and impact ionization. Per-table SHA256, threshold and range are
recorded. Raw LXCat downloads are not copied or published. Excitation energy
is removed and the electron scattered by WarpX; ionization uses WarpX's built-in
redistribution and creates an electron plus He+. No custom energy redistribution
or complete helium chemistry is claimed. Constant sigma 1e-21 m² is deliberately
artificial and used ONLY in the analytic relaxation benchmark.

## Outputs and contract

Latest compact output: `results/dynamics_validation/pic_dynamics.json`.
An identical immutable run copy, `summary.json`, `diagnostics.png`, `REPORT.md`,
`eedf_frames.png`, `status.json`, `resource_estimate.json`, `commands.json` and each case's
`config.json`, `meta.json`, `case.json`, small compressed `raw.npz` and exact
WarpX inputs are saved under the unique run directory.

Top-level JSON is exactly `{provenance, limits, cases}`. Each target/refinement
case has `{label, pressure_pa, e_peak_v_m, f_hz, dt_s, seed, time_s,
mean_energy_ev, mean_vz_m_s, relative_electron_count, energy_bins_ev,
eedf_probability, particles, validation}`. There are 65 frames, including
t=0. `particles` contains real fixed initial electron-ID samples, at most 200
per frame: `{time_s,z_m,energy_ev,weight,particle_id}`. Fixed IDs allow actual
trajectory tracking; births contribute to the full EEDF, not to these initial
sample trajectories. No interpolated or synthesized particle animation.

EEDF values are bin probability (not probability density); divide by energy
bin width for 1/eV. Normalization uses *all* live electron weight, with under/
overflow separately reported. Weights in this 1D geometry are electron count
per m² transverse area; density and multiplication come from actual summed
weights, not number of animation points. Raw `records` columns:
`time_s, mean_energy_ev, mean_vz_m_s, electron_weight, electron_macros,
ion_macros, ion_weight, overflow, underflow, escaped_e_macros,
escaped_e_weight, escaped_e_energy_ev_per_m2, escaped_i_macros,
escaped_i_weight, escaped_i_energy_ev_per_m2, nu_m_s, nu_iz_s,
excitation_loss_ev_s, elastic_loss_ev_s, max_energy_ev, energy_variance_ev2`.

Acceptance criteria: acceleration velocity/energy relative error <1%;
isotropic relaxation drift error relative to initial drift <10%; last-cycle
energy NRMSE <10%, drift NRMSE <15% for fluid comparison and numerical
refinements. Ionization with <100 actual MCC events is *not validated*; an
exact one-sided 95% Poisson upper bound divided by measured macro-time exposure
is reported, under an effective constant-hazard interpretation. Snapshot
IID energy SE and two-seed descriptive SE are qualified, not robust confidence
intervals. Energy work accounting uses coarse-frame quadrature and expected
EEDF collision sinks, not exact per-event conservation.

This is a local dilute imposed-field kinetic closure test, not global discharge,
sheath or whole-cavity S-parameter validation. Neutrals are a fixed thermal
reservoir. No gas heating/depletion, Coulomb collisions, attachment,
recombination, metastable populations, superelastic/stepwise ionization,
three-body/dimer/molecular chemistry, optical ionization or RF feedback.
No reduced c, mass ratio or density scaling of a dense experimental plasma is
used: the deliberately dilute seed is the swarm hypothesis. Grid and domain
size convergence are explicitly untested within this budget. Two RF periods
are not a long-time avalanche or periodic-equilibrium proof.

Representative geometry (not a device drawing):

```text
absorbing / phi=0       illustrative seed                absorbing / phi=-E(t)*gap
z=0 ------------------ z=0.45..0.55 mm ------------------ z=1 mm
          uniform imposed E_z(t); fixed thermal neutral He background
```

Drift is small compared with electron thermal speed: 2000-macro instantaneous
drift frames can be noise-dominated. A separate OLS fundamental/third-harmonic
phase-response diagnostic reports coefficient uncertainty but never replaces
the real frame histories. Failed instantaneous-drift and coarse energy-work
criteria are retained as warnings; energy agreement alone is not full dynamic
validation.

References:
- https://warpx.readthedocs.io/en/latest/usage/examples/capacitive_discharge/README.html
- https://warpx.readthedocs.io/en/latest/theory/multiphysics/collisions.html
