---
description: Reproducible openEMS RF-plasma benchmarks in this workspace.
applyTo: 'simulations/rf_plasma/**/*.py'
---

# openEMS RF-plasma learnings

- The native Windows openEMS run may change the process working directory.
  Resolve output, input, and calibration paths to absolute paths before
  invoking the solver, and set `CSXCAD_INSTALL_PATH` only for that process.
- For the bundled openEMS 0.37.0 / CSXCAD 0.7.0 build, the verified prescribed
  Drude benchmark uses a `CSPropLorentzMaterial` with `eps_plasma` as ordinary
  frequency in Hz, zero pole frequency, and `eps_relax=1/nu`. Check its
  frequency convention and a slab solution rather than assuming an arbitrary
  Drude API or phasor sign is equivalent.
- `CSPropConductingSheet` is documented for thin sheets, not automatically a
  validated skin-effect model for a thick 3D copper wall. Keep PEC solves and
  analytical copper estimates distinct; verify any 3D wall treatment against
  an independent attenuation or cavity-Q reference.
- Circular TE11 has two polarizations and TM01 also propagates in the
  illustrative 90 mm guide at 3 GHz. Cross-mode voltage ratios alone are not
  powers; check impedance normalization, direct E/H flux, reciprocity, and
  resolution before claiming a complete multimode balance.

## Learnings (coupled 3D FDTD + kinetics, `coupled3d.py` / `plasma3d.py`)

- Never use Maxwellian rate tables at low ionization degree: at 2-6 eV mean
  energy, Maxwellian Ar/He ionization rates exceed two-term Boltzmann rates by
  10^2 to >10^100 and produced a spurious RF-only filament breakdown. Blend
  `k = (1-w) kB + w kM` (8-row `GasModel.tabs`) with the local weight
  `w = nu_ee/(nu_ee + nu_eps)` (NRL nu_ee; `nu_eps = N(kexc/eps + kel +
  kiz Ip/eps)` from the Boltzmann rows). Do not use a fixed transition
  ionization degree: at 52 Td the criterion gives 1.7e-4 (Ar), 5.7e-4 (He),
  1.9e-4 (N2), and N2 rises to 1e-3-5e-3 at 2-6 eV, so one constant is wrong
  somewhere. Bound it with `--ee-scale 0.1` and `10`.
- `kinetic_cell` must deplete neutrals inside its substeps (nuclei conserved)
  and limit `dt` so the ionization cost cannot overdraw W. Never re-floor the
  electron energy to the gas temperature outside the ledger: without
  superelastic collisions a floor at Tg ~ 1e5 K silently drained ~8x the gas
  energy. Electrons may cool below Tg; the implicit elastic term reheats them.
- In the Lagrangian filament hydro, move boundaries with the time-centred
  velocity and book work as `p * 2 pi r_old ubar dt`; the geometric `p dA`
  form leaked ~2 % under strong shocks. Every run summary carries
  `ledger_filament_J.closure_residual_rel`; treat > 1e-6 (with real
  absorption) as a bug, not physics.
- Normalize energy residuals by the largest term including stored energy and
  a roundoff floor; dividing by a near-zero absorbed energy gives meaningless
  residuals. Numba raises ZeroDivisionError on scalar division by zero, and
  Gaussian tails give denormal 0/0: guard denominators and zero tiny intensities.
- In the taper, 10 MV/m at the 90 mm circular-guide focus requires ~15.7 MV/m
  in the WR-284 feed; if gas fills the feed, RF breakdown starts there first.
  Inspect max |E| over the gas region of the prefill, not only the focus.
- At 1 MPa, 300 K helium, the published secondary atomic-ion association
  coefficient gives a ~0.123 ns He+ -> molecular He2+ timescale; an atomic-only
  network cannot predict a dense-gas afterglow by assuming "only a few ns."
  Keep the prescribed-sheet EM animation separate from restricted CR examples
  until molecular recombination and pressure-broadened optical data are sourced.
- In a PEC-shorted cavity, input reflection is already unity in the lossless
  steady state. Demonstrate a plasma reflection plane through downstream field
  exclusion and a separate matched-guide slab reference; do not relabel an
  instantaneous cavity backward/forward ratio as a steady reflectivity.
- In `plasma3d.py`/`coupled3d.py`, `x2p`/`nx2p`/`cx2p` is the molecular ion
  He2+ (He+ + 2He conversion), NOT doubly ionized He2+. ~0.9 at 1 MPa is
  correct. Never label it "multiple ionization". The base kernels have no
  He(2+), superelastic or radiation-transport terms. `plasma3d_cr.py` adds
  them (lumped n=2 CR, He(2+) as `alpha`, Holstein-trapped resonance) plus
  bulk transport, gas heating and paraxial laser propagation, switched on by
  `coupled3d` flags. Run the `nocr`/`notransport`/`noheat`/`norefr`
  ablations to attribute effects; at He 1 MPa only CR changes |Gamma|.
- The base `_fil_diffuse` is not energy-conservative when Te varies radially:
  it rebuilds ue = 1.5 Te_local × new count. It also ignores the two electrons
  per He(2+). This leaked −1.3e-3 of filament energy in dense runs. The
  extended path uses `_fil_diffuse_cr`, which takes upwind face fluxes of
  count and mean energy. Always reconcile a ledger mismatch with an isolated
  per-kernel harness before trusting extreme temperatures.
- WarpX PIC/MCC swarm checks on Windows:
  - Install with a short `MAMBA_ROOT_PREFIX` (`D:\mm`), `--ssl-no-revoke`, and
    `matplotlib-base` (qt6 breaks long paths).
  - Make a uniform field with 1D Dirichlet plates: `warpx_potential_hi_z = -E0*L`,
    optionally `*sin(2*pi*f*t)`.
  - `UniformDistribution` needs `PseudoRandomLayout(n_macroparticles_per_cell=...)`.
  - Particle `ux` is γv.
  - The PIC EEDF tail in He is far below a same-mean-energy Maxwellian, so
    Maxwellian filament tables overestimate ionization.
  - Never validate RF drift from single snapshots: with N ~ 2000 macros, the
    snapshot SE v_th/sqrt(N) is ~30 % of the drift (seed-to-seed NRMSE ~0.44).
    Add `picmi.ReducedDiagnostic` `ParticleMomentum`/`ParticleEnergy`/
    `ParticleNumber` every ~20 steps, and average in RF phase windows. Normalize
    energy budgets by gross field work, never by the near-zero net (work - loss).
- At 1 MPa He, 532 nm needs 11 photons: 1–10 mJ alone gives a ~1e-12 J
  absorbed, micrometre seed. Even at ~5e24 m^-3, a 5–13 µm channel (~MΩ/m)
  cannot load a 3 GHz mode. The reflector is the RF-driven volume plasma, and
  it is lossy (nu >> omega: ~25 % reflected, ~70 % absorbed in the matched
  guide). Always run an RF-only control and a sub-threshold field (6 vs
  8 MV/m) before claiming that the laser triggers reflection.
- Numerical floors must never act as particle sources. `if ne < 1e-6: ne =
  1e-6; nxp = max(nxp, 1e-6)` in `kinetic_cell(_cr)` silently created
  unbooked electrons in every cell, and they avalanched into a "laser-
  triggered" breakdown that was really RF volume breakdown. Likewise any
  uniform fluid seed (even 1 m^-3) breaks down the whole volume. Keep
  `ne_true` for sources (`r_iz = kiz N ne_true`), floor only eps = w/ne, run
  with `seed_m3 = 0`, and seed locally (`laser_seed_m3`, k-order focal
  profile; pure-He 532 nm PPT and 10 ppm N2 give < 1 electron).
- With exactly-zero bulk/filament states, round-off rings (counts ~1e-40 /m)
  give eps = ue/cnt up to 1e44 eV. The free-electron (Allis–Rose) D then
  reaches 1e42 m^2/s, causing a zero pivot in `_thomas` (ZeroDivisionError).
  Capping eps in the energy flux is wrong: it leaked 0.3–1.3 % of filament
  energy at 200 kV/cm, where real Te exceeds the table ceiling. Instead,
  zero rings below `FIL_CNT_MIN` (1e-9 charges/m), cap Te only in the
  transport coefficients, and limit each face flux to the donor ring's count.
  Debug Numba kernels by swapping in `.py_func` with
  `np.seterr(all="raise", under="ignore")`, then read the frame locals.
- `scan_queue.py` now writes `.done` only after an `ok` RESULT line, and it
  hashes names longer than 160 characters (Windows MAX_PATH). Never set
  `NUMBA_NUM_THREADS` at module import (it breaks later imports in the same
  test process). A new `dx` needs `coupled3d prefill --dx` first.
- `scan_runner.py --tag` is appended to the auto run name (p, E, laser…), never
  substituted for it. A bare tag once made runs differing only in E share one
  output directory, and the later run silently overwrote the earlier one. Never
  put the same job in two concurrent plan files: the queue only skips jobs
  whose `.out` already has an `ok` RESULT, not jobs that are still running.
- Sparse kinetics masks require material-weighted edge-work deposition matching
  the Drude coefficient populations. Equal endpoint splitting sends boundary
  work into skipped empty nodes. Check edge, nodal and accepted kinetic work
  independently and reject current-memory work without supported electrons.
- CR cell kernels receive ground-state neutral density, not total nuclei.
  Handoffs must preserve heavy-gas internal energy and recover temperature with
  the actual heavy-particle count (molecules count once). Reject excess nuclei
  and incomplete kinetic timesteps; separate photon escape from total energy
  disposal. A tiny perturbation of ambient gas energy needs an explicit
  floating-point resolution bound, not a misleading relative-closure claim.
- `excited_avalanche` (radial reaction-diffusion with excited reservoirs): never put
  diagnostic counters or one-way global sources in the BDF state or its sparsity. Their
  dense rows force one finite-difference column group per unknown (10x slower). Exactly
  zero columns (e.g. molecular ions where n_e underflows) make scipy's adaptive num_jac
  overflow its step factor ("Factor is exactly singular"). Supply a fixed-step
  `approx_derivative(..., sparsity=(S, group_columns(S)))` Jacobian. Integrate counters
  afterwards with the log-mean rule, which is exact for exponential growth.
- He/Ar ordering lessons: in pure-gas dilute avalanches Ar/He stays at about 2 for every
  excited-state channel. Dense-regime (e-e) rates reverse it only at x >= 1e-4 and
  <= 10 Td. Dense-seed thresholds deposit 6-8 eV/atom (not self-consistent), so always
  report an energy ledger. A photon-feedback comparison must include He's own prompt
  He2(A) Hopfield emission (Y ~ 0.1); without it the result is biased toward Ar.
- Keep estimated laser seeding distinct from a modeled dense plasma: absence of a visible
  spark alone neither proves nor excludes an x ~ 1e-4 seed. Label dense-seed scans as
  hypothetical, compare them with the user's estimated density and its uncertainty, and map
  density to the seed volume before comparing total electrons.
- Match observables before claiming validation against the Technion S-band switch: the
  measured laser-to-light-onset delay is not the model's cavity field-collapse proxy, and the Stark
  density is not a peak-axis density. The published 16.6/22.3 kV/cm fields were calculated
  from power with peak/RMS convention unclear; an illustrative multiplier fitted to an
  assumed probability is not an independent field calibration.
