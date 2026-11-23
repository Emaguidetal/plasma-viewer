# streamer2d — axisymmetric channel-growth side-model (He, RF field ∥ laser channel)

This is a standalone 2D (r, z) drift-diffusion-Poisson model. It follows a laser-made plasma channel in
high-pressure helium driven by a uniform RF field E_z(t) = E0 sin(2π·2.95 GHz·t) parallel to the channel
axis. It is **not coupled** to the 3D FDTD solver.

All numbers below are **simulation results** of a reduced fluid model. They have not been validated against
helium channel-growth measurements, because none exist in the workspace. The only measured channel geometry
and speed are for **air**.

- Code: `src/rf_plasma/streamer2d.py` (model, CLI) and `src/rf_plasma/streamer2d_report.py` (aggregation, figures).
- Tests: `tests/test_streamer2d.py`.
- Outputs: `results/streamer2d/<run>/{summary.json, history.npz, snap_t*.npz}` (one directory per run),
  `results/streamer2d/summary.json` (cross-case table and comparisons) and `results/streamer2d/figures/`.

## Equations

```
∂n_e/∂t + ∇·Γ_e = ν_iz(|E|/N) n_e − α n_e n_i,       Γ_e = −μ_e n_e E − D_e ∇n_e
∂n_i/∂t + ∇·Γ_i = ν_iz(|E|/N) n_e − α n_e n_i,       Γ_i = +μ_i n_i E − D_i ∇n_i
−∇·(ε0 ∇φ) = e (n_i − n_e),                            E = −∇φ + E0 sin(2πft) ẑ
```

- **Geometry:** axisymmetric finite volumes, with face areas 2πr·Δz and π(r₊² − r₋²) and cell volumes π(r₊² − r₋²)Δz.
  The axis has zero area, so no flux crosses it. Do not read this as a Cartesian 2D model.
- **Boundaries:** grounded φ = 0 at r = 3 mm and z = ±z_half. The applied field is superposed as a uniform E_z,
  so the far ends act like ideal electrodes for the space-charge potential. Particles can only leave the
  domain (outflow-only drift, zero diffusive boundary flux). The cumulative outflow is recorded, and
  charge conservation is checked including it.
- **Initial state:** n_e = n_i = n0·exp(−r²/a² − z²/L²) + n_bg, quasi-neutral. Defaults: n0 = 1e23 m⁻³,
  a = 10 µm, L = 1 mm, n_bg = 1e9 m⁻³. At t = 0 the n_e ≥ 1e20 contour already spans ±2.63 mm, with radius 26 µm.

## Closures, assumptions and sources

| Item | Value or treatment | Source / status |
|---|---|---|
| μ_e·N, k_iz, mean energy vs E/N | Log-log interpolation of `cases/swarm/He_boltzmann.csv` (0.03–3000 Td) at the local \|E\|/N (local-field approximation, LFA). k_iz table noise (±1e-35) below 3.6 Td is set to 0. | BOLOS two-term Boltzmann with LXCat He cross-sections (existing workspace table) |
| D_e | μ_e·T_e with T_e = (2/3)·mean energy (Einstein relation) | as specified |
| Ions | All ions are treated as He2⁺. K0 = 16.7 cm²/Vs scaled by N0/N, so μ_i = 9.3e-4 m²/Vs at 0.2 MPa. D_i = μ_i·kT_g/e. | **Assumption.** Ellis 1976 value as quoted in `cases/gas_kinetics_model.json`, not re-verified. At 0.2 MPa, He⁺ → He2⁺ conversion takes about 4 ns, so the early ions are partly He⁺. |
| Recombination | α = 1e-14 m³/s, constant | **Unverified.** Saifutdinov 2023 gives 5e-15·(T_e/T_g)⁻¹ (dissociative) plus 5e-39·N·(T_e/T_g)⁻¹ m⁶/s (e + He2⁺ + He). These span about 1e-16 to 2e-13 m³/s at 0.2 MPa, depending on T_e. **This choice controls the on-axis n_e** (see the sensitivity section). |
| Photoionization | none | Pure He has no efficient photoionization. Radial and tip growth come only from drift, diffusion and the background seed. |
| Background pre-ionization | 1e9 m⁻³ | A 1e12 m⁻³ case was not run. With ν̄_iz ≤ 8e8 s⁻¹ (0.2 MPa, 46 Td), a background grows by ≤ e⁸ ≈ 3e3 in 10 ns, to ≤ 3e15 m⁻³. That is far below the 1e20 contour, so it cannot change the results. |
| Field amplitude | E0 is treated as **peak** | Technion "calculated" fields of 16.6 and 22.3 kV/cm, peak vs RMS **not stated**. If they are RMS, E/N is √2 higher. |
| Frequency | 2.95 GHz (as specified) | The experiment used 2.767 GHz. Under LFA, only the cycle timing changes (6% shift). |
| Not modelled | gas heating, electron energy equation, excited-state, metastable and stepwise ionization, Penning ionization by impurities, radiation, magnetic fields, circuit and waveguide feedback | — |

**LFA validity.** At 0.2 MPa the electron energy relaxation rate is about 2.6e10 s⁻¹ and ω = 1.85e10 s⁻¹,
so ωτ_ε ≈ 0.7. The LFA overestimates how closely T_e follows E(t). The error falls with pressure, roughly
∝ 1/N. Inside the screened channel the local |E| is small, so the LFA gives a near-thermal T_e and the
large low-field μ_e (0.45 m²/Vs at 0.2 MPa). An energy-equation model would give a hotter, less mobile
electron population there.

### Swarm-derived scales (cycle average over |sin|)

| Case | E0/N (Td) | ν_iz at E0 (s⁻¹) | ν̄_iz (s⁻¹) | 1/α_iz at E0 (µm) | KPP speed 2√(D̄ν̄) (m/s) | ν̄/α (m⁻³) |
|---|---|---|---|---|---|---|
| 0.2 MPa, 16.6 kV/cm | 34.4 | 7.7e8 | 3.1e8 | 103 | 1.6e4 | 3.1e22 |
| 0.2 MPa, 22.3 kV/cm | 46.2 | 1.9e9 | 8.0e8 | 57 | 2.8e4 | 8.0e22 |
| 0.446 MPa | 46.2 | 4.3e9 | 1.8e9 | 25 | 2.8e4 | 1.8e23 |
| 1.0 MPa | 46.2 | 9.6e9 | 4.0e9 | 11 | 2.8e4 | 4.0e23 |
| 1.825 MPa | 46.2 | 1.75e10 | 7.3e9 | 6.2 | 2.8e4 | 7.3e23 |
| 1.0 MPa, 80 kV/cm | 33.1 | 3.4e9 | 1.35e9 | 22 | 1.5e4 | 1.35e23 |

## Numerics

- **Fluxes:** Scharfetter–Gummel fluxes on face fields. The density update is explicit, conservative
  finite-volume forward Euler. dt = 0.45 × the exact per-cell positivity bound for SG outflow (this covers
  drift and diffusion). Also ν·dt ≤ 0.2 and dt ≤ T_RF/40. Snapshot and diagnostic times are hit exactly.
  If a density goes negative, the step is retried with dt halved. **No density floors are used.**
- **Poisson solve:** a direct sparse LU (SuperLU, `splu`, MMD_AT_PLUS_A ordering). The vacuum Laplacian is
  factorised once per grid.
- **Semi-implicit dielectric relaxation:** τ_d = ε0/(e·n·μ) ≈ 1e-14 s at 1e23 m⁻³, far below the transport
  dt (about 1e-13 to 1e-12 s), so a plain explicit Poisson step would be unstable. The charge flux is
  linearised as J^{n+1} ≈ J^n + Λ_f (E^{n+1} − E^n). Λ_f = h·dt·max(σ_f, σ_env), where σ_f is the exact
  derivative of the SG charge flux with respect to the face field. σ_env is an envelope built from the
  table's maximum μ_e and max(n_L, n_R); h = 1.5. Substituting this into Gauss's law gives the SPD operator
  −∇·((ε0 + Λ)∇). That operator is factorised and reused for up to 100 steps, and refactorised earlier if
  dt·σ_f > 1.5(ε0 + Λ_f) anywhere.
- **Why the semi-implicit fix is justified:** the 0D amplification factor is 1 − x/(1 + y), with
  x = dt·σ/ε0 and y = Λ/ε0. It is stable for x ≤ 2(1 + y). Over-estimating Λ only adds damping of the
  relaxation mode, and steady or quasi-steady (screened) states are unchanged.
- **Check against exact mode:** compared with "exact" mode (Λ = dt·σ_f refactorised every step) on a coarse
  grid at 1 ns, the L2 difference is 0.1% in |E| and 0.5% in n_e (max pointwise 1.9% and 0.5%). Exact mode
  is 10× slower.
- **Grid (production):**
  - r: 1 µm cells to 20 µm, ×1.05 growth to 8 µm, 8 µm cells to 0.5 mm, then ×1.15 growth to 3 mm (132 cells).
  - z: uniform 5 µm cells over ±4 mm (1600 cells) or ±8 mm (3200 cells, the `_z8` runs).
  - Total: 211k or 422k cells.

## Verification (`python -m unittest tests.test_streamer2d`: 5 tests, OK, about 4–10 s)

| Test | Result |
|---|---|
| Poisson vs manufactured axisymmetric solution on a stretched grid | max error 1.87e-4 → 4.68e-5 when refined 2× (order 2.0) |
| Poisson vs analytic charged-cylinder E_r(r) | passes (tolerance in test) |
| Charge conservation, ionization on, space charge and outflow both nontrivial | \|Q_i − Q_e + outflow − Q_0\|/Q_e ≤ 3.5e-17 over 137 steps |
| E0 = 0, no ionization: Gaussian diffusion vs analytic | ⟨r²⟩ 0.9998, ⟨z²⟩ 1.0000, peak 0.998, particle count 1.0008 of analytic |
| Swarm interpolation and Einstein relation | passes |

In production runs the maximum relative charge residual is ≤ 9e-15.

**Sensitivities** (`summary.json → sensitivity`, `domain_size_check`):

| Check | Largest change in r(1e20) mid / maxz, r_σ, length, peak n_e | Tip speed |
|---|---|---|
| Grid 2× coarser (dz 10 µm, dr_min 2 µm), 0.2 MPa at 3 ns | ≤ 1.1% | −2.7% |
| Grid 2× coarser, 1 MPa at 1.5 ns | ≤ 2.8% | −11% |
| Grid 2× coarser, 1.825 MPa at 0.96 ns | ≤ 2.3% | +3.6% |
| Domain ±4 mm vs ±8 mm (four cases, common window) | ≤ 2.3% (length), ≤ 1.9% otherwise | — |
| **α = 1e-15 vs 1e-14**, 0.2 MPa / 22.3 kV/cm at 3 ns | peak n_e **+177%** (1.5e23 vs 5.4e22, still rising), r_maxz +24%, r_mid −3%, r_σ −15% | **+36%** |

Two caveats on these checks:
- At 1–1.8 MPa the coarse-grid tip speed differs by up to 11%. After 2 ns, the outer parts of the
  high-pressure tip bulbs (r > 0.5 mm) sit on the stretched radial grid (cells ≥ 8 µm, versus 1/α_iz of
  6–11 µm), so they are under-resolved.
- The σ-weighted radius shows narrow periodic dips near the zero-crossings of E_app(t). These come from
  evaluating the instantaneous LFA mobility in the diagnostic; the transport is not affected.

## Results

A `_z8` run uses the ±8 mm domain. All runs use 1 thread each, with at most 3 running at once.

| Case | Run | Stop reason, t_f | Wall time |
|---|---|---|---|
| 0.2 MPa, 16.6 kV/cm (34 Td) | He0p2MPa_E16p6 (±4 mm; tips stay within ±3.5 mm) | t_end, 10 ns | 3639 s |
| 0.2 MPa, 22.3 kV/cm (46 Td) | He0p2MPa_E22p3_z8 | t_end, 10 ns | 12601 s |
| 0.446 MPa, 46 Td (49.7 kV/cm) | He0p446MPa_46Td_z8 | z-boundary, 5.90 ns | 3972 s |
| 1.0 MPa, 46 Td (111.5 kV/cm) | He1MPa_46Td_z8 | z-boundary, 3.04 ns | 2782 s |
| 1.825 MPa, 46 Td (203.5 kV/cm) | He1p825MPa_46Td_z8 | z-boundary, 1.82 ns | 2733 s |
| 1.0 MPa, 80 kV/cm (33 Td) | He1MPa_E80 (±4 mm) | z-boundary, 4.86 ns | 333 s |
| 0.2 MPa, 22.3 kV/cm, n0 = 1e24, a = 20 µm | He0p2MPa_E22p3_n1e24_a20_z8 | t_end, 5 ns | 8474 s |

The ±4 mm runs of the high-pressure cases (0.15–1.9 ks each) and the sensitivity runs are kept as well.
Total wall time over all runs is 39 ks (about 10.8 h), spread over at most 3 single-threaded processes.

### Channel radius, tip speed, n_e

Notes on the table:
- r_mid: radius of the n_e ≥ 1e20 m⁻³ contour at z = 0.
- r_maxz: maximum of that radius over z, which sits in the tip "bulbs".
- r_σ: conductivity-weighted RMS radius at z = 0.
- "—": the run ended before that time.
- v_tip is a linear fit of both 1e20 tips for t ≥ 1 ns. "late" is the fit over the last ns, showing acceleration.

| Case | r_mid 2 / 5 / 10 ns (µm) | r_maxz 2 / 5 / 10 ns (µm) | r_σ 2 / 5 / 10 ns (µm) | v_tip fit / late (m/s) | n_e(0,0) at t_f (m⁻³) | peak \|E\| (V/m) |
|---|---|---|---|---|---|---|
| 0.2 MPa, 16.6 | 42 / 56 / 75 | 54 / 93 / 172 | 17 / 24 / 33 | 9.9e4 / 1.4e5 | 2.35e22 | 8.3e6 |
| 0.2 MPa, 22.3 | 46 / 67 / 105 | 99 / 241 / 645 | 18 / 28 / 48 | 4.5e5 / 7.0e5 | 7.1e22 | 1.6e7 |
| 0.446 MPa | 40 / 86 / — (114 at 5.9 ns) | 135 / 515 / — (762) | 16 / 31 / — | 9.6e5 / 1.6e6 | 2.4e23 | 4.9e7 |
| 1.0 MPa | 55 / — / — (88 at 3.0 ns) | 266 / — / — (811) | 19 / — / — | 2.2e6 / 2.7e6 | 5.5e23 | 1.3e8 |
| 1.825 MPa | — (43 at 1 ns; 69 at 1.8 ns) | — (163 at 1 ns; 765) | — (16 at 1 ns; 29) | 4.3e6 (t ≥ 1 ns) | 1.3e24 | 2.7e8 |
| 1.0 MPa, 80 kV/cm | 34 / — / — (54 at 4.9 ns) | 68 / — / — (162) | 14 / — / — | 2.1e5 / 4.4e5 | 1.4e23 | 3.0e7 |
| 0.2 MPa, 22.3, n0 = 1e24, a = 20 µm | 71 / 92 / — | 286 / 490 / — | 31 / 42 / — | 5.5e5 / 6.1e5 | 6.5e22 | 1.9e7 |

Snapshot values of n_e(r = 0, z = 0) are in `summary.json → cases.*.ne_axis_z0_m3`. That list also gives the
average over r ≤ 100 µm, which is closer to a line-of-sight measurement.

### What the model shows

1. **The channel is screened, and growth starts at its ends.** The conducting channel polarises in the
   axial RF field. The field inside it collapses, and the field at its ends is enhanced by 4–13× over E0
   (peak |E| above). Ionization fronts start at the two tips, propagate along E in both directions, and
   swell into "bulbs". By 10 ns at 22.3 kV/cm the channel is dumbbell-shaped (see `figures/ne_map_*`): the
   mid-plane radius is 105 µm, while the 1e20 contour near the tips reaches 0.65 mm.
2. **Mid-plane radial growth is slow and diffusion/avalanche-limited.** It runs at 4–8e3 m/s at 0.2 MPa,
   below the unscreened Fisher–KPP estimate of 1.6–2.8e4 m/s, because the field is screened. The
   conductivity-weighted radius stays at 33–48 µm (diameter < 0.1 mm) at 10 ns.
3. **Higher pressure at the same E/N grows faster.** With n0 fixed, the seed is less screened relative to
   the gas (n0/N² is smaller), and ν̄/α rises ∝ N. The result is faster, accelerating tips and higher n_e.
   At 1–1.8 MPa the tips run at 2–4e6 m/s and reach the ±8 mm domain end within 2–3 ns.
4. **On-axis n_e tends toward an ionization–recombination balance near ν̄_iz/α.** At 0.2 MPa the values are
   2.35e22 at 16.6 kV/cm (ν̄/α = 3.1e22) and 5.4–7.1e22 at 22.3 kV/cm (ν̄/α = 8.0e22). The dense
   1e24 seed decays to the same 6–6.5e22 within about 2 ns. **This level scales as 1/α**, and α is unverified.

### Comparison with measurements (He data where available; air data are only an order-of-magnitude reference)

| Quantity | Measured | Model | Comment |
|---|---|---|---|
| He n_e, 0.2 MPa | ≈ 6e22 m⁻³ peak (Stark, line-of-sight through a 5 mm hole, about 300 kW output). About 3× lower at about 90 kW. | On axis at z = 0: 7.1e22 (22.3 kV/cm, 10 ns), 2.35e22 (16.6 kV/cm). Average over r ≤ 100 µm: 2.6e22 and 3.2e21. | Same order. The agreement is **set by α = 1e-14** (n_e ≈ ν̄/α); with α = 1e-15, n_e exceeds 1.5e23 by 3 ns. Not a validation. |
| Channel diameter | 0.2–0.6 mm (**air**, Technion imaging) | 1e20 diameter at z = 0: 0.15 / 0.21 mm at 10 ns (16.6 / 22.3). Tip bulbs: 0.34 / 1.3 mm. σ-weighted: 0.07 / 0.10 mm. | Comparable scale, but a different gas. The model's "diameter" depends strongly on its definition. |
| Growth along E | ≈ 5e7 cm/s = 5e5 m/s (**air**) | 0.2 MPa: 1.0–1.4e5 m/s (16.6) and 4.5–7.0e5 m/s (22.3). Higher pressures at 46 Td: 1–4e6 m/s. | 22.3 kV/cm at 0.2 MPa happens to match the air number. He and air differ (no photoionization, other transport). |
| Delay, low minus high output-power groups (He, 0.2 MPa) | +30 ± 15 ns (laser entry → light-emission onset; ±15 ns is the sum of jitters) | From the model's already-formed channel: Δ ≈ 7 ns for the midplane diameter proxy (including a linear extrapolation beyond 10 ns), 3.4 ns for the any-z diameter proxy, and > 3.2 ns for line-conductance doubling. | These are different observables in a different setup; their numerical differences are not a validation or a measured seed constraint. The model starts from a prescribed 1e23 m⁻³ channel and runs only 10 ns. |

## Limitations

- **Single fluid under LFA, no energy equation.** T_e in the screened core and the high-field tips is not
  self-consistent, and ωτ_ε ≈ 0.7 at 0.2 MPa.
- **Recombination coefficient is unverified.** It controls the on-axis n_e (about ν̄/α) and changes the tip
  speed by about 36% per decade of α. The ion species (He2⁺ only) and its mobility are also assumptions.
- **No gas heating or hydrodynamics.** At n_e ≈ 1e23–1e24 and |E| up to 2.7e8 V/m (1.8 MPa), Joule heating
  over nanoseconds is likely large, so the high-pressure late-time results lie outside the model's range.
- **Domain.** Grounded ends at ±4 mm or ±8 mm act like electrodes. Runs stop when the 1e20 contour comes
  within 0.25 mm of a boundary, so the high-pressure cases end at 1.8–5.9 ns and have no 5 or 10 ns values.
  The 1 MPa / 80 kV/cm case was run only on ±4 mm (to 4.86 ns). There is no waveguide wall, mode structure,
  or circuit feedback: E0 is fixed and uniform, and power absorbed by the plasma is not taken from the RF.
- **Resolution at high pressure.** dz = 5 µm is close to 1/α_iz (6 µm at 1.825 MPa). Tip speeds there have
  about 10% grid uncertainty, and the outer tip bulbs are on stretched radial cells.
- **Seed is prescribed.** The laser channel is a fixed Gaussian (n0, a, L); it does not come from a laser
  model. Absolute "delays" are therefore not predicted (see the delay row above).
- **Air comparison values are not helium data.**

## Reproduce

```powershell
cd simulations\rf_plasma
$env:PYTHONPATH = 'src'
.venv\Scripts\python.exe -m unittest tests.test_streamer2d
.venv\Scripts\python.exe -m rf_plasma.streamer2d list
.venv\Scripts\python.exe -m rf_plasma.streamer2d run He0p2MPa_E16p6 He1MPa_E80
.venv\Scripts\python.exe -m rf_plasma.streamer2d run He0p2MPa_E22p3 He0p446MPa_46Td He1MPa_46Td He1p825MPa_46Td --z-half 8e-3 --suffix _z8
.venv\Scripts\python.exe -m rf_plasma.streamer2d run He0p2MPa_E22p3_n1e24_a20 --z-half 8e-3 --t-end 5e-9 --suffix _z8
# +/-4 mm reference runs of the high-pressure cases (used for the domain check)
.venv\Scripts\python.exe -m rf_plasma.streamer2d run He0p2MPa_E22p3 He0p446MPa_46Td He1MPa_46Td He1p825MPa_46Td
# sensitivities
.venv\Scripts\python.exe -m rf_plasma.streamer2d run He0p2MPa_E22p3 --alpha 1e-15 --t-end 3e-9 --suffix _alpha1e-15
.venv\Scripts\python.exe -m rf_plasma.streamer2d run He0p2MPa_E22p3 --dz 10e-6 --dr-min 2e-6 --t-end 3e-9 --suffix _coarse
.venv\Scripts\python.exe -m rf_plasma.streamer2d run He1MPa_46Td --dz 10e-6 --dr-min 2e-6 --t-end 1.5e-9 --suffix _coarse
.venv\Scripts\python.exe -m rf_plasma.streamer2d run He1p825MPa_46Td --dz 10e-6 --dr-min 2e-6 --t-end 0.96e-9 --suffix _coarse
# aggregate table, summary.json and figures
.venv\Scripts\python.exe -m rf_plasma.streamer2d_report
```

Environment: Python 3.14.3, numpy 2.5.5, scipy 1.18.1, numba 0.68.1, matplotlib 3.11.2, Windows 11. When run
as the CLI, the BLAS/OpenMP thread counts are set to 1 per process. Importing the module does not change the
thread environment.

## Figures (`figures/`)

- `ne_map_<case>.png`: n_e(r, z) at 4–5 times, log scale, with the 1e20 contour in cyan and the r axis
  mirrored.
- `radius_length_vs_time.png`: r_1e20 (at z = 0 and max over z), r_σ, axial extent, and n_e(0, 0)(t), with
  the measured He 6e22 and the air radius band marked.
- `tip_velocity_vs_EN.png`: tip speed vs applied E0/N, with the air value 5e5 m/s marked.
- `diameter_vs_measured.png`: diameters at 1, 2, 5 and 10 ns for the three definitions, with the air band
  0.2–0.6 mm.
