# RF and laser-plasma simulation benchmarks

This project contains reproducible numerical benchmarks and a **coupled,
cycle-resolved 3D RF-plasma simulation** (`coupled3d`). The coupled model runs
3D FDTD with bulk plasma kinetics at every gas node, plus a sub-cell laser
filament driven by PPT photoionization. It has been run for He, N2 and Ar in
the iris cavity and the taper. It is numerical evidence under stated
assumptions (PEC walls, illustrative inputs, a mean-energy EEDF closure whose
Maxwellization point is derived from cited rates rather than fitted), not a
validated prediction. Copper is represented only by separate first-order
analytical loss estimates.

The [benchmark report](../../literature/rf-plasma-simulation-benchmarks.md)
covers the earlier benchmarks. The
[coupled 3D report](../../literature/coupled-3d-rf-laser-plasma-simulation.md)
covers the gas-specific matrix, validation, sensitivity and limitations.

## RF, optics, gas expansion, emission and real local PIC

The [dynamics viewer](results/dynamics/dynamics.html) separates a bounded
perpendicular-seed Maxwell/CR-fluid run, an assumed-RF-off local cylindrical
Euler/CR afterglow, and actual 1D3V WarpX/MCC particles. The optical panel is a
prescribed Gaussian, not a validated propagated laser field. See the
[current report](../../literature/standing-wave-helium-plasma-reflector.md)
and [refinement metrics](results/dynamics/validation.json).

- [cavity_dynamics](src/rf_plasma/cavity_dynamics.py) maps signed edge work to
  occupied material, carries exact heavy-gas internal energy, and exports the
  RF phase clock, dense focus probe and conservative handoff.
- [radial_afterglow](src/rf_plasma/radial_afterglow.py) conserves nuclei and total
  energy, including boundary flux and CR sinks. It resolves local gas motion,
  not the continuing 5-us driven cavity or a sub-cell optical filament.
- [radiation_diagnostics](src/rf_plasma/radiation_diagnostics.py) separates
  modeled CR radiative energy from unresolved recombination disposal.
  Visible free-free is conditional; it is not a calibrated camera signal.
- [PIC execution guide](../pic/DYNAMICS_README.md) documents actual particle
  histories and passed/failed validation. Mean-energy agreement does not
  validate drift, exact MCC energy balance, ionization at 50 kV/cm, or a
  whole-cavity kinetic model.

The Sod shock benchmark tests the Euler implementation; the helium example
shows subsonic expansion, not an established experimental shock. Velocity
peaks remain mesh-dependent. Low-field gas heating is below the conservative
floating-point resolution bound against ambient gas energy and is explicitly
flagged, rather than presented as resolved motion.

From this directory, use existing environments, new output paths, and check
each exit status:

```powershell
$env:PYTHONPATH = 'src'
$env:NUMBA_NUM_THREADS = '2'
& .\.venv\Scripts\python.exe -m unittest discover -s tests
& .\.venv\Scripts\python.exe -m rf_plasma.cavity_dynamics --output results\rerun-cavity --field-kv-cm 200 --block 1 --threads 2
& .\.venv\Scripts\python.exe -m rf_plasma.radial_afterglow --source results\rerun-cavity\handoff.npz --output results\rerun-afterglow --nr 160
```

The RF run needs the saved calibration and `prefill_cavity_dx1mm.npz` in
`results\coupled3d`. These are included in the public package's
`data\coupled3d`; restore them to that result directory before executing.

`--laser-axis x` (default) puts the laser along x, perpendicular to the TE10
field Ey, as in the experiment. `--laser-axis y` puts it along Ey (parallel,
through the broad wall). The matched comparison and 3D viewer are built from
`results\dynamics\orientation\E{50,200}_{perp,par}`:

```powershell
& .\.venv\Scripts\python.exe -m rf_plasma.cavity_dynamics --output results\dynamics\orientation\E200_par --field-kv-cm 200 --block 1 --threads 3 --laser-axis y
& .\.venv\Scripts\python.exe -m rf_plasma.orientation_visuals
```

`coupled3d` laser-trigger runs and `streamer2d` use the parallel (laser along
y) geometry.
Alternatively generate a new calibration/prefill with `coupled3d calibrate`
and `coupled3d prefill --geometry cavity --dx 0.001 --threads 2`, keeping the
new calibration separate from this published snapshot. The viewer builder
uses the explicitly named executed cases in `results\dynamics`, not an
arbitrary rerun directory.

## Dynamic standing-wave gas-sheet viewer

Open the [interactive 2D/3D animation](results/standing_wave_helium_1MPa/standing_wave.html)
in a browser, or view the [GIF](results/standing_wave_helium_1MPa/standing_wave.gif).
The [report](../../literature/standing-wave-helium-plasma-reflector.md) explains
the executed 3D WR-284 cavity comparison at 1 MPa absolute helium and 0.1 MV/m.
The transverse sheet is **ionized gas, not a metal insert**; its density history
is prescribed, not calculated from a laser pulse. Sequential ionization,
superelastic exchange and nonlocal reduced radiation transport are tested in a
separate restricted low-density atomic demonstration, **not coupled into a
claimed 1 MPa sheet lifetime**. Reproduction commands and missing dense-helium
data are documented in the report.

## What is here

- [cases/reference.json](cases/reference.json): illustrative SI inputs, not measured apparatus values.
- [src/rf_plasma](src/rf_plasma): analytical waveguide/copper and Drude references,
  openEMS PEC geometries, Gaussian optics, a prescribed-density resonator/current
  demonstration, and an abstract single-ion kinetics benchmark.
- [tests](tests): physics and numerical regression tests.
- [results](results): saved inputs, metrics, field data, and figures. Use the
  explicitly linked results in the report; directories named `prototype` or
  older superseded runs are historical, not validated final designs.
- [.tools/openEMS](.tools/openEMS) and `.venv`: project-local Windows solver and
  Python environment. No system-wide install or WSL is required here.

The reference pressure is **100 psig**, or 790800.7 Pa absolute assuming
101325 Pa ambient. The 3 GHz frequency, 10 MV/m peak empty-device focus field,
532 nm and 1 mJ at focus, and geometry and waist are illustrative. The pulse
endpoints are 300 ps and 3 ns intensity FWHM. The focus is not a predicted
plasma density or an apparatus operating specification.

## Versions and provenance

The executed environment is Python 3.14.3, openEMS 0.37.0, CSXCAD 0.7.0,
NumPy 2.5.5, SciPy 1.18.1, Matplotlib 3.11.2, and h5py 3.16.0.
The [official Windows package](https://github.com/thliebig/openEMS-Project/releases/tag/v0.37.0)
was verified against its published SHA-256:
`E4F14FD6E9BE7D0C1FD9C244768C21C1C2B62B8AA81F08253EBF6E06040CDA5A`.
The included Windows wheels target 64-bit CPython 3.13 or 3.14; the broader
Python requirement in [pyproject.toml](pyproject.toml) applies to the pure
Python code, not necessarily to the solver wheels.

If rebuilding an environment, use the matching official solver wheels and
the dependency ranges in `pyproject.toml`. The scripts set
`CSXCAD_INSTALL_PATH` in the current process when using the bundled solver.
Do not install or set it globally. The openEMS maintainers describe the
[Windows package](https://docs.openems.de/en/latest/install/package.html) as
a quick start and note possible performance limitations.

## Run on this Windows workspace

In PowerShell, start from this directory. Existing output directories are never
overwritten; each run below uses a new unique directory.

```powershell
Set-Location -LiteralPath 'D:\CST projects\Plasma research\simulations\rf_plasma'
$python = '.\.venv\Scripts\python.exe'
$run = Join-Path (Get-Location) ("results\reruns\run-" + [guid]::NewGuid())

# Pure Python tests and analytical inputs
& $python -B -m unittest discover -s tests -v
& $python -B -m rf_plasma.cli '.\cases\reference.json'

# Optical focus scan and prescribed-density transient, not ionization models
& $python -B -m rf_plasma.laser_study --output (Join-Path $run 'optics')
& $python -B -m rf_plasma.transient --pulse-ps 300 --output (Join-Path $run 'transient-300ps')
& $python -B -m rf_plasma.transient --pulse-ps 3000 --output (Join-Path $run 'transient-3000ps')

# PEC rectangular guide, iris cavity, closed cavity, and prescribed Drude slab
& $python -B -m rf_plasma.cold_rf guide --mesh-mm 2 --threads 2 --output (Join-Path $run 'guide')
& $python -B -m rf_plasma.cold_rf cavity --mesh-mm 2 --threads 2 --calibration (Join-Path $run 'guide\metrics.json') --output (Join-Path $run 'cavity')
& $python -B -m rf_plasma.closed_cavity --mesh-mm 2 --threads 2 --output (Join-Path $run 'closed-cavity')
& $python -B -m rf_plasma.cold_rf slab --mesh-mm 1 --threads 2 --output (Join-Path $run 'drude-slab')

# Exploratory forward and reverse taper excitation, with E/H flux planes
& $python -B -m rf_plasma.transition_rf taper --mesh-mm 2 --slices 20 --threads 2 --mode TE11-90 --flux-validation --output (Join-Path $run 'taper-forward')
& $python -B -m rf_plasma.transition_rf taper --mesh-mm 2 --slices 20 --threads 2 --mode RECT-TE10 --flux-validation --output (Join-Path $run 'taper-reverse')

# Checks the specific saved benchmark snapshot, not the new reruns above
& $python -B -m rf_plasma.validate_results --output (Join-Path $run 'saved-result-checks')
```

Check each command's exit status before using its output. The field-calibration
ratio is approximate, not an independently validated power requirement. A new
`validate_results` input tree must contain the same named benchmark directories
and metadata as the saved snapshot. It is not a validator for arbitrary cases.
It currently exits nonzero because the finest saved taper run still exceeds the
1% rectangular flux/modal discrepancy target; inspect its JSON and plots for
the individual passed and failed checks.

Runs use two threads by default and refuse a preliminary estimate larger than
the lesser of 8 GiB or one quarter of currently available RAM. That estimate
is **not** an enforced memory limit. Run large jobs sequentially on shared
computers, monitor actual resources, and do not assume the solver obeys the
estimate.

## Boundaries of interpretation

- The prescribed Drude slab has a fixed density and collision rate. It does not
  calculate laser seeding, avalanche, chemistry, or plasma evolution.
- The resonator demonstration prescribes density versus time, assumes coupling
  Qs, and accounts for collision and birth-mixing energy. It is not a reduction
  validated against either 3D device; the initial resonator is prefilled.
- The circular guide supports two TE11 polarizations and TM01 at 3 GHz. Taper
  scattering, mode-power normalization, reciprocity, and total flux require
  further checks; a modal deficit is not automatically plasma absorption.
- Copper skin depth, ideal cavity Q, and straight-guide attenuation are
  separately calculated approximations. They are not copper-loss results for
  the 3D iris or transition models.
- The openEMS benchmarks above contain no reaction network or laser
  ionization. Those are only in `coupled3d`, whose gas data carry per-entry
  verification flags (`cases/gas_kinetics_model.json`). Some entries are
  unverified or are assumptions; see the coupled report. No experimental
  comparison exists.

## Coupled 3D RF-plasma simulation (`coupled3d`)

```powershell
$py = '.\.venv\Scripts\python.exe'; $env:PYTHONPATH = 'src'
& $py -m unittest tests.test_plasma3d tests.test_gas_model
& $py -m rf_plasma.coupled3d calibrate
& $py -m rf_plasma.coupled3d prefill --geometry cavity
& $py -m rf_plasma.coupled3d prefill --geometry taper
& $py -m rf_plasma.coupled3d run --geometry cavity --gas Ar --energy 0.03 --threads 8
& $py -m rf_plasma.coupled3d matrix --geometry cavity taper --skip-existing --threads 8
& $py -m rf_plasma.plots3d
```

`run` options include:

- `--no-laser` and `--no-rf`;
- `--fwhm`, `--energy` and `--waist`;
- `--seed`;
- `--ee-scale`, a multiplier on ν_ee in the EEDF closure criterion
  (default 1; used only for sensitivity runs);
- `--block`, the FDTD steps per coupling step;
- `--rings`, the filament radial resolution;
- `--tag`.

Each case writes the following to `results/coupled3d/<case>/`:

- `summary.json`, which holds:
  - the provenance;
  - the RF balance;
  - the kinetic ledgers, with `ledger_filament_J.closure_residual_rel`;
  - the peak values;
- `timeseries.npz`, `rf_series.npz` and `snapshots.npz`;
- three figures.

`plots3d` writes `matrix_overview.png`, `laser_thresholds.png`,
`sensitivity.png` and `matrix_summary.json`.

The high illustrative field, laser pulse, and pressurized gases require an
independent institutional laser/RF/electrical/pressure safety assessment; this
project is not a construction or operating procedure.
