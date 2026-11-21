# Laser-induced plasma research package

This public package collects the literature reviews, research reports, PDFs, simulation source, input tables, tests, and selected scan outputs from the Plasma research workspace.

Start at the [research landing page](index.html) or open the [interactive viewer](../index.html) and [parameter-scan explorer](../scans.html).
The [new dynamics viewer](../dynamics/) separates RF/fluid, prescribed optics,
local assumed-quench gas dynamics and real imposed-field PIC.
The ZIP also contains a [standalone offline dynamics page](data/dynamics/dynamics.html).

## Reproducible project contents

- `source/src/rf_plasma/`: Python source, including the coupled 3D solver and the new breakdown-probability model.
- `source/cases/`: gas-kinetics, cross-section, and swarm input tables.
- `source/tests/`: project tests.
- `data/breakdown_probability/`: model output JSON and the probability/delay/focus plot.
- `data/coupled3d/`: selected calibration and scan tables, together with the interactive scan figures.
- `data/stochastic_seed/`: weak-seed jitter results.
- `data/dynamics/`: selected current cavity/afterglow outputs, dense RF probes,
  standalone viewer, shock and refinement diagnostics.
- `data/pic_dynamics/`: actual particle/EEDF histories, immutable run configs,
  validation summary and figures; passed and failed checks are retained.
- `pic_source/`: bounded WarpX/MCC drivers and execution guide.
- `agents/`: plasma, RF and optical AI reviewer profiles. These are
  instructions, not human credentials.
- `literature/` and `documents/`: research reports, literature reviews, and PDFs.

The original LXCat cross-section downloads are not redistributed here: their
headers state "All rights reserved." The derived swarm tables used by the
published model are included. See the reports for source attribution and
instructions for obtaining source data.
Original and converted LXCat cross-section tables are not included in the
PIC package; acquire licensed inputs separately and use the supplied export
driver. Actual computed particle histories and EEDF diagnostics are included.

## Large run archives not included

The workspace contains 3.23 GB of intermediate `.npz` field/frame archives, plus local runtime and solver files. GitHub Pages is not suitable for these large, numerous intermediate artifacts. They are not required to inspect the reports or reproduce the models from source and inputs; the published summary tables and figures cover the reported scans. This package does not claim to contain every saved transient field frame.

The breakdown-probability model and report distinguish measured observations from model assumptions and explicitly document the low-field ionization discrepancy. The figures should not be read as experimental validation.

The current dynamics study has 88 passing RF-plasma regressions and qualified
PIC mean-energy agreement. PIC drift/coarse-work checks fail, velocity peaks
are mesh-dependent, and the 50-kV/cm gas signal is below the conservative
floating-point resolution bound. Three AI reviews accepted only this restricted
exploratory scope. Radiation is conditional, not calibrated spectroscopy; optics
is prescribed, not validated propagation.

To reproduce the RF stage, restore `data/coupled3d/calibration.json` and
`prefill_cavity_dx1mm.npz` into `source/results/coupled3d`. See the source README
and PIC execution guide for exact commands and data-acquisition requirements.
