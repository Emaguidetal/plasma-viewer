# plasma-viewer
Topic-based research site for laser-induced plasma in gases and RF-plasma switching.
The homepage presents eight illustrated topics with links to model explanations,
interactive viewers, reports, data, and source.

- [Open the research homepage](https://emaguidetal.github.io/plasma-viewer/)
- [Read the model explanations](https://emaguidetal.github.io/plasma-viewer/topics.html)
- [Open the original coupled 3D viewer (1 MPa helium)](https://emaguidetal.github.io/plasma-viewer/cavity.html)
- [Open the parameter-scan explorer](https://emaguidetal.github.io/plasma-viewer/scans.html)
- [Explore RF, prescribed optics, real local PIC, gas expansion and emission dynamics](https://emaguidetal.github.io/plasma-viewer/dynamics/)
- [3D view: laser perpendicular vs parallel to the RF field](https://emaguidetal.github.io/plasma-viewer/dynamics/orientation3d.html)
- [Plasma switch in a WR-284 cavity: 3D rod growth along E, seed source in sealed helium, rod under continued RF drive](https://emaguidetal.github.io/plasma-viewer/switch/)
- [3D laser-induced plasma dynamics: seed → rods along E → bridge → glow → arc under continued drive](https://emaguidetal.github.io/plasma-viewer/switch/laser_plasma_3d.html)
- [Browse the research reports, source code, input data, and results](https://emaguidetal.github.io/plasma-viewer/research/)
- [Download the research package](https://emaguidetal.github.io/plasma-viewer/research-download.zip)

The dynamics page separates a 3D RF/CR-fluid calculation, a local assumed-RF-off
gas afterglow and real imposed-field WarpX/MCC diagnostics. Three read-only AI
reviewers accepted this restricted exploratory scope. It is not full-device or
experimental validation: PIC drift/work checks and gas velocity-peak convergence
remain incomplete, and physical optical attenuation/emission remain conditional.

## Site organization

- `index.html`: lightweight, responsive topic homepage; no embedded simulation data.
- `topics.html`: readable explanations for optics/geometry, gas/emission,
  local PIC validation, and parameter scans.
- `site.css`: shared homepage and explanation-page styling.
- `cavity.html`: the original root viewer, preserved without changes.
- `switch/`: current 1.48-MPa switch study; section anchors link directly to
  rod growth, seed sources, and continued drive.
- `research/`: reports, data, source, and download-scope notes.

The new homepage changes navigation, not the saved scientific results or viewer
calculations. Generated viewer files remain separate from the hand-authored portal.
