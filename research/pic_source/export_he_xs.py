"""Export the LoKI/LXCat He cross sections (the same file BOLOS uses for the fluid tables) to
WarpX MCC two-column files (energy eV, sigma m^2) plus a manifest with thresholds.

Run with the rf_plasma venv (needs bolos):
    ..\\rf_plasma\\.venv\\Scripts\\python.exe export_he_xs.py
"""
import json
from pathlib import Path

import numpy as np
from bolos import parser

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "rf_plasma" / "cases" / "cross_sections" / "LoKI_He_LXCat.txt"
OUT = HERE / "xs_he"


def main():
    OUT.mkdir(exist_ok=True)
    with SRC.open(encoding="utf-8") as fh:
        procs = [p for p in parser.parse(fh) if p["target"] == "He"]
    manifest = {"source": str(SRC.relative_to(HERE.parent.parent)), "processes": []}
    n_exc = 0
    for p in procs:
        d = np.asarray(p["data"], dtype=float)
        kind = p["kind"]
        if kind in ("ELASTIC", "EFFECTIVE"):
            name, energy = "elastic", None
        elif kind == "EXCITATION":
            n_exc += 1
            name, energy = f"excitation{n_exc}", float(p["threshold"])
        elif kind == "IONIZATION":
            name, energy = "ionization", float(p["threshold"])
        else:
            continue
        if energy is not None and d[0, 0] > energy:
            d = np.vstack([[energy, 0.0], d])
        if d[0, 0] > 0.0:
            d = np.vstack([[0.0, d[0, 1] if energy is None else 0.0], d])
        fname = f"{name}.dat"
        np.savetxt(OUT / fname, d, fmt="%.6e")
        manifest["processes"].append({"name": name, "kind": kind, "file": fname, "energy_eV": energy,
                                      "product": p.get("product"), "mass_ratio": p.get("mass_ratio")})
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"wrote {len(manifest['processes'])} processes ({n_exc} excitations) to {OUT}")


if __name__ == "__main__":
    main()
