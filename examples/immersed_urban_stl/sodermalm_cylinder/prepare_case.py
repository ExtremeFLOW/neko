#!/usr/bin/env python3
"""Build or verify the p7 sharp building cache for the corrected cylinder."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from check_cylinder_boundary_zones import check_boundary_zones


HERE = Path(__file__).resolve().parent
GENERATED = HERE / "generated_sharp_mask"
MESH = GENERATED / "sodermalm_cylinder_sharp_mask.nmsh"
CACHE = HERE / "cache_sharp_mask"
MANIFEST = CACHE / "inputs.json"


def fingerprint() -> dict:
    case = json.loads((HERE / "run" / "sharp_mask.case").read_text())["case"]
    if case["numerics"]["polynomial_order"] != 7:
        raise ValueError("This cache preparation requires polynomial order 7")
    if (HERE / "run" / case["mesh_file"]).resolve() != MESH:
        raise ValueError("Case does not reference the sharp-mask mesh")
    check_boundary_zones(MESH, 225.0, 220.0, 90.0)
    paths = [MESH, *sorted(GENERATED.glob("*.geojson")),
             *sorted(GENERATED.glob("*metadata.json")),
             HERE / "make_template_field_from_nmsh.py",
             HERE / "make_sodermalm_building_cache.py"]
    hashes = {}
    for path in paths:
        with path.open("rb") as handle:
            hashes[str(path.relative_to(HERE))] = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"polynomial_order": 7, "smooth_width_m": 0, "max_height_m": 65,
            "sha256": hashes}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify preparation without writing files")
    args = parser.parse_args()
    expected = fingerprint()
    if MANIFEST.exists():
        if json.loads(MANIFEST.read_text()) != expected:
            raise RuntimeError("Cache inputs changed; regenerate into a fresh cache_sharp_mask directory")
        for name in ("buildings0.f00000", "buildings0.nek5000"):
            if not (CACHE / name).is_file():
                raise FileNotFoundError(CACHE / name)
        print("Verified p7 cache inputs and the actual 220-degree inlet")
        return
    if args.check:
        raise RuntimeError("Run prepare_case.py before launching the solver")
    if CACHE.exists() and any(CACHE.iterdir()):
        raise RuntimeError("Refusing to reuse an incomplete or unverified cache_sharp_mask directory")
    CACHE.mkdir(parents=True, exist_ok=True)
    template = CACHE / "template.f00000"
    subprocess.run([sys.executable, str(HERE / "make_template_field_from_nmsh.py"),
                    str(MESH), str(template), "--polynomial-order", "7"], check=True)
    subprocess.run([sys.executable, str(HERE / "make_sodermalm_building_cache.py"),
                    "--generated", str(GENERATED), "--template-field", str(template),
                    "--cache-prefix", str(CACHE / "buildings"), "--smooth-width", "0",
                    "--max-height", "65"], check=True)
    MANIFEST.write_text(json.dumps(expected, indent=2) + "\n")
    print("Prepared sharp-mask cache")


if __name__ == "__main__":
    main()
