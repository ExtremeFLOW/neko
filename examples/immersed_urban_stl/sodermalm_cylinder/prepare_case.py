#!/usr/bin/env python3
"""Build or verify the p9 sharp building cache and same-mesh restart inputs."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys

from check_cylinder_boundary_zones import check_boundary_zones
from sample_terrain_relative_velocity_npz import parse_fld_header


HERE = Path(__file__).resolve().parent
GENERATED = HERE / "generated_sharp_mask"
MESH = GENERATED / "sodermalm_cylinder_sharp_mask.nmsh"
CACHE = HERE / "cache_sharp_mask"
MANIFEST = CACHE / "inputs.json"
BASELINE_MESH_SHA256 = "846b42f0fb4c8de0d196388e7f2dce77a7cd9b1850b5a7a25f5e8a7893c0098d"


def check_restart(case: dict, mesh_hash: str) -> None:
    if "restart_file" not in case:
        return
    if mesh_hash != BASELINE_MESH_SHA256:
        raise ValueError("The baseline restart requires the original mesh, not regenerated geometry")
    path = HERE / "run" / case["restart_file"]
    with path.open("rb") as handle:
        nel, dim, lx, flags, time = struct.unpack("<4id", handle.read(24))
    if (nel, dim, flags) != (347160, 3, 13) or lx not in (8, 10):
        raise ValueError("Expected a same-mesh p7 or p9 double-precision fluid checkpoint")
    # Four primary fields, six velocity histories, six AB histories, and time arrays.
    if path.stat().st_size != 184 + 16 * nel * lx**3 * 8:
        raise ValueError("Incomplete checkpoint: unexpected file size")
    if not math.isfinite(time) or not 0 <= time < case["time"]["end_time"]:
        raise ValueError("Checkpoint time must precede the requested end time")


def check_cache(path: Path) -> None:
    header = parse_fld_header(path)
    with MESH.open("rb") as handle:
        nel, = struct.unpack("<i", handle.read(4))
    if (header["lx"], header["ly"], header["lz"], header["nelv"],
            header["rdcode"], header["wdsz"]) != (10, 10, 10, nel, "XS01", 8):
        raise ValueError("Cache must contain p9 coordinates and the sharp scalar mask")
    if path.stat().st_size != 136 + 4 * nel + nel * 1000 * 4 * 8:
        raise ValueError("Incomplete p9 building cache")


def fingerprint() -> dict:
    case = json.loads((HERE / "run" / "sharp_mask.case").read_text())["case"]
    if case["numerics"]["polynomial_order"] != 9:
        raise ValueError("This cache preparation requires polynomial order 9")
    if (HERE / "run" / case["mesh_file"]).resolve() != MESH.resolve():
        raise ValueError("Case does not reference the sharp-mask mesh")
    obj = case["fluid"]["source_terms"][0]["objects"][0]
    if (HERE / "run" / obj["file_name"]).resolve() != (CACHE / "buildings0.fld").resolve():
        raise ValueError("Case does not reference the verified building cache")
    check_boundary_zones(MESH, 225.0, 180.0, 90.0)
    paths = [MESH, *sorted(GENERATED.glob("*.geojson")),
             *sorted(GENERATED.glob("*metadata.json")),
             HERE / "make_template_field_from_nmsh.py",
             HERE / "make_sodermalm_building_cache.py"]
    hashes = {}
    for path in paths:
        with path.open("rb") as handle:
            hashes[str(path.relative_to(HERE))] = hashlib.file_digest(handle, "sha256").hexdigest()
    check_restart(case, hashes[str(MESH.relative_to(HERE))])
    return {"polynomial_order": 9, "smooth_width_m": 0, "max_height_m": 65,
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
        check_cache(CACHE / "buildings0.f00000")
        print("Verified p9 cache inputs and the actual 180-degree inlet")
        return
    if args.check:
        raise RuntimeError("Run prepare_case.py before launching the solver")
    if CACHE.exists() and any(CACHE.iterdir()):
        raise RuntimeError("Refusing to reuse an incomplete or unverified cache_sharp_mask directory")
    CACHE.mkdir(parents=True, exist_ok=True)
    template = CACHE / "template.f00000"
    subprocess.run([sys.executable, str(HERE / "make_template_field_from_nmsh.py"),
                    str(MESH), str(template), "--polynomial-order", "9"], check=True)
    subprocess.run([sys.executable, str(HERE / "make_sodermalm_building_cache.py"),
                    "--generated", str(GENERATED), "--template-field", str(template),
                    "--cache-prefix", str(CACHE / "buildings"), "--smooth-width", "0",
                    "--max-height", "65"], check=True)
    check_cache(CACHE / "buildings0.f00000")
    MANIFEST.write_text(json.dumps(expected, indent=2) + "\n")
    print("Prepared sharp-mask cache")


if __name__ == "__main__":
    main()
