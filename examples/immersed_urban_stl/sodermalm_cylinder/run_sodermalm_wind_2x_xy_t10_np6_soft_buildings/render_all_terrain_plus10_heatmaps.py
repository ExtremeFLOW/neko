#!/usr/bin/env python3
"""Render terrain+10 m velocity heatmaps for every field snapshot."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
PYTHON = Path(sys.executable)
RENDER = HERE / "render_terrain_plus10_heatmap.py"
FIELDS = HERE / "fields"
OUT = HERE / "renders"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fields = sorted(FIELDS.glob("field0.f*"))
    if not fields:
        raise SystemExit(f"No field snapshots found in {FIELDS}")

    for field in fields:
        suffix = field.suffix[1:]
        out = OUT / f"sodermalm_velocity_terrain_plus10_{suffix}.png"
        if out.exists() and out.stat().st_mtime >= field.stat().st_mtime:
            print(f"Skipping up-to-date {out.name}", flush=True)
            continue
        env = os.environ.copy()
        env["FIELD_PATH"] = str(field)
        env["OUT_PATH"] = str(out)
        print(f"Rendering {field.name} -> {out.name}", flush=True)
        subprocess.run([str(PYTHON), str(RENDER)], check=True, env=env, cwd=HERE)


if __name__ == "__main__":
    main()
