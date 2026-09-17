#!/usr/bin/env python3
"""Render existing field outputs with one shared, uncapped velocity scale."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

from render_velocity_field_heatmap import color_maximum, upscale_field
from sample_terrain_relative_velocity_npz import parse_fld_header


HERE = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fields", nargs="+", type=Path)
    parser.add_argument("--generated", type=Path, default=HERE / "generated_sharp_mask")
    parser.add_argument("--geometry-field", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="New output directory")
    parser.add_argument("--lift", type=float, default=20.0)
    parser.add_argument("--fps", type=float, default=5.0, help="Frames/s; 5 reproduces the latest fast movie")
    parser.add_argument("--gif", action="store_true", help="Also encode a GIF from the same frames")
    parser.add_argument("--width", type=int, default=1400, help="Movie/GIF width; PNGs remain full size")
    args = parser.parse_args()
    if not np.isfinite(args.fps) or args.fps <= 0:
        parser.error("--fps must be positive")
    if args.width < 2 or args.width % 2:
        parser.error("--width must be a positive even number")
    fields = sorted(args.fields, key=lambda p: parse_fld_header(p)["time"])
    times = [parse_fld_header(p)["time"] for p in fields]
    if len(set(times)) != len(times):
        parser.error("Select only one output per simulation time")
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        import imageio_ffmpeg
        ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    args.out.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env.pop("SAMPLE_Z", None)
    env.update(GENERATED_PATH=str(args.generated.resolve()),
               GEOMETRY_FIELD=str(args.geometry_field.resolve()), LIFT_M=str(args.lift))
    vmax = 0.0
    samples = []
    for i, path in enumerate(fields):
        sample = args.out / f"sample_{i:05d}.npz"
        env.update(FIELD_PATH=str(path.resolve()), NPZ_OUT=str(sample.resolve()))
        subprocess.run([sys.executable, str(HERE / "sample_terrain_relative_velocity_npz.py")],
                       env=env, check=True)
        with np.load(sample) as data:
            valid = data["disk"].astype(bool) & data["valid"].astype(bool)
            if not valid.any() or not np.isfinite(data["speed"][valid]).all():
                raise ValueError(f"Invalid velocity samples in {sample}")
            field, mask = upscale_field(data["speed"], data["disk"], data["valid"])
            vmax = max(vmax, color_maximum(data["speed"], valid, field, mask))
        samples.append(sample)
    env["VMAX"] = str(max(vmax, 1e-6))
    for i, sample in enumerate(samples):
        subprocess.run([sys.executable, str(HERE / "render_velocity_field_heatmap.py"),
                        str(args.generated), str(sample), str(args.out / f"frame_{i:05d}.png")],
                       env=env, check=True)
    subprocess.run([ffmpeg, "-nostdin", "-n", "-framerate", str(args.fps),
                    "-i", str(args.out / "frame_%05d.png"), "-c:v", "libx264",
                    "-crf", "16", "-threads", "2", "-pix_fmt", "yuv420p",
                    "-vf", f"scale={args.width}:-2:flags=lanczos", "-movflags", "+faststart",
                    str(args.out / "velocity.mp4")], check=True)
    if args.gif:
        filters = (f"scale={args.width}:-2:flags=lanczos,split[a][b];"
                   "[a]palettegen=reserve_transparent=0[p];[b][p]paletteuse=dither=bayer")
        subprocess.run([ffmpeg, "-nostdin", "-n", "-framerate", str(args.fps),
                        "-i", str(args.out / "frame_%05d.png"), "-filter_complex", filters,
                        "-filter_complex_threads", "1", "-loop", "0",
                        str(args.out / "velocity.gif")], check=True)
    (args.out / "manifest.json").write_text(json.dumps({
        "times": times, "lift_m": args.lift, "vmax": vmax, "fps": args.fps,
        "width": args.width,
        "colormap": "Cool to Warm", "color_bands": 32,
    }, indent=2) + "\n")
    print(f"Rendered {times[0]:g} to {times[-1]:g} s; shared scale 0 to {vmax:g} m/s")


if __name__ == "__main__":
    main()
