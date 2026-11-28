#!/usr/bin/env python3
"""Plot the exact and computed three-body density profiles."""

from pathlib import Path
import csv

import numpy as np
from PIL import Image, ImageDraw, ImageFont


def load_plane():
    samples = {}
    paths = sorted(Path.cwd().glob("three_body_density_rank_*.csv"))
    if not paths:
        raise FileNotFoundError(
            "No three_body_density_rank_*.csv files found; run the case first."
        )

    for path in paths:
        with path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                x = round(float(row["x"]), 13)
                y = round(float(row["y"]), 13)
                key = (x, y)
                values = samples.setdefault(key, [0.0, 0.0, 0])
                values[0] += float(row["rho"]) - 1.0
                values[1] += float(row["rho_exact"]) - 1.0
                values[2] += 1

    x_values = sorted({key[0] for key in samples})
    y_values = sorted({key[1] for key in samples})
    x_index = {value: i for i, value in enumerate(x_values)}
    y_index = {value: i for i, value in enumerate(y_values)}
    computed = np.zeros((len(y_values), len(x_values)))
    exact = np.zeros_like(computed)
    for (x, y), (rho, rho_exact, count) in samples.items():
        computed[y_index[y], x_index[x]] = rho / count
        exact[y_index[y], x_index[x]] = rho_exact / count
    return exact, computed


def colorize(values, size=520):
    anchors = np.array(
        [
            [35, 24, 65],
            [55, 74, 140],
            [35, 137, 153],
            [104, 187, 126],
            [238, 205, 86],
        ],
        dtype=float,
    )
    scaled = np.clip((values + 0.12) / 1.24, 0.0, 1.0)
    position = scaled * (len(anchors) - 1)
    lower = np.floor(position).astype(int)
    upper = np.minimum(lower + 1, len(anchors) - 1)
    fraction = (position - lower)[..., None]
    rgb = anchors[lower] * (1.0 - fraction) + anchors[upper] * fraction
    image = Image.fromarray(np.uint8(rgb[::-1, :, :]))
    return image.resize((size, size), Image.Resampling.BILINEAR)


def main():
    exact, computed = load_plane()
    panel_size = 520
    margin = 46
    title_height = 58
    footer_height = 72
    gap = 28
    width = 2 * panel_size + gap + 2 * margin
    height = panel_size + title_height + footer_height
    canvas = Image.new("RGB", (width, height), (249, 248, 244))
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default(size=18)
    small_font = ImageFont.load_default(size=15)

    panels = [
        (exact, "Exact profile after one turn"),
        (computed, "Euler IDP, degree 5, C_R = 0.1, 24 x 24 x 1 elements"),
    ]
    for panel, (values, title) in enumerate(panels):
        left = margin + panel * (panel_size + gap)
        top = title_height
        canvas.paste(colorize(values, panel_size), (left, top))
        draw.rectangle(
            (left, top, left + panel_size, top + panel_size),
            outline=(45, 45, 45),
            width=2,
        )
        text_width = draw.textbbox((0, 0), title, font=font)[2]
        draw.text(
            (left + (panel_size - text_width) / 2, 20),
            title,
            fill=(30, 30, 30),
            font=font,
        )

    footer = "Shown field: rho - 1   |   color range: -0.12 to 1.12"
    text_width = draw.textbbox((0, 0), footer, font=small_font)[2]
    draw.text(
        ((width - text_width) / 2, height - 43),
        footer,
        fill=(45, 45, 45),
        font=small_font,
    )
    canvas.save("three_body_rotation.png")


if __name__ == "__main__":
    main()
