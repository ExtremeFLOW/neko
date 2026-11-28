#!/usr/bin/env python3
"""Plot the initial and computed KPP density fields."""

from pathlib import Path
import argparse
import csv

import numpy as np
from PIL import Image, ImageDraw, ImageFont


RHO_MINIMUM = 0.25 * np.pi
RHO_MAXIMUM = 3.5 * np.pi


def load_plane():
    samples = {}
    paths = sorted(Path.cwd().glob("kpp_density_rank_*.csv"))
    if not paths:
        raise FileNotFoundError(
            "No kpp_density_rank_*.csv files found; run the case first."
        )

    for path in paths:
        with path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                x = round(float(row["x"]), 13)
                y = round(float(row["y"]), 13)
                key = (x, y)
                values = samples.setdefault(key, [0.0, 0.0, 0])
                values[0] += float(row["rho"])
                values[1] += float(row["rho_initial"])
                values[2] += 1

    x_values = sorted({key[0] for key in samples})
    y_values = sorted({key[1] for key in samples})
    x_index = {value: i for i, value in enumerate(x_values)}
    y_index = {value: i for i, value in enumerate(y_values)}
    computed = np.zeros((len(y_values), len(x_values)))
    initial = np.zeros_like(computed)
    for (x, y), (rho, rho_initial, count) in samples.items():
        computed[y_index[y], x_index[x]] = rho / count
        initial[y_index[y], x_index[x]] = rho_initial / count
    return np.asarray(x_values), np.asarray(y_values), initial, computed


def resample_uniform(x_values, y_values, values, size=720):
    """Interpolate a rectilinear GLL grid onto uniform physical coordinates."""
    x_uniform = np.linspace(x_values[0], x_values[-1], size)
    y_uniform = np.linspace(y_values[0], y_values[-1], size)
    x_interpolated = np.vstack(
        [np.interp(x_uniform, x_values, row) for row in values]
    )
    return np.vstack(
        [
            np.interp(y_uniform, y_values, x_interpolated[:, column])
            for column in range(size)
        ]
    ).T


def colorize(values, size=520):
    anchors = np.array(
        [
            [43, 28, 103],
            [37, 102, 172],
            [38, 181, 167],
            [239, 207, 75],
            [195, 52, 62],
        ],
        dtype=float,
    )
    scaled = np.clip(
        (values - RHO_MINIMUM) / (RHO_MAXIMUM - RHO_MINIMUM), 0.0, 1.0
    )
    position = scaled * (len(anchors) - 1)
    lower = np.floor(position).astype(int)
    upper = np.minimum(lower + 1, len(anchors) - 1)
    fraction = (position - lower)[..., None]
    rgb = anchors[lower] * (1.0 - fraction) + anchors[upper] * fraction
    image = Image.fromarray(np.uint8(rgb[::-1, :, :]))
    return image.resize((size, size), Image.Resampling.BILINEAR)


def main(elements, result_label):
    x_values, y_values, initial, computed = load_plane()
    initial = resample_uniform(x_values, y_values, initial)
    computed = resample_uniform(x_values, y_values, computed)
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
        (initial, "Initial KPP density"),
        (
            computed,
            f"{result_label}, {elements} x {elements} x 1",
        ),
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

    footer = (
        "Shown field: rho   |   color range: pi/4 to 3.5 pi"
        f"   |   computed range: {computed.min():.4f} to {computed.max():.4f}"
    )
    text_width = draw.textbbox((0, 0), footer, font=small_font)[2]
    draw.text(
        ((width - text_width) / 2, height - 43),
        footer,
        fill=(45, 45, 45),
        font=small_font,
    )
    canvas.save("kpp.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--elements", type=int, default=24)
    parser.add_argument(
        "--result-label",
        default="Euler IDP at t = 1, degree 5, C_R = 10",
    )
    arguments = parser.parse_args()
    main(arguments.elements, arguments.result_label)
