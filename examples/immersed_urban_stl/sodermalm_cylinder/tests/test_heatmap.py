from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import render_velocity_field_heatmap as render


class HeatmapTest(unittest.TestCase):
    def test_32_colors_and_log_mapping(self):
        values = .025 * np.expm1(np.linspace(0, np.log1p(6 / .025), 1000))
        colors = render.color_velocity(values, 6)
        self.assertEqual(len(np.unique(colors, axis=0)), 32)
        self.assertAlmostEqual(float(render.velocity_fraction(0, 6)), 0)
        self.assertAlmostEqual(float(render.velocity_fraction(6, 6)), 1)
        self.assertGreater(float(render.velocity_fraction(.1, 6)), .1 / 6)

    def test_no_upper_clipping(self):
        speed = np.array([[1., 30.]])
        valid = np.ones_like(speed, dtype=bool)
        self.assertEqual(render.color_maximum(speed, valid, speed, valid), 30)
        with self.assertRaisesRegex(ValueError, "hide"):
            render.color_maximum(speed, valid, speed, valid, "6")
        with self.assertRaises(ValueError):
            render.color_maximum(speed, valid, speed, valid, "nan")

    def test_courtyard_and_overlapping_buildings(self):
        outer = [[-8, -8], [8, -8], [8, 8], [-8, 8], [-8, -8]]
        hole = [[-3, -3], [3, -3], [3, 3], [-3, 3], [-3, -3]]
        inner = [[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, -1]]
        feature = {"geometry": {"type": "Polygon", "coordinates": [outer, hole]}}
        with patch.object(render, "CANVAS_N", 101):
            mask = np.asarray(render.building_mask([feature], (0, 0), 10))
            self.assertEqual(mask[50, 50], 0)
            self.assertEqual(mask[50, 20], 255)
            self.assertEqual(mask[0, 0], 0)
            second = {"geometry": {"type": "MultiPolygon", "coordinates": [[inner]]}}
            mask = np.asarray(render.building_mask([second, feature], (0, 0), 10))
            self.assertEqual(mask[50, 50], 255)
            self.assertEqual(mask[50, 40], 0)

    def test_upscale_keeps_positive_extreme(self):
        speed = np.array([[0., 0.], [0., 100.]])
        valid = np.ones_like(speed, dtype=bool)
        with patch.object(render, "CANVAS_N", 64):
            field, mask = render.upscale_field(speed, valid, valid)
        self.assertGreaterEqual(field.min(), 0)
        self.assertGreaterEqual(render.color_maximum(speed, valid, field, mask), 100)


if __name__ == "__main__":
    unittest.main()
