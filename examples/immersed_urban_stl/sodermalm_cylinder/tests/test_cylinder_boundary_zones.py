from pathlib import Path
import re
import struct
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_sodermalm_cylinder as mesh
from check_cylinder_boundary_zones import check_boundary_zones


class BoundaryZonesTest(unittest.TestCase):
    def make_mesh(self, folder, width, bearing=225.0):
        with patch.multiple(mesh, INFLOW_ARC_WIDTH_DEG=width, INFLOW_FROM_DEG=bearing, NZ=2):
            geo = folder / "disk.geo"
            mesh.write_gmsh_geo(geo, 10.0, [], (0.0, 0.0))
            text = geo.read_text()
            outer = {}
            for idx, x, y in re.findall(r"Point\((\d+)\) = \{([^,]+), ([^,]+), 0,", text):
                if int(idx) != 1:
                    outer[int(idx)] = (float(x), float(y))
            nodes = dict(outer)
            nodes.update({k + 1000: (x * .5, y * .5) for k, (x, y) in outer.items()})
            labels = {}
            for label, curves in re.findall(r"Physical Curve\((\d+)\) = \{([^}]+)\}", text):
                labels.update({int(c): int(label) for c in curves.split(",")})
            quads, lines = [], []
            for idx, a, b in re.findall(r"Circle\((\d+)\) = \{(\d+), 1, (\d+)\}", text):
                idx, a, b = int(idx), int(a), int(b)
                quads.append((idx, [a + 1000, a, b, b + 1000]))
                lines.append((labels[idx], a, b))
            path = folder / "test.nmsh"
            with patch.object(mesh, "terrain_height", return_value=0.0):
                info = mesh.write_nmsh(path, nodes, quads, lines, np.empty((0, 3)), 0.0, [], (0., 0.), 10.)
        return path, info["boundary_validation"]

    def test_exact_endpoints(self):
        for width, bearing in [(180.,225.), (220.,225.), (170.,225.), (180.,0.), (220.,315.)]:
            with self.subTest(width=width, bearing=bearing), tempfile.TemporaryDirectory() as d:
                _, info = self.make_mesh(Path(d), width, bearing)
                self.assertAlmostEqual(info["inlet_width_degrees"], width, places=6)
                expected = sorted([(bearing - width / 2) % 360, (bearing + width / 2) % 360])
                np.testing.assert_allclose(info["inlet_endpoints_degrees"], expected, atol=1e-6)

    def test_old_180_mesh_rejected_when_220_requested(self):
        with tempfile.TemporaryDirectory() as d:
            path, _ = self.make_mesh(Path(d), 180.)
            with self.assertRaisesRegex(ValueError, "mismatch|crosses|differs"):
                check_boundary_zones(path, 225., 220., 90.)

    def test_corrupted_binary_label_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path, _ = self.make_mesh(Path(d), 220.)
            with path.open("r+b") as f:
                ne, = struct.unpack("<i", f.read(4))
                f.seek(8 + ne * 228 + 4 + 3 * 4)
                f.write(struct.pack("<i", 7))
            with self.assertRaisesRegex(ValueError, "Invalid cylindrical side label"):
                check_boundary_zones(path, 225., 220., 90.)

    def test_invalid_arc_step(self):
        for step in [0., -45., 180., 17.]:
            with patch.object(mesh, "CIRCLE_ARC_STEP_DEG", step), tempfile.TemporaryDirectory() as d:
                with self.assertRaises(ValueError):
                    mesh.write_gmsh_geo(Path(d) / "bad.geo", 10., [], (0., 0.))


if __name__ == "__main__":
    unittest.main()
