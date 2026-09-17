from pathlib import Path
import struct
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import make_sodermalm_building_cache as cache
import prepare_case as prepare


class P9PreparationTest(unittest.TestCase):
    def test_cache_order_and_completeness(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mesh = root / "test.nmsh"
            mesh.write_bytes(struct.pack("<2i", 1, 3))
            for lx in (8, 10):
                header = dict(lx=lx, ly=lx, lz=lx, nelv=1, nelgt=1)
                zeros = np.zeros(lx**3)
                cache.write_cache(root / "buildings", header, np.array([1]),
                                  zeros, zeros, zeros, zeros)
                field = root / "buildings0.f00000"
                with patch.object(prepare, "MESH", mesh):
                    if lx == 8:
                        with self.assertRaisesRegex(ValueError, "p9"):
                            prepare.check_cache(field)
                    else:
                        prepare.check_cache(field)
                        with field.open("r+b") as handle:
                            handle.truncate(136)
                        with self.assertRaisesRegex(ValueError, "Incomplete"):
                            prepare.check_cache(field)

    def test_restart_requires_original_mesh(self):
        with self.assertRaisesRegex(ValueError, "original mesh"):
            prepare.check_restart({"restart_file": "restart.chkp"}, "regenerated")

    def test_restart_order_size_and_time(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "run").mkdir()
            path = root / "run" / "restart.chkp"
            case = {"restart_file": "restart.chkp", "time": {"end_time": 500}}
            for lx in (8, 10):
                path.write_bytes(struct.pack("<4id", 347160, 3, lx, 13, 225.0))
                size = 184 + 16 * 347160 * lx**3 * 8
                with patch.object(prepare, "HERE", root):
                    with patch.object(Path, "stat", return_value=SimpleNamespace(st_size=size)):
                        prepare.check_restart(case, prepare.BASELINE_MESH_SHA256)
                        with self.assertRaisesRegex(ValueError, "time"):
                            prepare.check_restart({**case, "time": {"end_time": 200}},
                                                  prepare.BASELINE_MESH_SHA256)
                    with self.assertRaisesRegex(ValueError, "Incomplete"):
                        prepare.check_restart(case, prepare.BASELINE_MESH_SHA256)


if __name__ == "__main__":
    unittest.main()
