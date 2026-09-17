import json
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import build_sodermalm_cylinder as geometry
import make_sodermalm_building_cache as cache


class SharpMaskSettingsTest(unittest.TestCase):
    def test_geometry_and_case_defaults(self):
        case = json.loads((HERE / "run" / "sharp_mask.case").read_text())["case"]
        self.assertEqual(geometry.INFLOW_ARC_WIDTH_DEG, 180)
        self.assertEqual(geometry.INFLOW_FROM_DEG, 225)
        self.assertEqual(geometry.NZ, 10)
        self.assertEqual(geometry.DOMAIN_HEIGHT_ABOVE_LOWEST_M, 350)
        self.assertEqual(case["numerics"]["polynomial_order"], 9)
        self.assertAlmostEqual(geometry.LAND_MESH_SIZE_M / 9, 3.8888888889)
        self.assertTrue(case["numerics"]["dealias"])
        self.assertEqual(case["restart_file"], "restart.chkp")
        self.assertEqual(case["time"]["end_time"], 500)
        self.assertEqual(case["fluid"]["output_value"], 10)
        self.assertEqual(case["checkpoint_value"], 25)
        zones = {z: bc["type"] for bc in case["fluid"]["boundary_conditions"]
                 for z in bc["zone_indices"]}
        self.assertEqual(zones, {1: "user_velocity", 2: "outflow+dong",
                                 3: "outflow+dong", 5: "no_slip", 6: "normal_outflow"})
        brinkman, hpfrt = case["fluid"]["source_terms"]
        self.assertEqual(brinkman["type"], "brinkman")
        self.assertEqual(brinkman["filter"]["radius"], 10)
        self.assertEqual(brinkman["penalty"], 10)
        self.assertEqual(brinkman["ramp_time"], 30)
        self.assertNotIn("implicit", brinkman)
        self.assertEqual(hpfrt, {"type": "hpfrt", "filter_modes": 2, "filter_weight": 2})

    def test_sharp_raw_mask(self):
        distance = np.array([-10., -1e-12, 0., 1e-12, 10.])
        np.testing.assert_array_equal(cache.smooth_step(distance, 0), [1, 1, 1, 0, 0])

    def test_geometry_and_cache_building_defaults_match(self):
        for name in ("BUILDING_BASE_EMBED_M", "BUILDING_CONTACT_SAMPLE_SPACING_M",
                     "BUILDING_MIN_FOOTPRINT_AREA_M2", "BUILDING_MIN_EDGE_M",
                     "BUILDING_SIMPLIFY_M", "BUILDING_MIN_HEIGHT_M", "BUILDING_MAX_HEIGHT_M",
                     "SHORE_BLEND_M"):
            self.assertEqual(getattr(geometry, name), getattr(cache, name), name)
        self.assertEqual(cache.BUILDING_MIN_FOOTPRINT_AREA_M2, 0)
        self.assertEqual(cache.BUILDING_MIN_EDGE_M, 0)


if __name__ == "__main__":
    unittest.main()
