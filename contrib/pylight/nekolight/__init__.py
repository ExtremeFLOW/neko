# Copyright (c) 2026, The Neko Authors
# All rights reserved.  BSD-3-Clause: see any module in this package for the
# full licence text.
#
#     _  __  ____  __ __  ____
#    / |/ / / __/ / //_/ / __ \
#   /    / / _/  / ,<   / /_/ /
#  /_/|_/ /___/ /_/|_|  \____/
#
"""nekolight -- the shared library behind Neko's Python light tools.

Plain functions over numpy arrays; the only data structure is the
:class:`~nekolight.formats.Mesh` named tuple.  Core dependencies are numpy
(everywhere) and scipy (partitioning only); pymetis, pyamg, pyvista and
matplotlib are optional extras.
"""

from .formats import (Mesh, EL_DT, ZONE_DT, CURVE_DT, FACE_RE2, EDGE_RE2,
                      FACET_MAP, SIJK, MAX_ZLBLS, banner, atomic_output,
                      read_nmsh, iter_nmsh_elements, write_nmsh,
                      validate_zones, validate_curves, read_re2,
                      read_dist_csv, write_zone_indices_fld, bc_type_str)
from .topology import (pos_of_elid_map, dedup_points, periodic_min_merge,
                       periodic_replace_merge, compress_ids,
                       face_multiplicity, count_edges, skin, dual_graph)
from .geometry import gll_xyz, min_jacobian, facet_gll_mask, GLL3
from .partition import (neko_linear_sizes, spectral_partition,
                        metis_partition, geometric_partition, repair_sizes,
                        weighted_cut, reorder_and_write)

__all__ = [n for n in dir() if not n.startswith('_')]
__version__ = '0.1.0'
