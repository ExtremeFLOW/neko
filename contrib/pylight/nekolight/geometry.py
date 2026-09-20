# Copyright (c) 2026, The Neko Authors
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials provided
#     with the distribution.
#
#   * Neither the name of the authors nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
#     _  __  ____  __ __  ____
#    / |/ / / __/ / //_/ / __ \
#   /    / / _/  / ,<   / /_/ /
#  /_/|_/ /___/ /_/|_|  \____/
#
"""Hex-element geometry: the trilinear map at the 3x3x3 GLL nodes and its
Jacobian.  Everything is one einsum over element chunks; the Jacobian is
exact for straight-sided elements (curve records are not applied), the same
scope as Neko's mesh_checker.
"""

import numpy as np

from .formats import SIJK

# the 3-point GLL nodes on [-1, 1]
GLL3 = np.array([-1.0, 0.0, 1.0])
# corner signs of the 8 nmsh vertex slots (2 * SIJK - 1)
_RST = 2.0 * SIJK.astype(np.float64) - 1.0                     # (8, 3)

# precompute the 27 trilinear shape functions and their derivatives at the
# GLL nodes, in nek fld order (t outer, s, r inner)
_t, _s, _r = np.meshgrid(GLL3, GLL3, GLL3, indexing='ij')
_pts = np.stack([_r.ravel(), _s.ravel(), _t.ravel()], axis=1)  # (27, 3)

_one = 1.0 + _pts[:, None, :] * _RST[None, :, :]               # (27, 8, 3)
SHAPE_N = 0.125 * _one.prod(axis=2)                            # (27, 8)
SHAPE_D = np.empty((3, 27, 8))                                 # d/dr, d/ds, d/dt
for _c in range(3):
    _f = _one.copy()
    _f[:, :, _c] = _RST[None, :, _c]
    SHAPE_D[_c] = 0.125 * _f.prod(axis=2)
del _t, _s, _r, _pts, _one, _f, _c


def gll_xyz(corners):
    """GLL-node coordinates of the trilinear map: (m, 8, 3) -> (m, 27, 3)."""
    return np.einsum('pk,mkc->mpc', SHAPE_N, corners)


def min_jacobian(corners):
    """Minimum Jacobian determinant over the 27 GLL nodes per element:
    (m, 8, 3) -> (m,).  <= 0 flags an inverted/degenerate element."""
    J = np.einsum('dpk,mkc->mpdc', SHAPE_D, corners)   # (m, 27, 3, 3)
    det = (J[..., 0, 0] * (J[..., 1, 1] * J[..., 2, 2]
                           - J[..., 1, 2] * J[..., 2, 1])
           - J[..., 0, 1] * (J[..., 1, 0] * J[..., 2, 2]
                             - J[..., 1, 2] * J[..., 2, 0])
           + J[..., 0, 2] * (J[..., 1, 0] * J[..., 2, 1]
                             - J[..., 1, 1] * J[..., 2, 0]))
    return det.min(axis=1)


def facet_gll_mask():
    """(6, 27) bool: which of the 27 GLL nodes lie on facet 1..6."""
    idx = np.arange(27)
    ir, is_, it = idx % 3, (idx // 3) % 3, idx // 9
    m = np.zeros((6, 27), dtype=bool)
    m[0], m[1] = ir == 0, ir == 2
    m[2], m[3] = is_ == 0, is_ == 2
    m[4], m[5] = it == 0, it == 2
    return m
