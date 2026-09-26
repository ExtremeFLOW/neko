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
"""Mesh visualisation: extract the external surface (skin) as a quad mesh and
hand it to PyVista (interactive), a .vtu file (ParaView) or matplotlib (tiny
meshes).

The volume is never rendered: an n-element mesh has O(n^(2/3)) boundary
quads, so even a 3e8-hex mesh reduces to a few million quads -- comfortably
interactive on any workstation GPU once extracted.
"""

import base64
import struct
import sys

import numpy as np

from .formats import FACE_RE2, atomic_output
from .topology import skin


class Skin:
    """The external surface as an indexed quad mesh with per-face data."""

    def __init__(self, points, quads, celldata, elem_pos, facet):
        self.points = points        # (np, 3) f32
        self.quads = quads          # (nq, 4) int64 into points
        self.celldata = celldata    # name -> (nq,) array
        self.elem_pos = elem_pos    # (nq,) record position of the owner hex
        self.facet = facet          # (nq,) owner facet 0..5


def extract_skin(mesh, celldata_per_element=None):
    """Build the :class:`Skin` of a mesh (raw vertex ids, so periodic
    boundaries are included -- they are physical surfaces you want to see).

    ``celldata_per_element`` maps a name to an (nelv,) array; each skin quad
    inherits its owner element's value.
    """
    vidx = mesh.elems['v']['idx'].astype(np.int64)
    epos, fct = skin(vidx)
    slots = FACE_RE2[fct]                                # (nq, 4)
    ids = vidx[epos[:, None], slots]                     # (nq, 4)
    xyz = mesh.elems['v']['xyz'][epos[:, None], slots]   # (nq, 4, 3)
    uniq, inv = np.unique(ids.ravel(), return_inverse=True)
    points = np.empty((uniq.size, 3), dtype=np.float32)
    points[inv] = xyz.reshape(-1, 3).astype(np.float32)  # last write wins; ids agree
    quads = inv.reshape(-1, 4)
    data = {}
    if celldata_per_element:
        for name, arr in celldata_per_element.items():
            data[name] = np.asarray(arr)[epos]
    return Skin(points, quads, data, epos, fct)


# ---------------------------------------------------------------------------
# PyVista (interactive; ParaView's rendering engine via the vtk wheel)
# ---------------------------------------------------------------------------
def show_pyvista(sk, color_by=None, screenshot=None, title='meshview'):
    try:
        import pyvista as pv
    except ImportError:
        sys.exit('Error: interactive viewing needs pyvista '
                 '(pip install pyvista); use --export skin.vtu for ParaView '
                 'or --matplotlib for small meshes')
    nq = sk.quads.shape[0]
    cells = np.empty((nq, 5), dtype=np.int64)
    cells[:, 0] = 4
    cells[:, 1:] = sk.quads
    grid = pv.PolyData(sk.points.astype(np.float64), faces=cells.ravel())
    for name, arr in sk.celldata.items():
        grid.cell_data[name] = arr
    off = screenshot is not None
    p = pv.Plotter(off_screen=off, title=title)
    kw = {}
    if color_by and color_by in sk.celldata:
        kw = dict(scalars=color_by, cmap='tab20'
                  if np.issubdtype(sk.celldata[color_by].dtype, np.integer)
                  else 'viridis')
    p.add_mesh(grid, show_edges=nq < 200_000, **kw)
    p.add_axes()
    if off:
        p.show(screenshot=screenshot)
        print('  wrote %s' % screenshot)
    else:
        p.show()


# ---------------------------------------------------------------------------
# .vtu export (hand-rolled, base64-inline XML -- no extra dependencies)
# ---------------------------------------------------------------------------
def _b64(arr):
    raw = arr.tobytes()
    return base64.b64encode(struct.pack('<Q', len(raw)) + raw).decode('ascii')


def export_vtu(path, sk):
    """Write the skin as a .vtu (VTK XML unstructured grid of quads) that
    ParaView opens directly.  Inline base64, ~35 bytes per quad."""
    nq = sk.quads.shape[0]
    npts = sk.points.shape[0]
    conn = sk.quads.astype(np.int64).ravel()
    offs = (np.arange(1, nq + 1, dtype=np.int64) * 4)
    types = np.full(nq, 9, dtype=np.uint8)               # VTK_QUAD
    with atomic_output(path) as f:
        w = lambda s: f.write(s.encode('ascii'))
        w('<?xml version="1.0"?>\n')
        w('<VTKFile type="UnstructuredGrid" version="1.0" '
          'byte_order="LittleEndian" header_type="UInt64">\n')
        w('<UnstructuredGrid><Piece NumberOfPoints="%d" NumberOfCells="%d">\n'
          % (npts, nq))
        w('<Points><DataArray type="Float32" NumberOfComponents="3" '
          'format="binary">\n%s\n</DataArray></Points>\n'
          % _b64(sk.points.astype('<f4')))
        w('<Cells>\n')
        w('<DataArray type="Int64" Name="connectivity" format="binary">\n'
          '%s\n</DataArray>\n' % _b64(conn.astype('<i8')))
        w('<DataArray type="Int64" Name="offsets" format="binary">\n'
          '%s\n</DataArray>\n' % _b64(offs.astype('<i8')))
        w('<DataArray type="UInt8" Name="types" format="binary">\n'
          '%s\n</DataArray>\n' % _b64(types))
        w('</Cells>\n<CellData>\n')
        for name, arr in sk.celldata.items():
            if np.issubdtype(arr.dtype, np.integer):
                w('<DataArray type="Int32" Name="%s" format="binary">\n'
                  '%s\n</DataArray>\n' % (name, _b64(arr.astype('<i4'))))
            else:
                w('<DataArray type="Float32" Name="%s" format="binary">\n'
                  '%s\n</DataArray>\n' % (name, _b64(arr.astype('<f4'))))
        w('</CellData>\n</Piece></UnstructuredGrid></VTKFile>\n')


# ---------------------------------------------------------------------------
# matplotlib fallback (tiny meshes only)
# ---------------------------------------------------------------------------
def show_matplotlib(sk, color_by=None, screenshot=None):
    try:
        import matplotlib
        if screenshot:
            matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except ImportError:
        sys.exit('Error: --matplotlib needs matplotlib installed')
    nq = sk.quads.shape[0]
    if nq > 50_000:
        sys.exit('Error: %d boundary quads is too many for matplotlib -- '
                 'install pyvista or use --export' % nq)
    polys = sk.points[sk.quads]                          # (nq, 4, 3)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    pc = Poly3DCollection(polys, edgecolor='k', linewidths=0.2)
    if color_by and color_by in sk.celldata:
        vals = sk.celldata[color_by].astype(np.float64)
        pc.set_array(vals)
        pc.set_cmap('tab20')
        fig.colorbar(pc, ax=ax, label=color_by, shrink=0.7)
    else:
        pc.set_facecolor('lightsteelblue')
    ax.add_collection3d(pc)
    lo, hi = sk.points.min(axis=0), sk.points.max(axis=0)
    c, r = (lo + hi) / 2, (hi - lo).max() / 2 or 1.0
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_box_aspect((1, 1, 1))
    if screenshot:
        fig.savefig(screenshot, dpi=150, bbox_inches='tight')
        print('  wrote %s' % screenshot)
    else:
        plt.show()
