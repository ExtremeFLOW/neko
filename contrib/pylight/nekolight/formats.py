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
"""Every binary layout and topology table used by the Neko light tools.

This module is the single home for the format knowledge shared by all the
tools; nothing else in the package hard-codes a byte offset or a
vertex-ordering table.

The Neko binary mesh (`.nmsh`, little-endian, no record markers)::

    header    2 x int32                      nelv, gdim
    element   int32 + 8 x (int32 + 3 x f64)  el_idx, 8 x (v_idx, x, y, z)   228 B
    (nzones)  int32
    zone      9 x int32                       36 B
              e, f, p_e, p_f, glb_pt_ids(4), type   (5 = periodic, 7 = labelled)
    (ncurves) int32
    curve     int32 + 60 x f64 + 12 x int32  el, curve_data(5,12), type(12) 532 B

Real Neko meshes may carry trailing bytes past the curve section (an MPI-IO
no-truncate artifact); Neko's reader ignores them and so do we (the checker
reports them).

The NEKTON re2 mesh (little-endian)::

    header    80 ASCII chars: '#v001'..'#v004' + counts (fixed widths below)
    endian    f32 = 6.54321  (byte-swapped files are rejected)
    element   v2+: f64 group + 8 x f64 x/y/z          200 B   (v1: f32, 100 B)
    (ncurve)  v2+: f64  (v1: int32)
    curve     v2+: 2 x f64 + 5 x f64 + char(8)         64 B   (v1: 32 B)
    (nbc)     v2+: f64  (v1: int32)
    bc        v2+: 2 x f64 + 5 x f64 + char(8)         64 B   (v1: 32 B)

Topology tables (all validated byte-exact against Neko's own writers by the
golden tests): ``FACE_RE2`` gives the 4 corners of Neko facet 1..6 as nmsh
vertex slots -- Neko's ``face_nodes`` composed with the nmsh->mesh read swap
[1,2,4,3,5,6,8,7], which nets to identity on the file layout.  ``EDGE_RE2``
likewise for the 12 hex edges, ``FACET_MAP`` maps an re2 face to Neko's
symmetric facet, and ``SIJK`` maps an nmsh vertex slot to a corner of the
reference cube.
"""

import os
import sys
import tempfile
from typing import NamedTuple

import numpy as np

# ---------------------------------------------------------------------------
# nmsh record dtypes
# ---------------------------------------------------------------------------
EL_DT = np.dtype([('id', '<i4'),
                  ('v', [('idx', '<i4'), ('xyz', '<f8', (3,))], (8,))])
ZONE_DT = np.dtype([('e', '<i4'), ('f', '<i4'), ('p_e', '<i4'), ('p_f', '<i4'),
                    ('g', '<i4', (4,)), ('t', '<i4')])
CURVE_DT = np.dtype([('e', '<i4'), ('data', '<f8', (12, 5)),
                     ('type', '<i4', (12,))])
assert EL_DT.itemsize == 228 and ZONE_DT.itemsize == 36 \
    and CURVE_DT.itemsize == 532

# ---------------------------------------------------------------------------
# re2 record dtypes (v2+ = double precision body, v1 = single precision)
# ---------------------------------------------------------------------------
RE2_EL_DT = {True: np.dtype([('rg', '<f8'), ('x', '<f8', (8,)),
                             ('y', '<f8', (8,)), ('z', '<f8', (8,))]),
             False: np.dtype([('rg', '<f4'), ('x', '<f4', (8,)),
                              ('y', '<f4', (8,)), ('z', '<f4', (8,))])}
RE2_CURVE_DT = {True: np.dtype([('e', '<f8'), ('edge', '<f8'),
                                ('d', '<f8', (5,)), ('t', 'S8')]),
                False: np.dtype([('e', '<i4'), ('edge', '<i4'),
                                 ('d', '<f4', (5,)), ('t', 'S4')])}
RE2_BC_DT = {True: np.dtype([('e', '<f8'), ('f', '<f8'),
                             ('d', '<f8', (5,)), ('t', 'S8')]),
             False: np.dtype([('e', '<i4'), ('f', '<i4'),
                              ('d', '<f4', (5,)), ('t', 'S4')])}
RE2_ENDIAN_TEST = np.float32(6.54321)

# ---------------------------------------------------------------------------
# Topology tables (0-based numpy versions of the validated Fortran tables)
# ---------------------------------------------------------------------------
# corners of Neko facet 1..6 as nmsh vertex slots
FACE_RE2 = np.array([[1, 5, 8, 4], [2, 6, 7, 3], [1, 2, 6, 5],
                     [4, 3, 7, 8], [1, 2, 3, 4], [5, 6, 7, 8]],
                    dtype=np.int64) - 1
# the 12 hex edges as nmsh vertex slot pairs
EDGE_RE2 = np.array([[1, 2], [3, 4], [5, 6], [7, 8], [1, 4], [2, 3],
                     [5, 8], [6, 7], [1, 5], [2, 6], [4, 8], [3, 7]],
                    dtype=np.int64) - 1
# re2 face (1..6) -> Neko symmetric facet (1..6); index with [face - 1]
FACET_MAP = np.array([3, 2, 4, 1, 5, 6], dtype=np.int64)
# nmsh vertex slot (0..7) -> (ix, iy, iz) corner of the reference cube
SIJK = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                 [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=np.int64)
# labels outside [1, NEKO_MSH_MAX_ZLBLS] are rejected, exactly as Neko does
MAX_ZLBLS = 20

BANNER = r"""
     _  __  ____  __ __  ____
    / |/ / / __/ / //_/ / __ \
   /    / / _/  / ,<   / /_/ /
  /_/|_/ /___/ /_/|_|  \____/   %s
"""


def banner(name):
    """The NEKO banner with a tool name attached."""
    return BANNER % name


class Mesh(NamedTuple):
    """A whole .nmsh in memory: the raw structured record arrays.

    ``elems['id']`` are the global element ids, ``elems['v']['idx']`` the
    (nelv, 8) vertex ids and ``elems['v']['xyz']`` the (nelv, 8, 3) corner
    coordinates.  Zones and curves are kept in file layout so a write is a
    plain ``tofile``.
    """
    nelv: int
    elems: np.ndarray     # EL_DT (nelv,)
    zones: np.ndarray     # ZONE_DT
    curves: np.ndarray    # CURVE_DT
    trailing: int = 0     # bytes past the curve section (MPI-IO artifact)


# ---------------------------------------------------------------------------
# Atomic output writing (never clobber inputs, never leave partial files)
# ---------------------------------------------------------------------------
class atomic_output:
    """Write to a temporary file and rename it into place only on success.

    Guards against ``out == in`` (which would truncate the input) and
    guarantees no partial output survives an error -- any exception unlinks
    the temporary file and the destination is never touched.
    """

    def __init__(self, path, inputs=()):
        rp = os.path.realpath(path)
        for src in inputs:
            if os.path.exists(src) and os.path.realpath(src) == rp:
                sys.exit('Error: output %s would overwrite an input file'
                         % path)
        self.path = path
        d = os.path.dirname(rp) or '.'
        fd, self.tmp = tempfile.mkstemp(prefix=os.path.basename(path) + '.',
                                        suffix='.tmp', dir=d)
        self.f = os.fdopen(fd, 'wb')

    def __enter__(self):
        return self.f

    def __exit__(self, exc_type, exc, tb):
        self.f.close()
        if exc_type is None:
            os.replace(self.tmp, self.path)
        else:
            try:
                os.unlink(self.tmp)
            except OSError:
                pass
        return False


# ---------------------------------------------------------------------------
# nmsh read / write
# ---------------------------------------------------------------------------
def _count(f, what, path):
    a = np.fromfile(f, dtype='<i4', count=1)
    if a.size != 1:
        sys.exit('Error: %s: truncated file (missing %s count)' % (path, what))
    n = int(a[0])
    if n < 0:
        sys.exit('Error: %s: negative %s count (%d) -- corrupt file?'
                 % (path, what, n))
    return n


def read_nmsh(path):
    """Read a whole .nmsh into a :class:`Mesh` (validating as it goes)."""
    try:
        f = open(path, 'rb')
    except OSError as ex:
        sys.exit('Error: cannot open %s (%s)' % (path, ex.strerror))
    with f:
        hdr = np.fromfile(f, dtype='<i4', count=2)
        if hdr.size != 2:
            sys.exit('Error: %s is not a Neko .nmsh file (short header)'
                     % path)
        nelv, gdim = int(hdr[0]), int(hdr[1])
        if gdim != 3:
            sys.exit('Error: the light tools support 3D (hex) meshes only '
                     '(gdim=%d)' % gdim)
        if nelv < 1:
            sys.exit('Error: %s: non-positive element count (%d)'
                     % (path, nelv))
        elems = np.fromfile(f, dtype=EL_DT, count=nelv)
        if elems.size != nelv:
            sys.exit('Error: %s: truncated element section (%d of %d records)'
                     % (path, elems.size, nelv))
        nz = _count(f, 'zone', path)
        zones = np.fromfile(f, dtype=ZONE_DT, count=nz)
        if zones.size != nz:
            sys.exit('Error: %s: truncated zone section (%d of %d records)'
                     % (path, zones.size, nz))
        nc = _count(f, 'curve', path)
        curves = np.fromfile(f, dtype=CURVE_DT, count=nc)
        if curves.size != nc:
            sys.exit('Error: %s: truncated curve section (%d of %d records)'
                     % (path, curves.size, nc))
        here = f.tell()
        f.seek(0, os.SEEK_END)
        trailing = f.tell() - here
    return Mesh(nelv, elems, zones, curves, trailing)


def iter_nmsh_elements(path, chunk=1 << 21):
    """Yield (start, elems_chunk) over the element section of a .nmsh.

    The low-memory alternative to :func:`read_nmsh` for reductions that never
    need the whole mesh at once (bounding boxes, Jacobian scans, ...).
    """
    with open(path, 'rb') as f:
        hdr = np.fromfile(f, dtype='<i4', count=2)
        if hdr.size != 2 or int(hdr[1]) != 3:
            sys.exit('Error: %s is not a 3D Neko .nmsh file' % path)
        nelv = int(hdr[0])
        done = 0
        while done < nelv:
            n = min(chunk, nelv - done)
            e = np.fromfile(f, dtype=EL_DT, count=n)
            if e.size != n:
                sys.exit('Error: %s: truncated element section '
                         '(%d of %d records)' % (path, done + e.size, nelv))
            yield done, e
            done += n


def write_nmsh(path, elems, zone_arrays, curves, inputs=()):
    """Write a .nmsh atomically.  ``zone_arrays`` is a sequence of ZONE_DT
    arrays written in order (periodic first, then labelled, then any legacy
    types -- the order Neko's own writer uses)."""
    nz = sum(int(z.shape[0]) for z in zone_arrays)
    with atomic_output(path, inputs) as f:
        np.array([elems.shape[0], 3], dtype='<i4').tofile(f)
        elems.tofile(f)
        np.array([nz], dtype='<i4').tofile(f)
        for z in zone_arrays:
            z.tofile(f)
        np.array([curves.shape[0]], dtype='<i4').tofile(f)
        curves.tofile(f)


def validate_zones(nelv, zones, path=''):
    """Refuse malformed zone records instead of silently repairing them.

    Checks: plausible type, element/facet ranges, periodic partner ranges and
    labelled-zone label range.  Any violation is a hard error -- no tool in
    this package drops or rewrites a bad record while exiting 0.
    """
    where = (' in %s' % path) if path else ''
    if zones.size == 0:
        return
    t = zones['t']
    if ((t < 1) | (t > 7)).any():
        sys.exit('Error: zone record with implausible type (valid: 1..7)%s '
                 '-- corrupt or mis-framed zone section?' % where)
    bad = (zones['e'] < 1) | (zones['e'] > nelv) \
        | (zones['f'] < 1) | (zones['f'] > 6)
    if bad.any():
        i = int(np.flatnonzero(bad)[0])
        sys.exit('Error: zone record %d references element %d facet %d, '
                 'outside [1,%d]x[1,6]%s'
                 % (i + 1, int(zones['e'][i]), int(zones['f'][i]), nelv,
                    where))
    z5 = zones[t == 5]
    if z5.size:
        bad = (z5['p_e'] < 1) | (z5['p_e'] > nelv) \
            | (z5['p_f'] < 1) | (z5['p_f'] > 6)
        if bad.any():
            sys.exit('Error: periodic zone record references partner '
                     'element/facet out of range%s' % where)
    z7 = zones[t == 7]
    if z7.size:
        lbl = z7['p_f']       # the label lives in the p_f field
        if ((lbl < 1) | (lbl > MAX_ZLBLS)).any():
            sys.exit('Error: labelled zone with label outside [1,%d]%s'
                     % (MAX_ZLBLS, where))


def validate_curves(nelv, curves, path=''):
    """Refuse malformed curve records (element range, known curve types)."""
    where = (' in %s' % path) if path else ''
    if curves.size == 0:
        return
    if ((curves['e'] < 1) | (curves['e'] > nelv)).any():
        sys.exit('Error: curve record references element outside [1,%d]%s'
                 % (nelv, where))
    ct = curves['type']
    if ((ct != 0) & (ct != 3) & (ct != 4)).any():
        sys.exit('Error: curve record with unknown edge type (valid: 0, '
                 '3 = circle, 4 = midside)%s' % where)


# ---------------------------------------------------------------------------
# re2 read
# ---------------------------------------------------------------------------
class Re2(NamedTuple):
    """A whole .re2 in memory (3D only): raw coordinates + curve/BC records."""
    nelv: int
    version: str
    xyz: np.ndarray       # (nelv, 8, 3) f64 corner coordinates
    curves: np.ndarray    # RE2_CURVE_DT records (as read)
    bcs: np.ndarray       # RE2_BC_DT records (as read)


def read_re2(path, chunk=1 << 21):
    """Read a NEKTON .re2 (versions #v001..#v004, little-endian, 3D)."""
    try:
        f = open(path, 'rb')
    except OSError as ex:
        sys.exit('Error: cannot open %s (%s)' % (path, ex.strerror))
    with f:
        hdr = f.read(80)
        if len(hdr) != 80:
            sys.exit('Error: %s is not a .re2 file (short header)' % path)
        ver = hdr[:5].decode('latin-1')
        try:
            if ver == '#v004':
                nel = int(hdr[5:21]); ndim = int(hdr[21:24])
                nelv = int(hdr[24:40])
            elif ver in ('#v001', '#v002', '#v003'):
                nel = int(hdr[5:14]); ndim = int(hdr[14:17])
                nelv = int(hdr[17:26])
            else:
                sys.exit('Error: unknown re2 version %r' % ver)
        except ValueError:
            sys.exit('Error: cannot parse re2 header of %s' % path)
        v2 = ver != '#v001'
        endian = np.fromfile(f, dtype='<f4', count=1)
        if endian.size != 1 or abs(float(endian[0]) - 6.54321) > 1e-4:
            sys.exit('Error: byte-swapped or corrupt re2 (endian tag)')
        if ndim != 3:
            sys.exit('Error: the light tools support 3D (hex) meshes only')
        del nel

        # elements, chunked (200 B/record dp; v1 is f32 and upcast)
        xyz = np.empty((nelv, 8, 3), dtype=np.float64)
        done = 0
        while done < nelv:
            n = min(chunk, nelv - done)
            rec = np.fromfile(f, dtype=RE2_EL_DT[v2], count=n)
            if rec.size != n:
                sys.exit('Error: truncated or corrupt .re2 file (element '
                         'section, record %d of %d)' % (done + rec.size, nelv))
            xyz[done:done + n, :, 0] = rec['x']
            xyz[done:done + n, :, 1] = rec['y']
            xyz[done:done + n, :, 2] = rec['z']
            done += n
        if not np.isfinite(xyz).all():
            sys.exit('Error: non-finite coordinate in %s' % path)

        ncurve = _re2_count(f, v2, 'curve')
        curves = np.fromfile(f, dtype=RE2_CURVE_DT[v2], count=ncurve)
        if curves.size != ncurve:
            sys.exit('Error: truncated or corrupt .re2 file (curve section)')
        nbc = _re2_count(f, v2, 'boundary-condition')
        bcs = np.fromfile(f, dtype=RE2_BC_DT[v2], count=nbc)
        if bcs.size != nbc:
            sys.exit('Error: truncated or corrupt .re2 file (BC section)')

    # bounds validation up front (validate-and-refuse; nothing written yet)
    ce = curves['e'].astype(np.int64)
    cz = curves['edge'].astype(np.int64)
    if curves.size and (((ce < 1) | (ce > nelv)).any()):
        sys.exit('Error: curve record references element out of range')
    if curves.size and (((cz < 1) | (cz > 12)).any()):
        sys.exit('Error: curve record edge index out of [1,12]')
    be = bcs['e'].astype(np.int64)
    bf = bcs['f'].astype(np.int64)
    if bcs.size and ((be < 1) | (be > nelv)).any():
        sys.exit('Error: BC record references element out of range')
    if bcs.size and ((bf < 1) | (bf > 6)).any():
        sys.exit('Error: BC record face out of [1,6]')
    return Re2(nelv, ver, xyz, curves, bcs)


def _re2_count(f, v2, what):
    if v2:
        a = np.fromfile(f, dtype='<f8', count=1)
    else:
        a = np.fromfile(f, dtype='<i4', count=1)
    if a.size != 1:
        sys.exit('Error: truncated or corrupt .re2 file (missing %s count)'
                 % what)
    n = int(a[0])
    if n < 0:
        sys.exit('Error: negative %s count in .re2' % what)
    return n


def bc_type_str(raw):
    """The BC/curve type field as Neko sees it: leading/trailing blanks and
    NUL padding stripped (Neko reads into a blank-padded CHARACTER)."""
    return raw.decode('latin-1').replace('\x00', ' ').strip()


# ---------------------------------------------------------------------------
# Distribution CSV files (genmeshbox), matching Neko's csv reader semantics
# ---------------------------------------------------------------------------
def read_dist_csv(path, n):
    """Read n+1 grid coordinates from a genmeshbox distribution file.

    Exactly Neko's ``csv_file_read_vector``: if the file has more than one
    line the first is treated as a header and skipped; commas, spaces and
    newlines all act as separators.
    """
    try:
        with open(path, 'r') as f:
            lines = f.read().splitlines()
    except OSError as ex:
        sys.exit('Error: cannot open distribution file %s (%s)'
                 % (path, ex.strerror))
    body = lines[1:] if len(lines) > 1 else lines
    toks = ' '.join(body).replace(',', ' ').split()
    if len(toks) < n + 1:
        sys.exit('Error: expected %d grid values in %s (Neko treats the '
                 'first line as a header when the file has >1 line)'
                 % (n + 1, path))
    try:
        return np.array([float(t) for t in toks[:n + 1]], dtype=np.float64)
    except ValueError:
        sys.exit('Error: cannot parse grid values in %s' % path)


# ---------------------------------------------------------------------------
# fld writing (single precision NEKTON/Neko field file, lx=ly=lz=3)
# ---------------------------------------------------------------------------
def write_zone_indices_fld(out_base, elids, gll_xyz, scal):
    """Write ``<out_base>.fld`` (+ ``.nek5000`` companion): the straight-sided
    trilinear geometry at the 3x3x3 GLL nodes plus one scalar field.

    ``elids`` are the ACTUAL global element ids in record order (a valid
    .nmsh may store its records in any order -- the id list is what maps a
    block to an element), ``gll_xyz`` is (nelv, 27, 3) f64 and ``scal`` is
    (nelv, 27) f32-compatible.
    """
    nel = len(elids)
    hdr = ('#std %1d %2d %2d %2d %10d %10d %20.13E %9d %6d %6d %-10s'
           % (4, 3, 3, 3, nel, nel, 0.0, 1, 1, 1, 'XS01')).ljust(132)
    assert len(hdr) == 132
    gx = gll_xyz.astype(np.float32)
    sc = np.asarray(scal, dtype=np.float32)
    with atomic_output(out_base + '.fld') as f:
        f.write(hdr.encode('ascii'))
        np.float32(6.54321).tofile(f)
        np.asarray(elids, dtype='<i4').tofile(f)
        # geometry block: per element gx(27), gy(27), gz(27)
        np.ascontiguousarray(gx.transpose(0, 2, 1)).tofile(f)
        sc.tofile(f)
        # per-element geometry bounding-box metadata (xmin,xmax,...,zmin,zmax)
        bb = np.empty((nel, 6), dtype=np.float32)
        bb[:, 0::2] = gx.min(axis=1)
        bb[:, 1::2] = gx.max(axis=1)
        bb.tofile(f)
        # per-element (min, max) metadata for the scalar field
        mm = np.stack([sc.min(axis=1), sc.max(axis=1)], axis=1)
        mm.astype(np.float32).tofile(f)
    with atomic_output(out_base + '.nek5000') as f:
        f.write(('filetemplate: %s.fld\nfirsttimestep: 1\nnumtimesteps: 1\n'
                 % os.path.basename(out_base)).encode('ascii'))
