#!/usr/bin/env python3
# Copyright (c) 2026, The Neko Authors
# All rights reserved.  BSD-3-Clause: see nekolight/formats.py for the full
# licence text.
#
#     _  __  ____  __ __  ____
#    / |/ / / __/ / //_/ / __ \
#   /    / / _/  / ,<   / /_/ /
#  /_/|_/ /___/ /_/|_|  \____/
#
"""genmeshbox -- generate a box mesh as a Neko .nmsh.

Every point id, element id and coordinate is a closed-form function of the
grid indices, evaluated with numpy broadcasting and streamed to the file in
fixed-size chunks, so memory is O(chunk + nelx + nely + nelz).  The output
is byte-identical to Neko's contrib/genmeshbox (except the labelled-zone
p_e/glb_pt_ids fields, which genmeshbox leaves uninitialised and this tool
writes as deterministic zeros).

Usage (matches genmeshbox):
  genmeshbox.py x0 x1 y0 y1 z0 z1 nelx nely nelz \\
      [px py pz (.true./.false.)] [dist_x dist_y dist_z] [out.nmsh] [--direct]

px/py/pz make that direction periodic.  dist_* is 'uniform' or a text file
of nel+1 grid coordinates read exactly as Neko's csv reader does (first line
skipped as a header when the file has more than one line).  --direct uses
the closed-form x0 + k*(x1-x0)/n coordinates (pre-2024 genmeshbox; what the
shipped example boxes store) instead of the current cumulative-sum ones.
"""

import sys

import numpy as np

from nekolight import (banner, EL_DT, ZONE_DT, SIJK, FACE_RE2, atomic_output,
                       read_dist_csv)

CHUNK = 1 << 20
OPP = np.array([2, 1, 4, 3, 6, 5], dtype=np.int64)   # facet -> opposite facet


def log(msg):
    print(msg, flush=True)


def parse_args(argv):
    pos, direct = [], False
    for a in argv:
        if a == '--direct':
            direct = True
        else:
            pos.append(a)
    na = len(pos)
    if na < 9:
        sys.exit(__doc__)
    if na not in (9, 10, 12, 13, 15, 16):
        sys.exit('Error: invalid number of positional arguments (%d); valid: '
                 '9, 12 or 15, optionally followed by an output filename '
                 '(10, 13, 16)' % na)
    try:
        box = [float(v) for v in pos[:6]]
        n = [int(v) for v in pos[6:9]]
    except ValueError as ex:
        sys.exit('Error: cannot parse box arguments (%s)' % ex)
    per = [False, False, False]
    dist = ['uniform'] * 3
    out = 'box.nmsh'
    if na >= 12:
        for k in range(3):
            s = pos[9 + k].lower()
            if s not in ('.true.', '.false.', 'true', 'false', 't', 'f'):
                sys.exit('Error: cannot parse "%s" as a logical '
                         '(.true./.false.)' % pos[9 + k])
            per[k] = s in ('.true.', 'true', 't')
    if na >= 15:
        dist = pos[12:15]
    if na in (10, 13, 16):
        out = pos[na - 1]
    return box, n, per, dist, out, direct


def build_axis(n, a0, a1, fname, direct):
    """Grid-line coordinates (n+1 values).  Default reproduces the CURRENT
    genmeshbox exactly: element lengths first, then cumulative summation;
    --direct is the closed form / verbatim file values (the pre-2024
    behaviour, stored by the shipped example boxes)."""
    fv = None
    if fname != 'uniform':
        fv = read_dist_csv(fname, n)
    if direct:
        if fv is None:
            k = np.arange(n + 1, dtype=np.float64)
            return a0 + ((a1 - a0) / float(n)) * k
        return fv
    ellen = (np.full(n, (a1 - a0) / float(n)) if fv is None
             else np.diff(fv))
    # cumulative sum SEEDED with a0 -- the identical floating-point chain to
    # genmeshbox's g(k+1) = g(k) + ellen(k), which a0 + cumsum(ellen) is not
    return np.cumsum(np.concatenate(([a0], ellen)))


def main():
    box, (nx, ny, nz), per, dist, out, direct = parse_args(sys.argv[1:])
    x0, x1, y0, y1, z0, z1 = box
    px, py, pz = per
    if nx < 1 or ny < 1 or nz < 1:
        sys.exit('Error: nelx/nely/nelz must all be >= 1 (got %d %d %d)'
                 % (nx, ny, nz))
    nel = nx * ny * nz
    npts = (nx + 1) * (ny + 1) * (nz + 1)
    if nel > np.iinfo(np.int32).max or npts > np.iinfo(np.int32).max:
        sys.exit('Error: mesh too large -- element/point counts must fit '
                 '32-bit .nmsh ids')

    log(banner('genmeshbox'))
    gx = build_axis(nx, x0, x1, dist[0], direct)
    gy = build_axis(ny, y0, y1, dist[1], direct)
    gz = build_axis(nz, z0, z1, dist[2], direct)

    def pt_id(a, b, c):
        """Raw lexicographic 1-based point id of grid node (a, b, c)."""
        return (1 + a + b * (nx + 1) + c * (nx + 1) * (ny + 1)).astype(np.int32)

    def merged_id(a, b, c):
        """Periodic-merged id: wrap the high boundary of each periodic
        direction to 0 (min-id merge), matching create_periodic_ids."""
        a = np.where(a == nx, 0, a) if px else a
        b = np.where(b == ny, 0, b) if py else b
        c = np.where(c == nz, 0, c) if pz else c
        return pt_id(a, b, c)

    def el_id(ex, ey, ez):
        return (1 + ex + ey * nx + ez * nx * ny).astype(np.int32)

    log('  [1/2] writing %d elements ...' % nel)
    with atomic_output(out) as f:
        np.array([nel, 3], dtype='<i4').tofile(f)
        for s in range(0, nel, CHUNK):
            idx = np.arange(s, min(s + CHUNK, nel), dtype=np.int64)
            ez, ey, ex = idx // (nx * ny), (idx // nx) % ny, idx % nx
            rec = np.empty(idx.size, dtype=EL_DT)
            rec['id'] = (idx + 1).astype(np.int32)
            for sl in range(8):
                a, b, c = ex + SIJK[sl, 0], ey + SIJK[sl, 1], ez + SIJK[sl, 2]
                rec['v']['idx'][:, sl] = pt_id(a, b, c)
                rec['v']['xyz'][:, sl, 0] = gx[a]
                rec['v']['xyz'][:, sl, 1] = gy[b]
                rec['v']['xyz'][:, sl, 2] = gz[c]
            rec.tofile(f)

        # ---- zones (boundary-sized): periodic first (x, y, z), then
        #      labelled by ascending label; loop order matches genmeshbox ----
        log('  [2/2] writing zones ...')

        def face_elems(face):
            """(ex, ey, ez) of the elements on box face 1..6, in genmeshbox
            marking order (outer/inner loops per direction)."""
            if face <= 2:            # x: e_y outer, e_z inner
                a = np.repeat(np.arange(ny), nz)
                b = np.tile(np.arange(nz), ny)
                ex = np.full(a.size, 0 if face == 1 else nx - 1)
                return ex, a, b
            if face <= 4:            # y: e_x outer, e_z inner
                a = np.repeat(np.arange(nx), nz)
                b = np.tile(np.arange(nz), nx)
                ey = np.full(a.size, 0 if face == 3 else ny - 1)
                return a, ey, b
            a = np.repeat(np.arange(nx), ny)     # z: e_x outer, e_y inner
            b = np.tile(np.arange(ny), nx)
            ez = np.full(a.size, 0 if face == 5 else nz - 1)
            return a, b, ez

        def periodic_records(face):
            ex, ey, ez = face_elems(face)
            p_ex, p_ey, p_ez = ex.copy(), ey.copy(), ez.copy()
            if face <= 2:
                p_ex[:] = nx - 1 if face == 1 else 0
            elif face <= 4:
                p_ey[:] = ny - 1 if face == 3 else 0
            else:
                p_ez[:] = nz - 1 if face == 5 else 0
            z = np.zeros(ex.size, dtype=ZONE_DT)
            z['e'] = el_id(ex, ey, ez)
            z['f'] = face
            z['p_e'] = el_id(p_ex, p_ey, p_ez)
            z['p_f'] = OPP[face - 1]
            for k in range(4):
                sl = FACE_RE2[face - 1, k]
                z['g'][:, k] = merged_id(ex + SIJK[sl, 0], ey + SIJK[sl, 1],
                                         ez + SIJK[sl, 2])
            z['t'] = 5
            return z

        def labelled_records(face):
            ex, ey, ez = face_elems(face)
            z = np.zeros(ex.size, dtype=ZONE_DT)
            z['e'] = el_id(ex, ey, ez)
            z['f'] = face
            z['p_f'] = face          # the label; p_e/glb_pt_ids stay 0
            z['t'] = 7
            return z

        zparts = []
        for d, flag in ((0, px), (1, py), (2, pz)):
            if flag:
                zparts.append(periodic_records(2 * d + 1))
                zparts.append(periodic_records(2 * d + 2))
        for lbl in range(1, 7):
            d = (lbl - 1) // 2
            if not per[d]:
                zparts.append(labelled_records(lbl))

        nzones = sum(z.shape[0] for z in zparts)
        np.array([nzones], dtype='<i4').tofile(f)
        for z in zparts:
            z.tofile(f)
        np.array([0], dtype='<i4').tofile(f)             # ncurves

    nper = sum(z.shape[0] for z in zparts if z['t'][0] == 5) if zparts else 0
    log('genmeshbox: %dx%dx%d = %d' % (nx, ny, nz, nel))
    log('Wrote %d periodic + %d labelled zone facets'
        % (nper, nzones - nper))
    log('Done. Wrote %s' % out)


if __name__ == '__main__':
    main()
