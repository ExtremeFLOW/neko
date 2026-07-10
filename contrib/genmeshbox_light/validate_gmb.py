#!/usr/bin/env python3
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
"""Byte-exactness validation for genmeshbox_light against every genmeshbox-made
box mesh shipped in the Neko tree.

For each candidate `.nmsh` we detect whether it is a uniform tensor-product box
(distinct x/y/z coordinates, uniform spacing, nelv = nx*ny*nz, genmeshbox element
ordering), infer its parameters and per-direction periodicity (from the presence
of type-5 zones on each face pair), regenerate it with genmeshbox_light, and
byte-compare:

  * header + every element record   -> must be byte-identical
  * periodic (type-5) zone records  -> must be byte-identical (incl. merged ids)
  * labeled (type-7) zone records   -> identical except p_e / glb_pt_ids, which
                                       Neko's writer leaves UNINITIALISED (stale
                                       heap); the tool writes deterministic zeros

Shipped box files often carry trailing dead bytes after the valid structure (a
non-truncating writer artifact, exactly like turb_pipe); the harness compares the
structured records by offset and ignores that tail.

Usage:  python3 validate_gmb.py [path/to/genmeshbox_light] [neko_root]
"""
import struct, subprocess, glob, os, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
EXE  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'genmeshbox_light')
NEKO = sys.argv[2] if len(sys.argv) > 2 else os.path.abspath(os.path.join(HERE, '..', '..'))


def read_box(path):
    d = open(path, 'rb').read()
    nelv, gdim = struct.unpack('<ii', d[:8])
    if gdim != 3:
        return None
    off = 8
    elems = []
    xs, ys, zs = set(), set(), set()
    for _ in range(nelv):
        struct.unpack('<i', d[off:off + 4]); off += 4
        vs = []
        for _ in range(8):
            vi = struct.unpack('<i', d[off:off + 4])[0]
            xyz = struct.unpack('<3d', d[off + 4:off + 28]); off += 28
            vs.append((vi, xyz)); xs.add(xyz[0]); ys.add(xyz[1]); zs.add(xyz[2])
        elems.append(vs)
    nz = struct.unpack('<i', d[off:off + 4])[0]; off += 4
    zbase = off
    faces = set()
    for i in range(nz):
        z = struct.unpack('<9i', d[off + i * 36:off + i * 36 + 36])
        if z[8] == 5:
            faces.add(z[1])
    off += nz * 36
    ncurves = struct.unpack('<i', d[off:off + 4])[0]
    return dict(d=d, nelv=nelv, elems=elems, xs=sorted(xs), ys=sorted(ys),
                zs=sorted(zs), nz=nz, zbase=zbase, ncurves=ncurves,
                valid_end=off + 4, per_faces=faces)


def is_uniform(vals):
    if len(vals) < 2:
        return False
    d0 = vals[1] - vals[0]
    if d0 <= 0:
        return False
    return all(abs((vals[k + 1] - vals[k]) - d0) < 1e-9 * max(1.0, abs(d0)) for k in range(len(vals) - 1))


def check(path, tmp):
    b = read_box(path)
    if b is None:
        return None
    nx, ny, nz = len(b['xs']) - 1, len(b['ys']) - 1, len(b['zs']) - 1
    if nx < 1 or ny < 1 or nz < 1 or nx * ny * nz != b['nelv']:
        return None                      # not a tensor-product box
    if b['ncurves'] != 0:
        return None
    # Is this genmeshbox output? Its point ids follow the lexicographic scheme
    # id = 1 + ix + iy*(nx+1) + iz*(nx+1)*(ny+1); other generators (nekbone/tgv/
    # poisson benchmark meshes) number points sequentially per element instead.
    def gid(ix, iy, iz):
        return 1 + ix + iy * (nx + 1) + iz * (nx + 1) * (ny + 1)
    exp1 = [gid(0,0,0), gid(1,0,0), gid(1,1,0), gid(0,1,0),
            gid(0,0,1), gid(1,0,1), gid(1,1,1), gid(0,1,1)]
    if [vi for vi, _ in b['elems'][0]] != exp1:
        return None                      # different generator -> not a genmeshbox box
    px = (1 in b['per_faces'] or 2 in b['per_faces'])
    py = (3 in b['per_faces'] or 4 in b['per_faces'])
    pz = (5 in b['per_faces'] or 6 in b['per_faces'])
    x0, x1 = b['xs'][0], b['xs'][-1]; y0, y1 = b['ys'][0], b['ys'][-1]
    z0, z1 = b['zs'][0], b['zs'][-1]
    outp = os.path.join(tmp, 'g.nmsh')
    # The exact command-line x0/x1 (and rounding order) that made each shipped box
    # can't be recovered from the mesh, so feed the EXACT grid lines as dist files
    # with --direct (verbatim). This validates the ordering / point ids / zone /
    # periodic-merge logic byte-exact, independent of the coordinate-formula detail.
    # Write single-line CSV distribution files (comma-separated, no header) -- the
    # format genmeshbox and the tool read identically (Neko's csv reader skips a
    # header only when the file has more than one line). A separate test below
    # exercises the header-skip path explicitly.
    dx = os.path.join(tmp, 'dx.csv'); dy = os.path.join(tmp, 'dy.csv'); dz = os.path.join(tmp, 'dz.csv')
    open(dx, 'w').write(','.join(repr(v) for v in b['xs']))
    open(dy, 'w').write(','.join(repr(v) for v in b['ys']))
    open(dz, 'w').write(','.join(repr(v) for v in b['zs']))
    args = [EXE, repr(x0), repr(x1), repr(y0), repr(y1), repr(z0), repr(z1),
            str(nx), str(ny), str(nz),
            str(px).lower(), str(py).lower(), str(pz).lower(), dx, dy, dz, outp, '--direct']
    p = subprocess.run(args, capture_output=True, text=True)
    if p.returncode != 0:
        return (path, 'RUN FAIL: ' + (p.stderr.strip()[:80]))
    ref = b['d']; out = open(outp, 'rb').read()
    e1 = 8 + b['nelv'] * 228
    elem_ok = (ref[:e1] == out[:e1])
    # zones by type
    per_ok = lab_only_uninit = True
    z0off = e1 + 4
    for i in range(b['nz']):
        r = struct.unpack('<9i', ref[z0off + i * 36:z0off + i * 36 + 36])
        o = struct.unpack('<9i', out[z0off + i * 36:z0off + i * 36 + 36])
        if r[8] == 5 and r != o:
            per_ok = False
        if r[8] == 7:
            # allowed to differ only in p_e (idx 2) and glb_pt_ids (idx 4..7)
            for k in (0, 1, 3, 8):
                if r[k] != o[k]:
                    lab_only_uninit = False
    ncur_ok = (out[e1 + 4 + b['nz'] * 36:e1 + 8 + b['nz'] * 36] ==
               ref[e1 + 4 + b['nz'] * 36:e1 + 8 + b['nz'] * 36])
    size_ok = (len(out) == b['valid_end'])          # tool writes no trailing tail
    ok = elem_ok and per_ok and lab_only_uninit and ncur_ok and size_ok
    return (path, nx, ny, nz, (px, py, pz), b['nz'], elem_ok, per_ok,
            lab_only_uninit, size_ok, ok)


def csv_header_test(tmp):
    """The distribution file follows Neko's csv reader: a single-line CSV has no
    header, but a file with >1 line has its first line treated as a header and
    skipped. Confirm all three layouts that must agree produce the SAME box:
    single-line CSV, header + one-value-per-line, header + CSV values."""
    xs = [0.0, 0.3, 0.55, 0.8, 1.0]          # nx = 4, non-uniform
    base = None
    layouts = {
        'single-line CSV        ': ','.join(repr(v) for v in xs),
        'header + one-per-line  ': 'xcoord\n' + '\n'.join(repr(v) for v in xs),
        'header + CSV values    ': 'x\n' + ','.join(repr(v) for v in xs),
    }
    allok = True
    print('csv distribution-file header handling (Neko-compatible):')
    for name, text in layouts.items():
        dx = os.path.join(tmp, 'h.csv'); open(dx, 'w').write(text)
        outp = os.path.join(tmp, 'h.nmsh')
        args = [EXE, '0', '1', '0', '1', '0', '1', '4', '1', '1',
                'false', 'false', 'false', dx, 'uniform', 'uniform', outp, '--direct']
        p = subprocess.run(args, capture_output=True, text=True)
        ok = (p.returncode == 0)
        data = open(outp, 'rb').read() if ok else b''
        if base is None and ok:
            base = data
        same = ok and data == base
        allok = allok and same
        print('  %s -> %s' % (name, 'OK (same box)' if same else 'FAIL: ' + p.stderr.strip()[:60]))
    return allok


def main():
    if not os.path.exists(EXE):
        print('error: genmeshbox_light not found at', EXE); return 1
    paths = sorted(glob.glob(os.path.join(NEKO, 'examples', '**', '*.nmsh'), recursive=True) +
                   glob.glob(os.path.join(NEKO, 'tests', '**', '*.nmsh'), recursive=True))
    print('%-46s %10s %8s %8s %6s %5s' %
          ('box mesh', 'grid', 'periodic', 'elem=', 'zones', 'ok'))
    n = nok = 0
    allok = True
    with tempfile.TemporaryDirectory() as tmp:
        for path in paths:
            r = check(path, tmp)
            if r is None:
                continue                 # not a genmeshbox-style box; skip
            if len(r) == 2:              # the tool crashed/errored on a box
                print('%-46s  %s' % (os.path.relpath(r[0], NEKO), r[1]))
                allok = False
                n += 1
                continue
            (p, nx, ny, nz, per, nzc, elem_ok, per_ok, lab_ok, size_ok, ok) = r
            n += 1; nok += ok; allok = allok and ok
            print('%-46s %10s %8s %8s %6d %5s' %
                  (os.path.relpath(p, NEKO), '%dx%dx%d' % (nx, ny, nz),
                   ''.join('xyz'[i] for i in range(3) if per[i]) or '-',
                   elem_ok, nzc, 'OK' if ok else 'FAIL'))
            if not ok:
                print('     elem=%s per=%s lab_only_uninit=%s size=%s' %
                      (elem_ok, per_ok, lab_ok, size_ok))
    print('\ngenmeshbox boxes tested : %d  (others use a different point numbering' % n)
    print('                          -> not genmeshbox output, skipped)')
    print('byte-exact              : %d / %d  (header + every element record and every' % (nok, n))
    print('                          periodic/type-5 zone identical; type-7 labeled zones')
    print('                          differ only in Neko-uninitialised p_e/glb_pt_ids)')
    print()
    with tempfile.TemporaryDirectory() as tmp:
        hdr_ok = csv_header_test(tmp)
    allok = allok and hdr_ok
    print('\nRESULT:', 'PASS' if allok and n > 0 else ('FAIL' if n else 'NO BOXES FOUND'))
    return 0 if (allok and n > 0 and hdr_ok) else 1


if __name__ == '__main__':
    sys.exit(main())
