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
"""Curved-element validation for rea2nbin_light across every curved mesh in Neko.

Curve records are what re2 -> nmsh must carry through so the solver can bend
element faces onto circular arcs ('C') or place explicit midside nodes ('m')
in dofmap_generate_xyz. Neko reads each re2 curve record (elem, edge, 5 values,
type char) into curve_data(5,12,elem)/curve_type(12,elem), then writes one nmsh
curve record per curved element (e, curve_data(5,12), type(12)) in ascending
element order. rea2nbin_light must reproduce that exactly.

There is no shipped .re2 with supported ('C'/'m') curves -- hemi's are 's'
(sphere), which Neko itself drops -- so for each curved .nmsh we reconstruct the
equivalent .re2 curve records from the mesh's own nmsh curve section (edge e's
5 values <-> re2 record point(1:5); nmsh type 3->'C', 4->'m'), run the tool, and
compare the resulting nmsh curve section byte-for-byte against the reference.

Note: a few shipped curved .nmsh (turb_pipe, cylinder) carry extra trailing
bytes after their curve section that Neko's reader never reads (it consumes
exactly ncurves records); this harness parses the structured curve records by
offset, so that dead tail is correctly ignored.

Usage:  python3 validate_curved.py [path/to/rea2nbin_light] [neko_root]
"""
import struct, subprocess, glob, os, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
EXE  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'rea2nbin_light')
NEKO = sys.argv[2] if len(sys.argv) > 2 else os.path.abspath(os.path.join(HERE, '..', '..'))

TYPE_CHAR = {3: b'C', 4: b'm'}                          # nmsh curve type -> re2 char


def read_mesh(path):
    d = open(path, 'rb').read()
    nelv, gdim = struct.unpack('<ii', d[:8])
    if gdim != 3:
        return None
    off = 8
    els = []
    for _ in range(nelv):
        el = struct.unpack('<i', d[off:off + 4])[0]; off += 4
        vs = []
        for _ in range(8):
            vi = struct.unpack('<i', d[off:off + 4])[0]
            xyz = struct.unpack('<3d', d[off + 4:off + 28]); off += 28
            vs.append((vi, xyz))
        els.append((el, vs))
    nz = struct.unpack('<i', d[off:off + 4])[0]; off += 4
    off += nz * 36
    nc = struct.unpack('<i', d[off:off + 4])[0]; off += 4
    curves = []                                          # (e, [60 doubles], [12 ints])
    for _ in range(nc):
        e = struct.unpack('<i', d[off:off + 4])[0]
        cd = struct.unpack('<60d', d[off + 4:off + 484])
        ty = struct.unpack('<12i', d[off + 484:off + 532])
        off += 532
        curves.append((e, cd, ty))
    return nelv, els, curves


def write_re2(path, nelv, els, curves):
    hdr = '#v002%9d%3d%9d' % (nelv, 3, nelv); hdr = hdr + ' ' * (80 - len(hdr))
    # explode nmsh per-element curve records into per-edge re2 curve records
    re2curves = []                                       # (elem, edge, [5 dbl], char)
    for (e, cd, ty) in curves:
        for edge in range(12):
            if ty[edge] != 0:
                pts = cd[edge * 5:edge * 5 + 5]          # curve_data(5,12) col-major
                re2curves.append((e, edge + 1, pts, TYPE_CHAR[ty[edge]]))
    with open(path, 'wb') as f:
        f.write(hdr.encode('ascii'))
        f.write(struct.pack('<f', 6.54321))
        for _, vs in els:
            f.write(struct.pack('<d', 1.0))
            for k in range(3):
                f.write(struct.pack('<8d', *[vs[j][1][k] for j in range(8)]))
        f.write(struct.pack('<d', float(len(re2curves))))          # ncurv
        for (elem, edge, pts, ch) in re2curves:                    # v2 curve records
            f.write(struct.pack('<d', float(elem)))
            f.write(struct.pack('<d', float(edge)))
            f.write(struct.pack('<5d', *pts))
            f.write(ch + b' ' * 7)
        f.write(struct.pack('<d', 0.0))                            # nbcs = 0


def check(path, tmp):
    r = read_mesh(path)
    if r is None:
        return None
    nelv, els, curves = r
    if not curves:
        return None
    re2 = os.path.join(tmp, 'c.re2'); out = os.path.join(tmp, 'c.nmsh')
    write_re2(re2, nelv, els, curves)
    p = subprocess.run([EXE, re2, out], capture_output=True, text=True)
    if p.returncode != 0:
        return (path, 'RUN FAIL: ' + p.stderr.strip()[:70])
    _, _, curves2 = read_mesh(out)
    same_n = (len(curves) == len(curves2))
    # element ids + curve_data + type, record by record
    e_ok = cd_ok = ty_ok = same_n
    for a, b in zip(curves, curves2):
        if a[0] != b[0]:
            e_ok = False
        if a[1] != b[1]:
            cd_ok = False
        if a[2] != b[2]:
            ty_ok = False
    types = sorted(set(t for _, _, ty in curves for t in ty if t))
    return (path, len(curves), same_n, e_ok, cd_ok, ty_ok, types)


def main():
    if not os.path.exists(EXE):
        print('error: rea2nbin_light not found at', EXE, '(build it first: make)')
        return 1
    paths = sorted(glob.glob(os.path.join(NEKO, 'examples', '**', '*.nmsh'), recursive=True) +
                   glob.glob(os.path.join(NEKO, 'tests', '**', '*.nmsh'), recursive=True))
    print('%-50s %8s %6s %6s %8s %6s %s' %
          ('curved mesh', 'ncurve', 'count', 'elem', 'crv_data', 'type', 'types'))
    ok = True; n = 0
    with tempfile.TemporaryDirectory() as tmp:
        for path in paths:
            r = check(path, tmp)
            if r is None:
                continue
            if len(r) == 2:
                print('%-50s %s' % (os.path.relpath(r[0], NEKO), r[1])); ok = False; continue
            p, nc, sn, eo, co, to, types = r; n += 1
            print('%-50s %8d %6s %6s %8s %6s %s' %
                  (os.path.relpath(p, NEKO), nc, sn, eo, co, to, types))
            if not (sn and eo and co and to):
                ok = False
    print('\ncurved meshes tested:', n, '  RESULT:',
          'PASS (curve records reproduced bit-exact)' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
