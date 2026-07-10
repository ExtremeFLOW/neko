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
"""Periodic-BC validation for rea2nbin_light across every periodic 3D mesh in Neko.

Neko encodes periodicity in two places in a .nmsh:
  * the element records use *un-merged* dedup point ids (via reset_periodic_ids),
  * each periodic (type-5) zone carries the *merged* point ids in glb_pt_ids,
    which is what tells the solver "these two facets are the same nodes".

rea2nbin_light must reproduce that merge from the raw re2 'P' BC records alone
(a 3-pass coordinate-matched union-find, min-id as root -- exactly Neko's
mesh_create_periodic_ids). There is no shipped .re2 with periodic BCs, so for
each periodic .nmsh we synthesize the equivalent .re2 ('P' records rebuilt from
the periodic zones via the neko-facet->re2-face map), run the tool, and check:

  struct     : the periodic zone records (e, f, p_e, p_f) match position-by-position
  partition  : the *merge partition* of physical coordinates is identical -- i.e.
               grouping coordinates by their merged id yields the same set of
               groups as the reference (numbering-independent: robust to any
               relabeling of the merged ids)
  consist    : every physical coordinate maps to exactly one merged id in both
               the reference and the tool output (no split/contradictory merges)

Usage:  python3 validate_periodic.py [path/to/rea2nbin_light] [neko_root]
"""
import struct, subprocess, glob, os, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
EXE  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'rea2nbin_light')
# default: this script lives in contrib/rea2nbin_light/, so neko root is ../../
NEKO = sys.argv[2] if len(sys.argv) > 2 else os.path.abspath(os.path.join(HERE, '..', '..'))

INV = {1: 4, 2: 2, 3: 1, 4: 3, 5: 5, 6: 6}                 # neko facet -> re2 face
FACE_RE2 = {1: [1, 5, 8, 4], 2: [2, 6, 7, 3], 3: [1, 2, 6, 5],
            4: [4, 3, 7, 8], 5: [1, 2, 3, 4], 6: [5, 6, 7, 8]}


def read(path):
    d = open(path, 'rb').read()
    nelv, gdim = struct.unpack('<ii', d[:8])
    if gdim != 3:
        return None
    off = 8; els = []
    for _ in range(nelv):
        el = struct.unpack('<i', d[off:off + 4])[0]; off += 4; vs = []
        for _ in range(8):
            vi = struct.unpack('<i', d[off:off + 4])[0]
            xyz = struct.unpack('<3d', d[off + 4:off + 28]); off += 28
            vs.append((vi, xyz))
        els.append((el, vs))
    nz = struct.unpack('<i', d[off:off + 4])[0]; off += 4
    Z = [struct.unpack('<9i', d[off + i * 36:off + i * 36 + 36]) for i in range(nz)]
    return nelv, els, Z


def write_re2_periodic(path, nelv, els, per):
    hdr = '#v002%9d%3d%9d' % (nelv, 3, nelv); hdr = hdr + ' ' * (80 - len(hdr))
    with open(path, 'wb') as f:
        f.write(hdr.encode()); f.write(struct.pack('<f', 6.54321))
        for _, vs in els:
            f.write(struct.pack('<d', 1.0))
            for k in range(3):
                f.write(struct.pack('<8d', *[vs[j][1][k] for j in range(8)]))
        f.write(struct.pack('<d', 0.0))                    # ncurv = 0
        f.write(struct.pack('<d', float(len(per))))        # nbcs
        for (e, ff, pe, pf) in per:                        # 'P' records
            ty = b'P' + b' ' * 7
            f.write(struct.pack('<7d', float(e), float(INV[ff]),
                                float(pe), float(INV[pf]), 0.0, 0.0, 0.0))
            f.write(ty)


def partition(els, Z):
    # coord -> merged_id from periodic (type-5) zones; count contradictory assigns
    c2m = {}; multi = 0
    for z in Z:
        if z[8] != 5:
            continue
        e, ff = z[0], z[1]; ids = z[4:8]
        _, vs = els[e - 1]
        for k in range(4):
            coord = vs[FACE_RE2[ff][k] - 1][1]
            if coord in c2m and c2m[coord] != ids[k]:
                multi += 1
            c2m[coord] = ids[k]
    groups = {}
    for coord, m in c2m.items():
        groups.setdefault(m, set()).add(coord)
    return frozenset(frozenset(g) for g in groups.values()), multi, len(c2m)


def check(path, tmp):
    r = read(path)
    if r is None:
        return None
    nelv, els, Z = r
    per = [(z[0], z[1], z[2], z[3]) for z in Z if z[8] == 5]
    if not per:
        return None
    re2 = os.path.join(tmp, 'p.re2'); out = os.path.join(tmp, 'p.nmsh')
    write_re2_periodic(re2, nelv, els, per)
    p = subprocess.run([EXE, re2, out], capture_output=True, text=True)
    if p.returncode != 0:
        return (path, 'RUN FAIL: ' + p.stderr.strip()[:70])
    _, els2, Z2 = read(out)
    per2 = [(z[0], z[1], z[2], z[3]) for z in Z2 if z[8] == 5]
    struct_ok = (per == per2)
    part_o, mo, no = partition(els, Z); part_m, mm, nm = partition(els2, Z2)
    return (path, len(per), struct_ok, (part_o == part_m), mo, mm, no)


def main():
    if not os.path.exists(EXE):
        print('error: rea2nbin_light not found at', EXE, '(build it first: make)')
        return 1
    paths = sorted(glob.glob(os.path.join(NEKO, 'examples', '**', '*.nmsh'), recursive=True) +
                   glob.glob(os.path.join(NEKO, 'tests', '**', '*.nmsh'), recursive=True))
    print('%-48s %8s %8s %10s %8s' %
          ('periodic mesh', '#per', 'struct', 'partition', 'consist'))
    ok = True; n = 0
    with tempfile.TemporaryDirectory() as tmp:
        for path in paths:
            r = check(path, tmp)
            if r is None:
                continue
            if len(r) == 2:
                print('%-48s %s' % (os.path.relpath(r[0], NEKO), r[1])); ok = False; continue
            p, npf, sk, pk, mo, mm, no = r; n += 1
            cons = (mo == 0 and mm == 0)
            print('%-48s %8d %8s %10s %8s' % (os.path.relpath(p, NEKO), npf, sk, pk, cons))
            if not (sk and pk and cons):
                ok = False
    print('\nperiodic meshes tested:', n, '  RESULT:',
          'PASS (periodic identification matches Neko)' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
