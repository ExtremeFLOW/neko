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
"""Round-trip validation for rea2nbin_light across every 3D .nmsh in the Neko tree.

For each mesh:  nmsh  ->  synthesize a .re2  ->  rea2nbin_light  ->  nmsh'
and check the two invariants that define a correct mesh import:
  * coordinates are bit-exact          (geometry preserved)
  * the point dedup partition matches   (topology: which vertices are shared)

The point-id *values* additionally match exactly for meshes whose reference
numbering follows the re2/rea2nbin first-appearance convention; genmeshbox-
generated box meshes use a different (lexicographic) numbering, so their ids
differ by a clean bijection (a benign relabeling, reported separately).

Usage:  python3 validate_roundtrip.py [path/to/rea2nbin_light] [neko_root]
"""
import struct, subprocess, glob, os, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
EXE  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'rea2nbin_light')
# default: this script lives in contrib/rea2nbin_light/, so neko root is ../../
NEKO = sys.argv[2] if len(sys.argv) > 2 else os.path.abspath(os.path.join(HERE, '..', '..'))


def read_nmsh_elems(path):
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
    return nelv, els


def write_re2(path, nelv, els):
    hdr = '#v002%9d%3d%9d' % (nelv, 3, nelv)
    hdr = hdr + ' ' * (80 - len(hdr))
    with open(path, 'wb') as f:
        f.write(hdr.encode('ascii'))
        f.write(struct.pack('<f', 6.54321))            # endian test
        for _, vs in els:
            f.write(struct.pack('<d', 1.0))            # rgroup (unused)
            for k in range(3):
                f.write(struct.pack('<8d', *[vs[j][1][k] for j in range(8)]))
        f.write(struct.pack('<d', 0.0))                # ncurv = 0
        f.write(struct.pack('<d', 0.0))                # nbcs  = 0


def check(path, tmp):
    r = read_nmsh_elems(path)
    if r is None:
        return None
    nelv, els = r
    # is the reference itself a proper coordinate dedup?
    c2i, i2c, ref_ok = {}, {}, True
    for _, vs in els:
        for vi, xyz in vs:
            if xyz in c2i and c2i[xyz] != vi:
                ref_ok = False
            if vi in i2c and i2c[vi] != xyz:
                ref_ok = False
            c2i[xyz] = vi; i2c[vi] = xyz
    re2 = os.path.join(tmp, 'rt.re2'); out = os.path.join(tmp, 'rt.nmsh')
    write_re2(re2, nelv, els)
    p = subprocess.run([EXE, re2, out], capture_output=True, text=True)
    if p.returncode != 0:
        return (path, nelv, 'RUN FAILED: ' + (p.stderr.strip()[:60] or 'nonzero exit'))
    nelv2, els2 = read_nmsh_elems(out)
    coords_ok = (nelv == nelv2)
    vidx_exact, part_ok = True, True
    lc2i, fmap, gmap = {}, {}, {}
    for e in range(nelv):
        for j in range(8):
            if els[e][1][j][1] != els2[e][1][j][1]:
                coords_ok = False
            vi_ref, xyz = els[e][1][j][0], els2[e][1][j][1]
            vi_out = els2[e][1][j][0]
            if vi_ref != vi_out:
                vidx_exact = False
            if xyz in lc2i and lc2i[xyz] != vi_out:
                part_ok = False
            lc2i[xyz] = vi_out
            # bijection between reference ids and light ids
            if vi_ref in fmap and fmap[vi_ref] != vi_out:
                part_ok = False
            if vi_out in gmap and gmap[vi_out] != vi_ref:
                part_ok = False
            fmap[vi_ref] = vi_out; gmap[vi_out] = vi_ref
    uniq = len(set(v for _, vs in els2 for v, _ in vs))
    return (path, nelv, uniq, ref_ok, coords_ok, vidx_exact, part_ok)


def main():
    if not os.path.exists(EXE):
        print('error: rea2nbin_light not found at', EXE, '(build it first: make)')
        return 1
    paths = sorted(glob.glob(os.path.join(NEKO, 'examples', '**', '*.nmsh'), recursive=True) +
                   glob.glob(os.path.join(NEKO, 'tests', '**', '*.nmsh'), recursive=True))
    print('%-52s %8s %8s  %-5s %-6s %-6s %-6s' %
          ('mesh', 'nelv', 'uniqpts', 'refOK', 'coords', 'vidx=', 'part='))
    n3d = n_coords = n_part = n_vidx = 0
    ok = True
    with tempfile.TemporaryDirectory() as tmp:
        for path in paths:
            r = check(path, tmp)
            if r is None:
                continue
            if len(r) == 3:
                print('%-52s %8d  %s' % (os.path.relpath(r[0], NEKO), r[1], r[2]))
                ok = False; continue
            p, nelv, uniq, refok, cok, vx, pk = r
            n3d += 1; n_coords += cok; n_part += pk; n_vidx += vx
            print('%-52s %8d %8d  %-5s %-6s %-6s %-6s' %
                  (os.path.relpath(p, NEKO), nelv, uniq, refok, cok, vx, pk))
            if not (cok and pk):
                ok = False
    print()
    print('3D meshes tested          : %d' % n3d)
    print('coordinates bit-exact     : %d / %d' % (n_coords, n3d))
    print('topology (partition) exact: %d / %d' % (n_part, n3d))
    print('v_idx values exact        : %d / %d  (rest: benign relabeling on genmeshbox meshes)' % (n_vidx, n3d))
    print('\nRESULT:', 'PASS (geometry + topology preserved on all meshes)' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
