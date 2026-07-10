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
"""True head-to-head oracle: rea2nbin_light vs a REAL Neko-produced .nmsh.

The round-trip validators synthesize their .re2 from a .nmsh, so they can share a
blind spot with the tool. This one does not: give it a genuine (.re2, .nmsh) pair
where the .nmsh was produced by the real Neko rea2nbin, and it byte-diffs the
tool's output against that reference, section by section, and classifies every
differing byte. The ONLY differences that count as a PASS are in the type-7
(labeled) zone p_e / glb_pt_ids fields, which Neko's writer leaves UNINITIALIZED
(nondeterministic heap) and the tool writes as deterministic zeros -- fields the
nmsh reader never reads. Anything else is a real discrepancy and fails.

Usage:  python3 validate_golden.py <rea2nbin_light> <mesh.re2> <reference.nmsh>
"""
import sys, struct, subprocess, tempfile, os


def sections(d):
    nelv, gdim = struct.unpack('<ii', d[:8])
    e1 = 8 + nelv * 228
    nz = struct.unpack('<i', d[e1:e1 + 4])[0]; z0 = e1 + 4; z1 = z0 + nz * 36
    nc = struct.unpack('<i', d[z1:z1 + 4])[0]; c0 = z1 + 4; c1 = c0 + nc * 532
    return dict(nelv=nelv, gdim=gdim, e0=8, e1=e1, nz=nz, z0=z0, z1=z1,
               nc=nc, c0=c0, c1=c1)


def main():
    if len(sys.argv) != 4:
        print(__doc__); return 2
    exe, re2, ref_path = sys.argv[1:4]
    ref = open(ref_path, 'rb').read()
    with tempfile.TemporaryDirectory() as t:
        out_path = os.path.join(t, 'out.nmsh')
        p = subprocess.run([exe, re2, out_path], capture_output=True, text=True)
        if p.returncode != 0:
            print('tool FAILED:', (p.stderr or p.stdout).strip()[:200]); return 1
        out = open(out_path, 'rb').read()

    print('reference:', os.path.basename(ref_path))
    print('sizes: reference=%d  light=%d  equal=%s' % (len(ref), len(out), len(ref) == len(out)))
    if len(ref) != len(out):
        print('SIZE MISMATCH -> structural difference (curves? zones?). FAIL'); return 1
    s = sections(ref)
    print('nelv=%d gdim=%d nzones=%d ncurves=%d' % (s['nelv'], s['gdim'], s['nz'], s['nc']))

    dc = lambda lo, hi: sum(1 for k in range(lo, hi) if ref[k] != out[k])
    hdr_d = dc(0, 8)
    elem_d = dc(s['e0'] + 0, s['e1'])
    nzint_d = dc(s['e1'], s['z0'])
    curv_d = dc(s['c0'], s['c1'])
    ncint_d = dc(s['z1'], s['c0'])

    # zone diffs, split by type and field
    FIELD = ['e', 'f', 'p_e', 'p_f', 'gp0', 'gp1', 'gp2', 'gp3', 'type']
    benign = real = 0
    per_type_field = {5: [0] * 9, 7: [0] * 9}
    for i in range(s['nz']):
        b = s['z0'] + i * 36
        r = struct.unpack('<9i', ref[b:b + 36]); o = struct.unpack('<9i', out[b:b + 36])
        zt = r[8]
        for f in range(9):
            if r[f] != o[f]:
                per_type_field.setdefault(zt, [0] * 9)[f] += 1
                # benign iff type-7 and field is p_e (2) or glb_pt_ids (4..7)
                if zt == 7 and (f == 2 or 4 <= f <= 7):
                    benign += 4   # 4 bytes per int32 field
                else:
                    real += 4

    print('\n--- differing bytes by section ---')
    print('  header        : %d' % hdr_d)
    print('  elements      : %d / %d' % (elem_d, s['e1'] - s['e0']))
    print('  nzones int    : %d' % nzint_d)
    print('  zones         : %d  (benign type-7 uninit: %d, REAL: %d)' %
          (benign + real, benign, real))
    print('  ncurves int   : %d' % ncint_d)
    print('  curves        : %d / %d' % (curv_d, s['c1'] - s['c0']))
    print('\n  zone field-diff counts (records differing per field):')
    print('    %-6s %s' % ('type', ' '.join('%5s' % f for f in FIELD)))
    for zt in sorted(per_type_field):
        print('    %-6d %s' % (zt, ' '.join('%5d' % c for c in per_type_field[zt])))

    nonzone_real = hdr_d + elem_d + nzint_d + ncint_d + curv_d
    ok = (nonzone_real == 0 and real == 0)
    print('\nRESULT:', 'PASS (byte-exact modulo Neko-uninitialised type-7 zone fields)'
          if ok else 'FAIL -- real discrepancy outside the known type-7 uninit fields')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
