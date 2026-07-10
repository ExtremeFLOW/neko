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
"""Validation for prepart_light across the 3D meshes in a Neko tree.

For each mesh x backend it runs prepart_light.py and checks, independently of
the tool's internals:

  1. PURE PERMUTATION -- the output's element payloads (vertex ids + coords) are
     exactly a permutation of the input's; el_idx is contiguous 1..nelv.
  2. CURVES CARRIED + BOUND -- every curve payload survives verbatim and is
     attached to the *geometrically identical* element (matched by the element's
     8 (id,x,y,z) tuples), not merely preserved as a multiset.
  3. ZONES CARRIED + BOUND -- every periodic (type-5) zone keeps f/p_f/glb_pt_ids
     and lands on the geometrically identical element; labeled (type-7) zones
     keep their label and binding.
  4. CONTIGUOUS BLOCKS -- elements appear grouped by partition (derived by
     re-running the backend in-process), i.e. a linear read reproduces the
     partition (the prepart contract).
  5. BALANCE + CUT -- part sizes differ by <= 1 (spectral/geometric; METIS is
     allowed its ~3% imbalance tolerance) and the weighted edge cut is reported
     against a naive linear split.

Usage:  python3 validate_prepart.py [neko_root] [nparts]
"""
import collections
import glob
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prepart_light as pl                                     # noqa: E402

NEKO = sys.argv[1] if len(sys.argv) > 1 else os.path.abspath(
    os.path.join(HERE, '..', '..'))
NPARTS = int(sys.argv[2]) if len(sys.argv) > 2 else 8
BACKENDS = (sys.argv[3].split(',') if len(sys.argv) > 3
            else ['spectral', 'metis', 'geometric'])
SPECTRAL_MAX = 40000        # keep the (slowest) spectral backend tractable
ANY_MAX = 300000


def el_keys(elems):
    """Hashable per-element geometry keys (the 228-byte payload minus el_idx)."""
    return [e['v'].tobytes() for e in elems]


def check_mesh(path, backend, nparts):
    nelv, elems, zones, curves = pl.read_nmsh(path)
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, 'out.nmsh')
        r = subprocess.run([sys.executable,
                            os.path.join(HERE, 'prepart_light.py'),
                            path, str(nparts), out, '--backend', backend],
                           capture_output=True, text=True)
        if r.returncode != 0:
            return dict(fail='RUN: ' + r.stderr.strip()[:100])
        nelv2, elems2, zones2, curves2 = pl.read_nmsh(out)

    res = {}
    res['contig'] = (nelv2 == nelv and
                     np.array_equal(elems2['id'], np.arange(1, nelv + 1)))
    ki, ko = el_keys(elems), el_keys(elems2)
    res['perm'] = collections.Counter(ki) == collections.Counter(ko)

    # geometry-binding maps: element key -> positions
    pos_in = collections.defaultdict(list)
    for p, k in enumerate(ki):
        pos_in[k].append(p)
    pos_out = collections.defaultdict(list)
    for p, k in enumerate(ko):
        pos_out[k].append(p)

    pos_of_elid_in = np.full(nelv + 1, -1, dtype=np.int64)
    pos_of_elid_in[elems['id'].astype(np.int64)] = np.arange(nelv)

    # curves: payload verbatim + bound to the geometrically identical element
    ok = True
    for c_in in curves:
        p_in = pos_of_elid_in[int(c_in['e'])]
        k = ki[p_in]
        match = False
        for c_out in curves2:
            if (c_out['data'].tobytes() == c_in['data'].tobytes() and
                    c_out['type'].tobytes() == c_in['type'].tobytes()):
                p_out = int(c_out['e']) - 1              # output ids contiguous
                if ko[p_out] == k:
                    match = True
                    break
        if not match:
            ok = False
            break
    res['curves'] = ok and (curves2.shape[0] == curves.shape[0])

    # zones: binding by (element geometry key, f) with payload carried
    def zone_sig(zs, keys, pos_of, out_side):
        sig = collections.Counter()
        for z in zs:
            e = int(z['e'])
            if not (1 <= e <= nelv):
                sig[('oor', bytes(z.tobytes()))] += 1   # passed through verbatim
                continue
            p = (e - 1) if out_side else int(pos_of[e])
            if z['t'] == 5:
                sig[(5, keys[p], int(z['f']), int(z['p_f']),
                     z['g'].tobytes())] += 1
            elif z['t'] == 7:
                sig[(7, keys[p], int(z['f']), int(z['p_f']))] += 1
            else:
                # legacy types: whole payload (minus e) verbatim + binding
                sig[(int(z['t']), keys[p], int(z['f']), int(z['p_e']),
                     int(z['p_f']), z['g'].tobytes())] += 1
        return sig
    si = zone_sig(zones, ki, pos_of_elid_in, False)
    so = zone_sig(zones2, ko, None, True)
    res['zones'] = (si == so)

    # contiguous blocks: recompute the partition in-process and compare
    cell, pos_of_elid = pl.merged_cell(nelv, elems, zones)
    A = pl.dual_graph(cell)
    if backend == 'spectral':
        part = pl.spectral_partition(A, nelv, nparts)
    elif backend == 'metis':
        part = pl.metis_partition(A, nelv, nparts)
    else:
        part = pl.geometric_partition(elems, nelv, nparts)
    order = np.argsort(part, kind='stable')
    expect = el_keys(elems[order])
    res['blocks'] = (expect == ko)

    sizes = np.bincount(part, minlength=nparts)
    if backend == 'metis':
        # METIS balances the LARGEST part against the ideal size within its
        # ufactor tolerance (~3%); with tiny parts integer granularity adds +-1.
        # max/min is the wrong metric there (a [3,...,3,4] split of 25 into 8 is
        # optimal yet max/min = 1.33). What matters for solver load is the max.
        ideal = nelv / nparts
        res['balance'] = sizes.max() <= np.ceil(1.03 * ideal) + 1
    else:
        res['balance'] = (sizes.max() - sizes.min()) <= 1
    cut, total = pl.weighted_cut(A, part)
    lin = np.repeat(np.arange(nparts),
                    np.diff(np.linspace(0, nelv, nparts + 1).astype(int)))
    cut_lin, _ = pl.weighted_cut(A, lin[:nelv])
    res['cut'] = cut
    res['cut_lin'] = cut_lin
    res['nelv'] = nelv
    return res


def main():
    paths = sorted(
        glob.glob(os.path.join(NEKO, 'examples', '**', '*.nmsh'),
                  recursive=True) +
        glob.glob(os.path.join(NEKO, 'tests', '**', '*.nmsh'), recursive=True))
    print('prepart_light validation  (nparts = %d)\n' % NPARTS)
    print('%-42s %-9s %7s %9s %9s %5s' %
          ('mesh', 'backend', 'nelv', 'cut', 'cut_lin', 'ok'))
    allok = True
    n = 0
    for path in paths:
        try:
            nelv, _, _, _ = pl.read_nmsh(path)
        except SystemExit:
            continue
        if nelv < NPARTS or nelv > ANY_MAX:
            continue
        for backend in BACKENDS:
            if backend == 'spectral' and nelv > SPECTRAL_MAX:
                continue
            r = check_mesh(path, backend, NPARTS)
            n += 1
            if 'fail' in r:
                print('%-42s %-9s  FAIL %s' %
                      (os.path.relpath(path, NEKO), backend, r['fail']))
                allok = False
                continue
            checks = ('contig', 'perm', 'curves', 'zones', 'blocks', 'balance')
            ok = all(r[c] for c in checks)
            allok = allok and ok
            note = '' if ok else str({c: r[c] for c in checks if not r[c]})
            print('%-42s %-9s %7d %9d %9d %5s %s' %
                  (os.path.relpath(path, NEKO), backend, r['nelv'],
                   r['cut'], r['cut_lin'], 'OK' if ok else 'FAIL', note))
    print('\nmesh x backend combinations tested: %d' % n)
    print('checks: permutation + contiguous el_idx + curve/zone payloads bound')
    print('        to the geometrically identical elements + contiguous blocks')
    print('        + balance.')
    print('\nRESULT:', 'PASS' if allok and n > 0 else
          ('FAIL' if n else 'NO MESHES FOUND'))
    return 0 if (allok and n > 0) else 1


if __name__ == '__main__':
    sys.exit(main())
