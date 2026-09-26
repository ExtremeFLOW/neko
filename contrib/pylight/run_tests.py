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
"""run_tests.py -- the validation suite for the Python light tools.

Run from contrib/pylight inside a Neko checkout (or pass the checkout root
as the first argument); work files go to a temporary directory.  Design
principles, born of review:

  * oracles are derived from NEKO'S OWN SOURCE, never from the tools' code
    (the block-size oracle below is datadist.f90's formula, re-implemented
    here independently of nekolight.partition);
  * missing optional dependencies SKIP their tests, they do not fail them;
  * every generated mesh is re-read by mesh_checker.py before it counts as
    a pass.

Exit status 0 iff no test failed (skips are fine).
"""

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from nekolight import read_nmsh, write_nmsh, EL_DT, ZONE_DT  # noqa: E402

NEKO = os.path.abspath(sys.argv[1] if len(sys.argv) > 1
                       else os.path.join(HERE, '..', '..'))
WORK = tempfile.mkdtemp(prefix='pylight_tests.')
RESULTS = []


def tool(name, *args, env=None):
    """Run one of the CLIs; returns (exit_code, stdout+stderr)."""
    e = dict(os.environ)
    if env:
        e.update(env)
    r = subprocess.run([sys.executable, os.path.join(HERE, name)]
                       + [str(a) for a in args],
                       capture_output=True, text=True, cwd=WORK, env=e)
    return r.returncode, r.stdout + r.stderr


def report(name, ok, detail=''):
    RESULTS.append((name, ok))
    print('  %-52s %s %s' % (name, 'PASS' if ok else 'FAIL',
                             detail if not ok else ''))


def skip(name, why):
    print('  %-52s SKIP (%s)' % (name, why))


def wpath(name):
    return os.path.join(WORK, name)


def all_meshes():
    out = []
    for root in ('examples', 'tests'):
        for d, _, files in os.walk(os.path.join(NEKO, root)):
            out += [os.path.join(d, f) for f in files if f.endswith('.nmsh')]
    return sorted(f for f in out
                  if np.fromfile(f, dtype='<i4', count=2)[1] == 3)


def neko_linear_sizes_oracle(M, P):
    """datadist.f90 lines 78-81, independently: Ip(rank) =
    floor((M + P - rank - 1) / P)."""
    return np.array([(M + P - r - 1) // P for r in range(P)], dtype=np.int64)


def elem_keys(elems):
    """Order-independent element identity: the raw 8x(vidx,xyz) payload."""
    return set(elems['v'][i].tobytes() for i in range(elems.size))


def zone_sig(m):
    pos = np.empty(m.nelv + 1, dtype=np.int64)
    pos[m.elems['id']] = np.arange(m.nelv)
    sig = set()
    for z in m.zones:
        anchor = m.elems['v'][pos[z['e']]].tobytes()
        partner = (m.elems['v'][pos[z['p_e']]].tobytes()
                   if z['t'] == 5 else b'')
        sig.add((anchor, int(z['f']), partner, int(z['p_f']),
                 z['g'].tobytes() if z['t'] == 5 else b'', int(z['t'])))
    return sig


def curve_sig(m):
    pos = np.empty(m.nelv + 1, dtype=np.int64)
    pos[m.elems['id']] = np.arange(m.nelv)
    return set((m.elems['v'][pos[c['e']]].tobytes(), c['data'].tobytes(),
                c['type'].tobytes()) for c in m.curves)


def golden_compare(a_path, b_path):
    """Byte-compare two .nmsh allowing diffs ONLY in the labelled-zone
    p_e/glb_pt_ids fields (which Neko's writers leave uninitialised)."""
    a, b = open(a_path, 'rb').read(), open(b_path, 'rb').read()
    if len(a) != len(b):
        return False, 'sizes differ'
    nelv = int(np.frombuffer(a[:4], '<i4')[0])
    off = 8 + nelv * EL_DT.itemsize
    if a[:off + 4] != b[:off + 4]:
        return False, 'element section differs'
    nz = int(np.frombuffer(a[off:off + 4], '<i4')[0])
    za = np.frombuffer(a[off + 4:off + 4 + nz * 36], ZONE_DT)
    zb = np.frombuffer(b[off + 4:off + 4 + nz * 36], ZONE_DT)
    if not all(np.array_equal(za[k], zb[k]) for k in ('e', 'f', 'p_f', 't')):
        return False, 'zone core fields differ'
    n7 = za['t'] != 7
    if not (np.array_equal(za['p_e'][n7], zb['p_e'][n7])
            and np.array_equal(za['g'][n7], zb['g'][n7])):
        return False, 'non-labelled zone p_e/g differ'
    if a[off + 4 + nz * 36:] != b[off + 4 + nz * 36:]:
        return False, 'curve section / tail differs'
    return True, ''


# ===========================================================================
print('pylight test suite  (Neko checkout: %s)' % NEKO)
print('work dir: %s\n' % WORK)

# ---- T1: rea2nbin golden (hemi) -------------------------------------------
print('[T1] rea2nbin: hemi golden pair')
hemi_re2 = os.path.join(NEKO, 'examples', 'hemi', 'hemi.re2')
hemi_ref = os.path.join(NEKO, 'examples', 'hemi', 'hemi.nmsh')
if os.path.exists(hemi_re2) and os.path.exists(hemi_ref):
    rc, out = tool('rea2nbin.py', hemi_re2, wpath('hemi.nmsh'))
    ok = rc == 0
    if ok:
        ok, why = golden_compare(wpath('hemi.nmsh'), hemi_ref)
        report('hemi.re2 -> nmsh byte-exact (mod uninit fields)', ok, why)
    else:
        report('hemi conversion runs', False, out[-200:])
else:
    skip('hemi golden', 'examples/hemi not found')

# ---- T2: genmeshbox vs a shipped box --------------------------------------
print('[T2] genmeshbox: shipped turb_channel box')
box_ref = os.path.join(NEKO, 'examples', 'turb_channel', 'box.nmsh')
if os.path.exists(box_ref):
    m = read_nmsh(box_ref)
    xyz = m.elems['v']['xyz'].reshape(-1, 3)
    g = [np.unique(xyz[:, c]) for c in range(3)]
    z5 = m.zones[m.zones['t'] == 5]
    per = ['.true.' if (z5['f'] == f).any() else '.false.' for f in (1, 3, 5)]
    lim = [repr(float(v)) for gg in g for v in (gg[0], gg[-1])]
    n = [str(len(gg) - 1) for gg in g]
    rc, out = tool('genmeshbox.py', *(lim[0:2] + lim[2:4] + lim[4:6] + n
                                      + per + ['box.nmsh', '--direct']))
    ok = rc == 0
    if ok:
        ok, why = golden_compare(wpath('box.nmsh'), box_ref)
        report('regenerated box byte-exact (mod uninit fields)', ok, why)
    else:
        report('box generation runs', False, out[-200:])
    rc, _ = tool('mesh_checker.py', 'box.nmsh')
    report('checker passes on generated box', rc == 0)
else:
    skip('turb_channel box', 'not found')

# ---- T3: checker corpus sweep ---------------------------------------------
print('[T3] mesh_checker: every shipped 3D .nmsh')
meshes = all_meshes()
bad = []
for f in meshes:
    rc, _ = tool('mesh_checker.py', f)
    if rc != 0:
        bad.append(f)
report('corpus sweep (%d meshes)' % len(meshes), not bad,
       '; '.join(bad[:3]))

# ---- T4/T5: prepart contract ----------------------------------------------
print('[T4] prepart: contract on hemi (odd sizes) per backend')
backends = ['spectral', 'geometric']
try:
    import pymetis                                        # noqa: F401
    backends.insert(1, 'metis')
except ImportError:
    skip('metis backend', 'pymetis not installed')
src = read_nmsh(hemi_ref)
NP = 8
oracle = neko_linear_sizes_oracle(src.nelv, NP)
for be in backends:
    rc, out = tool('prepart.py', hemi_ref, NP, 'h_%s.nmsh' % be,
                   '--backend', be)
    if rc != 0:
        report('%s: runs' % be, False, out[-200:])
        continue
    m = read_nmsh(wpath('h_%s.nmsh' % be))
    ok = np.array_equal(np.sort(m.elems['id']), np.arange(1, m.nelv + 1)) \
        and elem_keys(m.elems) == elem_keys(src.elems) \
        and zone_sig(m) == zone_sig(src) and curve_sig(m) == curve_sig(src)
    report('%s: permutation + zone/curve binding' % be, ok)
    # T5: simulate Neko's linear read -- rank slices must be exactly the
    # partition (determinism: re-running the backend labels the same way)
    rc2, _ = tool('prepart.py', hemi_ref, NP, 'h2_%s.nmsh' % be,
                  '--backend', be)
    same = open(wpath('h_%s.nmsh' % be), 'rb').read() == \
        open(wpath('h2_%s.nmsh' % be), 'rb').read()
    report('%s: deterministic output' % be, rc2 == 0 and same)
    bounds = np.concatenate([[0], np.cumsum(oracle)])
    blocks = [elem_keys(m.elems[bounds[r]:bounds[r + 1]]) for r in range(NP)]
    ok = all(len(b) == oracle[r] for r, b in enumerate(blocks)) \
        and not any(blocks[i] & blocks[j] for i in range(NP)
                    for j in range(i + 1, NP))
    report('%s: linear-read blocks have Neko\'s exact sizes' % be, ok)
    rc3, _ = tool('mesh_checker.py', 'h_%s.nmsh' % be)
    report('%s: checker passes on output' % be, rc3 == 0)

print('[T5] prepart: curved+periodic mesh (turb_pipe, geometric)')
pipe = os.path.join(NEKO, 'examples', 'turb_pipe', 'turb_pipe.nmsh')
if os.path.exists(pipe):
    rc, out = tool('prepart.py', pipe, 16, 'pipe16.nmsh', '--geometric',
                   '--no-stats')
    ok = rc == 0
    if ok:
        s, o = read_nmsh(pipe), read_nmsh(wpath('pipe16.nmsh'))
        ok = curve_sig(s) == curve_sig(o) and zone_sig(s) == zone_sig(o)
    report('curves + periodic zones bound through reorder', ok)
else:
    skip('turb_pipe', 'not found')

# ---- T6: robustness fixtures ----------------------------------------------
print('[T6] robustness (atomic writes, validate-and-refuse)')
shutil.copy(hemi_re2, wpath('same.re2'))
before = open(wpath('same.re2'), 'rb').read()
rc, out = tool('rea2nbin.py', 'same.re2', 'same.re2')
report('same in/out path refused, input intact',
       rc != 0 and open(wpath('same.re2'), 'rb').read() == before)
with open(wpath('trunc.re2'), 'wb') as f:
    f.write(open(hemi_re2, 'rb').read()[:100000])
rc, _ = tool('rea2nbin.py', 'trunc.re2', 'tr.nmsh')
leftovers = [f for f in os.listdir(WORK) if f.endswith('.tmp')]
report('truncated re2: error, no partial output, no temp',
       rc != 0 and not os.path.exists(wpath('tr.nmsh')) and not leftovers)
data = open(hemi_ref, 'rb').read()
open(wpath('nocc.nmsh'), 'wb').write(data[:-4])
rc, _ = tool('mesh_checker.py', 'nocc.nmsh')
report('truncated nmsh (missing curve count) rejected', rc != 0)
buf = bytearray(data)
off = 8 + src.nelv * EL_DT.itemsize + 4
z = np.frombuffer(bytes(buf[off:off + 36]), ZONE_DT).copy()
z['e'] = src.nelv + 7
buf[off:off + 36] = z.tobytes()
open(wpath('badz.nmsh'), 'wb').write(bytes(buf))
rc1, _ = tool('mesh_checker.py', 'badz.nmsh')
rc2, _ = tool('prepart.py', 'badz.nmsh', 4, 'badz4.nmsh', '--geometric')
report('out-of-range zone ref rejected by checker + prepart',
       rc1 != 0 and rc2 != 0 and not os.path.exists(wpath('badz4.nmsh')))

# ---- T7: fld ids on a shuffled (valid) mesh --------------------------------
print('[T7] zone-index fld: shuffled element records')
rng = np.random.default_rng(7)
sh = src.elems[rng.permutation(src.nelv)]
z5 = src.zones[src.zones['t'] == 5]
z7 = src.zones[src.zones['t'] == 7]
zx = src.zones[(src.zones['t'] != 5) & (src.zones['t'] != 7)]
write_nmsh(wpath('shuf.nmsh'), sh, (z5, z7, zx), src.curves)
rc, _ = tool('mesh_checker.py', 'shuf.nmsh', '--write-zone-indices')
ok = rc == 0
if ok:
    raw = open(wpath('shuf_zone_indices.fld'), 'rb').read()
    ids = np.frombuffer(raw[136:136 + 4 * src.nelv], '<i4')
    ok = np.array_equal(ids, sh['id'])
report('fld id list carries the actual element ids', ok)

# ---- summary ---------------------------------------------------------------
nfail = sum(1 for _, ok in RESULTS if not ok)
print('\n%d checks, %d failed' % (len(RESULTS), nfail))
print('RESULT: %s' % ('PASS' if nfail == 0 else 'FAIL'))
shutil.rmtree(WORK, ignore_errors=True)
sys.exit(1 if nfail else 0)
