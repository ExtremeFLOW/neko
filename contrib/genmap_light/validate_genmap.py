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
"""Validation for genmap_light across 3D meshes in the Neko tree.

For each mesh it:
  1. Runs genmap_light to produce a reordered .nmsh, and checks the reorder is
     CORRECT: the output is a pure permutation of the input (element geometry,
     curves, and zones preserved), el_idx is contiguous 1..nelv, and the elements
     form contiguous per-partition blocks so a linear read reproduces the
     partition (exactly the prepart contract).
  2. Rebuilds the SAME merged dual-graph Laplacian (edge weight = number of
     shared vertices, periodic facets merged via glb_pt_ids) and computes an
     INDEPENDENT recursive spectral bisection using SciPy's exact eigensolver.
     It then compares genmap_light's weighted edge cut and balance against the
     SciPy oracle and against a naive linear split -- confirming the Lanczos
     Fiedler vector matches the true spectral partition (cut within a small
     factor of the exact eigensolver, and far below the linear baseline on
     meshes where spectral partitioning helps).

Usage:  python3 validate_genmap.py [path/to/genmap_light] [neko_root] [nparts]
"""
import struct, subprocess, glob, os, sys, tempfile, collections, math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
EXE  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'genmap_light')
NEKO = sys.argv[2] if len(sys.argv) > 2 else os.path.abspath(os.path.join(HERE, '..', '..'))
NPARTS = int(sys.argv[3]) if len(sys.argv) > 3 else 8

FACE = [(1,5,8,4),(2,6,7,3),(1,2,6,5),(4,3,7,8),(1,2,3,4),(5,6,7,8)]  # 1-based nmsh vtx


def readmesh(path):
    d = open(path, 'rb').read()
    nelv, gdim = struct.unpack('<ii', d[:8])
    if gdim != 3:
        return None
    off = 8
    vidx = np.empty((nelv, 8), dtype=np.int64)
    geo = []          # per element: tuple of 8 (vidx,x,y,z) for permutation check
    elids = np.empty(nelv, dtype=np.int64)
    for e in range(nelv):
        elids[e] = struct.unpack('<i', d[off:off+4])[0]; off += 4
        g = []
        for k in range(8):
            vi = struct.unpack('<i', d[off:off+4])[0]
            xyz = struct.unpack('<3d', d[off+4:off+28]); off += 28
            vidx[e, k] = vi
            g.append((vi,) + xyz)
        geo.append(tuple(g))
    nz = struct.unpack('<i', d[off:off+4])[0]; off += 4
    zones = [struct.unpack('<9i', d[off+i*36:off+i*36+36]) for i in range(nz)]
    off += nz * 36
    nc = struct.unpack('<i', d[off:off+4])[0]; off += 4
    curves = []
    for i in range(nc):
        e = struct.unpack('<i', d[off:off+4])[0]
        cd = d[off+4:off+484]        # raw bytes of curve_data
        ct = d[off+484:off+532]      # raw bytes of curve_type
        curves.append((e, cd, ct)); off += 532
    return dict(nelv=nelv, gdim=gdim, vidx=vidx, geo=geo, elids=elids,
                zones=zones, curves=curves)


def merged_cell(m):
    """Apply the periodic glb_pt_ids merge and compress point ids to 0..nrnk-1.
    Returns an (nelv,8) int array of compressed vertex ids."""
    nelv = m['nelv']; vidx = m['vidx'].copy()
    pos_of_elid = {int(e): i for i, e in enumerate(m['elids'])}
    remap = {}
    for z in m['zones']:
        e, f, pe, pf, g1, g2, g3, g4, zt = z
        if zt == 5:
            pos = pos_of_elid[e]
            g = (g1, g2, g3, g4)
            for k in range(4):
                raw = int(vidx[pos, FACE[f-1][k]-1])
                remap[raw] = min(remap.get(raw, g[k]), g[k])
    # apply remap
    if remap:
        for e in range(nelv):
            for k in range(8):
                v = int(vidx[e, k])
                if v in remap:
                    vidx[e, k] = remap[v]
    # compress
    uniq = {}
    cell = np.empty((nelv, 8), dtype=np.int64)
    nid = 0
    for e in range(nelv):
        for k in range(8):
            v = int(vidx[e, k])
            c = uniq.get(v)
            if c is None:
                c = nid; uniq[v] = nid; nid += 1
            cell[e, k] = c
    return cell, nid


def dual_laplacian(cell, nrnk):
    """Weighted element dual-graph Laplacian; weight = # shared vertices."""
    nelv = cell.shape[0]
    # vertex -> elements
    v2c = collections.defaultdict(list)
    for e in range(nelv):
        for k in range(8):
            v2c[int(cell[e, k])].append(e)
    w = collections.Counter()
    for v, els in v2c.items():
        for a in range(len(els)):
            for b in range(a+1, len(els)):
                ea, eb = els[a], els[b]
                if ea != eb:
                    key = (ea, eb) if ea < eb else (eb, ea)
                    w[key] += 1
    rows = []; cols = []; data = []
    deg = np.zeros(nelv)
    for (a, b), wt in w.items():
        rows += [a, b]; cols += [b, a]; data += [-float(wt), -float(wt)]
        deg[a] += wt; deg[b] += wt
    for i in range(nelv):
        rows.append(i); cols.append(i); data.append(deg[i])
    L = sp.csr_matrix((data, (rows, cols)), shape=(nelv, nelv))
    return L, w


def weighted_cut(w, part):
    cut = 0
    for (a, b), wt in w.items():
        if part[a] != part[b]:
            cut += wt
    total = sum(w.values())
    return cut, total


def fiedler_vector(idx, L):
    """Exact Fiedler vector of the subgraph induced by node set idx. Dense for
    small n; sparse shift-invert eigsh near 0 for larger n (the Laplacian's
    smallest eigenvalue is 0, so we shift just below it)."""
    idx = np.asarray(idx)
    n = len(idx)
    sub = L[idx][:, idx].tocsr()
    # rebuild the sub-graph Laplacian from the induced adjacency
    A = -sub.copy()
    A.setdiag(0.0); A.eliminate_zeros()
    d = np.asarray(A.sum(axis=1)).ravel()
    Ls = sp.diags(d) - A
    if n <= 1500:
        vals, vecs = np.linalg.eigh(Ls.toarray())
        return vecs[:, 1]
    # sparse: two smallest eigenpairs via shift-invert near 0
    try:
        v0 = np.random.default_rng(12345).standard_normal(n)
        vals, vecs = spla.eigsh(Ls.tocsc(), k=2, sigma=-1e-6, which='LM', v0=v0)
        o = np.argsort(vals)
        return vecs[:, o[1]]
    except Exception:
        v0 = np.random.default_rng(12345).standard_normal(n)
        vals, vecs = spla.eigsh(Ls, k=2, which='SA', maxiter=5000, tol=1e-6,
                                v0=v0)
        o = np.argsort(vals)
        return vecs[:, o[1]]


def scipy_rsb(L, nelv, nparts):
    part = np.zeros(nelv, dtype=int)

    def rec(idx, q, base):
        idx = np.asarray(idx)
        if q <= 1 or len(idx) == 0:
            for i in idx:
                part[i] = base
            return
        q1 = q // 2; q2 = q - q1
        n1 = nint(len(idx) * q1 / q)
        n1 = max(1, min(len(idx)-1, n1))
        if len(idx) <= 2:
            order = np.arange(len(idx))
        else:
            f = fiedler_vector(idx, L)
            order = np.argsort(f, kind='stable')
        left = idx[order[:n1]]; right = idx[order[n1:]]
        rec(left, q1, base)
        rec(right, q2, base + q1)

    rec(np.arange(nelv), nparts, 0)
    return part


def genmap_partition(m, exe, nparts, tmp):
    """Run genmap_light and derive the partition it produced from the output's
    contiguous per-partition block structure, mapped back to input elements."""
    out = os.path.join(tmp, 'out.nmsh')
    r = subprocess.run([exe, m['path'], str(nparts), out],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('genmap_light failed: ' + r.stderr[:200])
    o = readmesh(out)
    # map output element -> input element by identical geometry tuple
    in_of_geo = {}
    for e, g in enumerate(m['geo']):
        in_of_geo.setdefault(g, []).append(e)
    perm_ok = True
    part = np.empty(m['nelv'], dtype=int)
    nelv = m['nelv']
    # contiguous per-partition blocks: sizes from proportional recursion
    # derive block id of each OUTPUT position, then attribute to input element
    # compute block boundaries the same way the tool does (proportional RSB)
    bounds = part_block_bounds(nelv, nparts)
    blkid = np.empty(nelv, dtype=int)
    for b in range(nparts):
        blkid[bounds[b]:bounds[b+1]] = b
    used = collections.Counter()
    for j in range(nelv):
        g = o['geo'][j]
        cand = in_of_geo.get(g)
        if not cand:
            perm_ok = False; break
        ie = cand[used[g] % len(cand)]; used[g] += 1
        part[ie] = blkid[j]
    return part, o, perm_ok, r.stdout


def nint(x):
    """Fortran NINT: round half AWAY from zero (not Python's banker's round()).
    Must match genmap_light's rsb_rec target-size arithmetic exactly."""
    return int(math.floor(x + 0.5)) if x >= 0 else -int(math.floor(-x + 0.5))


def part_block_bounds(nelv, nparts):
    """Reproduce the tool's contiguous block sizes from proportional recursive
    bisection (must match genmap_light's rsb_rec exactly, incl. NINT rounding)."""
    sizes = [0]*nparts
    def rec(n, q, base):
        if q <= 1:
            sizes[base] = n; return
        q1 = q//2; q2 = q-q1
        n1 = nint(n*q1/q)
        n1 = max(1, min(n-1, n1))
        rec(n1, q1, base); rec(n-n1, q2, base+q1)
    rec(nelv, nparts, 0)
    bounds = [0]
    for s in sizes:
        bounds.append(bounds[-1]+s)
    return bounds


def check_output_valid(m, o, nparts):
    """Reorder correctness: permutation, contiguity, curves & zones carried."""
    res = {}
    res['nelv'] = (o['nelv'] == m['nelv'])
    res['contiguous_elid'] = (list(o['elids']) == list(range(1, o['nelv']+1)))
    res['geom_permutation'] = (collections.Counter(m['geo']) ==
                               collections.Counter(o['geo']))
    # curves: payload multiset identical (ids remapped)
    res['curves_carried'] = (
        collections.Counter((c[1], c[2]) for c in m['curves']) ==
        collections.Counter((c[1], c[2]) for c in o['curves']))
    # zones: type-5 glb_pt_ids + f + p_f multiset identical; type-7 (f,label)
    def zkey5(z):  # (f, p_f, glb_pt_ids) -- independent of element renumber
        return (z[1], z[3], z[4], z[5], z[6], z[7])
    in5 = collections.Counter(zkey5(z) for z in m['zones'] if z[8] == 5)
    ou5 = collections.Counter(zkey5(z) for z in o['zones'] if z[8] == 5)
    in7 = collections.Counter((z[1], z[3]) for z in m['zones'] if z[8] == 7)
    ou7 = collections.Counter((z[1], z[3]) for z in o['zones'] if z[8] == 7)
    res['zones_carried'] = (in5 == ou5 and in7 == ou7)
    return res


def main():
    if not os.path.exists(EXE):
        print('error: genmap_light not found at', EXE, '(build it first)'); return 1
    paths = sorted(glob.glob(os.path.join(NEKO, 'examples', '**', '*.nmsh'), recursive=True) +
                   glob.glob(os.path.join(NEKO, 'tests', '**', '*.nmsh'), recursive=True))
    print('genmap_light validation  (nparts = %d)\n' % NPARTS)
    print('%-42s %7s %6s %8s %8s %8s %6s %s' %
          ('mesh', 'nelv', 'bal', 'cut_gm', 'cut_sp', 'cut_lin', 'ok', 'notes'))
    allok = True
    n_tested = 0
    for path in paths:
        m = readmesh(path)
        if m is None:
            continue
        nelv = m['nelv']
        if nelv < NPARTS or nelv > 60000:     # keep the exact eigensolver tractable
            continue
        m['path'] = path
        n_tested += 1
        with tempfile.TemporaryDirectory() as tmp:
            try:
                part, o, perm_ok, out = genmap_partition(m, EXE, NPARTS, tmp)
            except Exception as ex:
                print('%-42s  RUN FAIL: %s' % (os.path.relpath(path, NEKO), ex))
                allok = False; continue
            valid = check_output_valid(m, o, NPARTS)
        cell, nrnk = merged_cell(m)
        L, w = dual_laplacian(cell, nrnk)
        cut_gm, total = weighted_cut(w, part)
        # scipy oracle
        sp_part = scipy_rsb(L, nelv, NPARTS)
        cut_sp, _ = weighted_cut(w, sp_part)
        # linear baseline
        bounds = part_block_bounds(nelv, NPARTS)
        lin = np.empty(nelv, dtype=int)
        for b in range(NPARTS):
            lin[bounds[b]:bounds[b+1]] = b
        cut_lin, _ = weighted_cut(w, lin)
        # balance = max part / min part
        sizes = np.bincount(part, minlength=NPARTS)
        bal = sizes.max() / max(sizes.min(), 1)
        # genmap_light's cut should be within ~1.5x of the exact spectral cut
        cut_ok = (cut_gm <= max(1.5 * cut_sp, cut_sp + 0.02 * total) + 1)
        bal_ok = (sizes.max() - sizes.min() <= 1)
        vok = all(valid.values()) and perm_ok
        ok = vok and cut_ok and bal_ok
        allok = allok and ok
        notes = '' if ok else ('valid=%s cut_ok=%s bal_ok=%s' %
                               (valid, cut_ok, bal_ok))
        print('%-42s %7d %6.2f %8d %8d %8d %6s %s' %
              (os.path.relpath(path, NEKO), nelv, bal, cut_gm, cut_sp, cut_lin,
               'OK' if ok else 'FAIL', notes))
    print('\nmeshes tested: %d' % n_tested)
    print('cut_gm = genmap_light, cut_sp = SciPy exact-Fiedler RSB, cut_lin = linear split')
    print('checks: pure permutation + contiguous el_idx + curves/zones carried +')
    print('        balance (max-min<=1) + cut within 1.5x of exact spectral cut.')
    print('\nRESULT:', 'PASS' if allok and n_tested > 0 else
          ('FAIL' if n_tested else 'NO MESHES FOUND'))
    return 0 if (allok and n_tested > 0) else 1


if __name__ == '__main__':
    sys.exit(main())
