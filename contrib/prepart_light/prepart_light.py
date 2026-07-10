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
"""prepart_light -- a lightweight, CPU-only mesh partitioner for Neko `.nmsh`
files, in pure scientific Python (numpy/scipy).

Like Neko's `contrib/prepart` it reads an `.nmsh`, partitions it into `nparts`
parts, and writes a reordered `.nmsh` whose elements are grouped into contiguous
per-partition blocks, so a later *linear* read reproduces the partition. Curves
and zones are first-class: every curve record and every periodic/labeled zone is
carried through with only its element reference renumbered; point ids and
coordinates are never changed. It does no Neko initialisation and needs no MPI,
no ParMETIS build, and no Fortran compiler -- just Python with numpy/scipy.

Three partitioning backends (all share the same weighted element dual graph:
edge weight = number of shared vertices, with periodic facets merged through the
stored glb_pt_ids so periodic neighbours are graph-adjacent):

  --backend spectral   (default) recursive spectral bisection: the few lowest
                       eigenvectors of the dual-graph Laplacian per bisection,
                       best-of-modes balanced split with a connectivity check
                       and BFS region-growing fallback. Deterministic. This is
                       Nek5000 genmap's algorithm with a robust eigensolver
                       (scipy shift-invert Lanczos; pyamg-preconditioned LOBPCG
                       on large sub-problems when pyamg is installed).
  --backend metis      multilevel k-way via pymetis (pip install pymetis).
                       Fed the same weighted dual graph, so it minimises the
                       same edge-cut objective; usually the best cut.
  --backend geometric  recursive coordinate bisection on element centroids:
                       no graph, no eigensolve, near-instant even on very large
                       meshes. Ignores periodic wrap-around (a known quality
                       trade-off on periodic meshes).

Usage:
  prepart_light.py mesh.nmsh nparts [out.nmsh] [--backend spectral|metis|geometric]

Default output is <base>_<nparts>.nmsh (the prepart convention).

Tested against the meshes shipped with Neko and cross-checked between backends,
but it remains your responsibility to confirm the result is correct for your own
mesh -- run the bundled validator on a case you trust first. See
prepart_light.md for the algorithm, validation, and scope. 3D hex meshes only.
"""

import argparse
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components, breadth_first_order

# ---------------------------------------------------------------------------
# nmsh binary layout (packed little-endian records; see nmsh format reference)
# ---------------------------------------------------------------------------
EL_DT = np.dtype([('id', '<i4'),
                  ('v', [('idx', '<i4'), ('xyz', '<f8', (3,))], (8,))])   # 228 B
ZONE_DT = np.dtype([('e', '<i4'), ('f', '<i4'), ('p_e', '<i4'), ('p_f', '<i4'),
                    ('g', '<i4', (4,)), ('t', '<i4')])                     # 36 B
CURVE_DT = np.dtype([('e', '<i4'), ('data', '<f8', (12, 5)),
                     ('type', '<i4', (12,))])                              # 532 B
assert EL_DT.itemsize == 228 and ZONE_DT.itemsize == 36 and CURVE_DT.itemsize == 532

# Neko facet (1..6) corners as nmsh vertex slots (0-based here) -- face_nodes
# composed with the nmsh->mesh read swap [1,2,4,3,5,6,8,7]; the same table
# validated byte-exact for periodic glb_pt_ids in the sibling *_light tools.
FACE_RE2 = np.array([[1, 5, 8, 4], [2, 6, 7, 3], [1, 2, 6, 5],
                     [4, 3, 7, 8], [1, 2, 3, 4], [5, 6, 7, 8]], dtype=np.int64) - 1

BANNER = r"""
     _  __  ____  __ __  ____
    / |/ / / __/ / //_/ / __ \
   /    / / _/  / ,<   / /_/ /
  /_/|_/ /___/ /_/|_|  \____/   prepart_light  (mesh partitioner)
"""


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def read_nmsh(path):
    def count(f, what):
        a = np.fromfile(f, dtype='<i4', count=1)
        if a.size != 1:
            sys.exit('Error: truncated file (missing %s count)' % what)
        n = int(a[0])
        if n < 0:
            sys.exit('Error: negative %s count (%d) -- corrupt file?' % (what, n))
        return n

    try:
        f = open(path, 'rb')
    except OSError as ex:
        sys.exit('Error: cannot open %s (%s)' % (path, ex.strerror))
    with f:
        hdr = np.fromfile(f, dtype='<i4', count=2)
        if hdr.size != 2:
            sys.exit('Error: %s is not a Neko .nmsh file (short header)' % path)
        nelv, gdim = int(hdr[0]), int(hdr[1])
        if gdim != 3:
            sys.exit('Error: prepart_light supports 3D (hex) meshes only '
                     '(gdim=%d).' % gdim)
        elems = np.fromfile(f, dtype=EL_DT, count=nelv)
        if elems.size != nelv:
            sys.exit('Error: truncated element section (%d of %d records)'
                     % (elems.size, nelv))
        nz = count(f, 'zone')
        zones = np.fromfile(f, dtype=ZONE_DT, count=nz)
        if zones.size != nz:
            sys.exit('Error: truncated zone section (%d of %d records)'
                     % (zones.size, nz))
        nc = count(f, 'curve')
        curves = np.fromfile(f, dtype=CURVE_DT, count=nc)
        if curves.size != nc:
            sys.exit('Error: truncated curve section (%d of %d records)'
                     % (curves.size, nc))
        # NOTE: bytes may remain past the curve section (an MPI-IO no-truncate
        # artifact seen in real meshes, e.g. turb_pipe); Neko's reader ignores
        # them and so do we.
    # semantic sanity -- catches mis-framed sections that happen to parse
    # (valid Neko zone types are 1..7; 5=periodic, 7=labeled, 1..4 legacy)
    if zones.size and ((zones['t'] < 1) | (zones['t'] > 7)).any():
        sys.exit('Error: zone record with implausible type (valid: 1..7) -- '
                 'corrupt or mis-framed zone section?')
    if curves.size and ((curves['e'] < 1) | (curves['e'] > nelv)).any():
        sys.exit('Error: curve element id out of [1,nelv] -- corrupt or '
                 'mis-framed curve section?')
    return nelv, elems, zones, curves


def write_nmsh(path, elems_out, zones5, zones7, zones_other, curves_out):
    nelv = elems_out.shape[0]
    nz = zones5.shape[0] + zones7.shape[0] + zones_other.shape[0]
    try:
        f = open(path, 'wb')
    except OSError as ex:
        sys.exit('Error: cannot write %s (%s)' % (path, ex.strerror))
    with f:
        np.array([nelv, 3], dtype='<i4').tofile(f)
        elems_out.tofile(f)
        np.array([nz], dtype='<i4').tofile(f)
        zones5.tofile(f)          # periodic (type 5) first, then labeled --
        zones7.tofile(f)          # the order Neko's own writer uses --
        zones_other.tofile(f)     # then any legacy zone types, carried through
        np.array([curves_out.shape[0]], dtype='<i4').tofile(f)
        curves_out.tofile(f)


# ---------------------------------------------------------------------------
# Dual graph (periodic-merged, shared-vertex weighted)
# ---------------------------------------------------------------------------
def merged_cell(nelv, elems, zones):
    """Apply the periodic glb_pt_ids merge (exactly what Neko's reader does via
    apply_periodic_facet) and compress point ids to 0..nrnk-1. Returns the
    (nelv, 8) compressed connectivity and pos_of_elid (0-based)."""
    elids = elems['id'].astype(np.int64)
    if elids.min() < 1 or elids.max() > nelv:
        sys.exit('Error: element id out of [1,nelv]')
    pos_of_elid = np.full(nelv + 1, -1, dtype=np.int64)
    pos_of_elid[elids] = np.arange(nelv)
    if (pos_of_elid[1:] < 0).any():
        sys.exit('Error: element ids are not a permutation of 1..nelv')

    vidx = elems['v']['idx'].astype(np.int64)          # (nelv, 8)
    z5 = zones[zones['t'] == 5]
    if z5.size:
        ok = (z5['e'] >= 1) & (z5['e'] <= nelv) & (z5['f'] >= 1) & (z5['f'] <= 6)
        z5 = z5[ok]
    if z5.size:
        # Merge over the set of ids actually present (memory O(#unique ids)),
        # not a dense table sized by max id -- sparse/large ids stay cheap.
        uniq = np.unique(vidx)                                # sorted
        merge = uniq.copy()                                   # identity default
        pos = pos_of_elid[z5['e'].astype(np.int64)]
        slots = FACE_RE2[z5['f'].astype(np.int64) - 1]        # (nz5, 4)
        raw = vidx[pos[:, None], slots]                       # (nz5, 4)
        ridx = np.searchsorted(uniq, raw.ravel())
        np.minimum.at(merge, ridx, z5['g'].astype(np.int64).ravel())
        vidx = merge[np.searchsorted(uniq, vidx.ravel())].reshape(nelv, 8)

    _, cell = np.unique(vidx.ravel(), return_inverse=True)
    return cell.reshape(nelv, 8), pos_of_elid


def dual_graph(cell):
    """Weighted element dual graph A = E.Ep^T with the diagonal removed:
    A[i,j] = number of shared (merged) vertices between elements i and j."""
    nelv = cell.shape[0]
    npts = int(cell.max()) + 1
    rows = np.repeat(np.arange(nelv, dtype=np.int64), 8)
    E = sp.csr_matrix((np.ones(nelv * 8), (rows, cell.ravel())),
                      shape=(nelv, npts))
    A = (E @ E.T).tocsr()
    A.setdiag(0)
    A.eliminate_zeros()
    return A


def weighted_cut(A, part):
    C = A.tocoo()
    mask = part[C.row] != part[C.col]
    return float(C.data[mask].sum()) / 2.0, float(C.data.sum()) / 2.0


def nint(x):
    """Fortran NINT (round half away from zero) -- same proportional-split
    arithmetic as the Fortran genmap_light, so block models agree."""
    return int(np.floor(x + 0.5)) if x >= 0 else -int(np.floor(-x + 0.5))


# ---------------------------------------------------------------------------
# Backend: spectral (recursive spectral bisection, best of the low modes)
# ---------------------------------------------------------------------------
def _low_modes(subA, n, kmodes, rng_seed=12345):
    """The kmodes lowest non-trivial eigenvectors of the sub-graph Laplacian.
    Dense for small n; shift-invert Lanczos for medium; pyamg-preconditioned
    LOBPCG for large (falls back to shift-invert if pyamg is missing)."""
    from scipy.sparse.csgraph import laplacian
    L = laplacian(subA).tocsr()
    k = min(kmodes + 1, n - 1)          # +1 for the trivial constant mode
    if n <= 800:
        vals, vecs = np.linalg.eigh(L.toarray())
        return vecs[:, 1:1 + kmodes]
    v0 = np.random.default_rng(rng_seed).standard_normal(n)
    if n > 50_000:
        # shift-invert LU on a large 3D dual graph suffers fill-in explosion;
        # AMG-preconditioned LOBPCG is the scalable path (needs pyamg).
        try:
            import pyamg
            from scipy.sparse.linalg import lobpcg
            ml = pyamg.smoothed_aggregation_solver(L.tocsr().astype(np.float64))
            M = ml.aspreconditioner()
            X0 = np.random.default_rng(rng_seed).standard_normal((n, k))
            X0[:, 0] = 1.0
            vals, vecs = lobpcg(L, X0, M=M, largest=False, tol=1e-6,
                                maxiter=200)
            order = np.argsort(vals)
            return vecs[:, order[1:1 + kmodes]]
        except ImportError:
            log('  note: pyamg not installed -- falling back to shift-invert '
                'Lanczos on a %d-element subset (may be slow; pip install '
                'pyamg, or use --backend metis/geometric)' % n)
        except Exception as ex:
            log('  note: pyamg/LOBPCG failed on a %d-element subset (%s) -- '
                'falling back to shift-invert Lanczos' % (n, ex))
    from scipy.sparse.linalg import eigsh
    vals, vecs = eigsh(L.tocsc(), k=k, sigma=-1e-6, which='LM', v0=v0)
    order = np.argsort(vals)
    return vecs[:, order[1:1 + kmodes]]


def _split_connected(subA, side):
    for s in (0, 1):
        idx = np.flatnonzero(side == s)
        if idx.size == 0:
            continue
        ncomp, _ = connected_components(subA[idx][:, idx], directed=False)
        if ncomp > 1:
            return False
    return True


def _region_grow(subA, n, n1):
    """BFS from node 0 until n1 nodes are claimed; guarantees a connected first
    half on a connected subgraph. Top up deterministically if disconnected."""
    order = breadth_first_order(subA, 0, directed=False,
                                return_predecessors=False)
    side = np.ones(n, dtype=np.int8)
    take = order[:min(n1, order.size)]
    side[take] = 0
    short = n1 - take.size
    if short > 0:
        rest = np.flatnonzero(side == 1)[:short]
        side[rest] = 0
    return side


def spectral_partition(A, nelv, nparts, kmodes=4):
    part = np.zeros(nelv, dtype=np.int64)

    def rec(idx, q, base):
        n = idx.size
        if q <= 1 or n == 0:
            part[idx] = base
            return
        q1 = q // 2
        n1 = max(1, min(n - 1, nint(n * q1 / q)))
        if n <= 2:
            left, right = idx[:n1], idx[n1:]
        else:
            subA = A[idx][:, idx].tocsr()
            F = _low_modes(subA, n, min(kmodes, n - 1))
            best = None
            for c in range(F.shape[1]):
                f = F[:, c]
                j = int(np.argmax(np.abs(f)))
                if f[j] < 0:
                    f = -f                              # sign convention
                order = np.argsort(f, kind='stable')
                side = np.ones(n, dtype=np.int8)
                side[order[:n1]] = 0
                C = subA.tocoo()
                cut = float(C.data[side[C.row] != side[C.col]].sum()) / 2
                conn = _split_connected(subA, side)
                score = cut if conn else cut + 1e12     # prefer connected
                if best is None or score < best[0]:
                    best = (score, side, conn)
            side = best[1]
            if not best[2]:
                side = _region_grow(subA, n, n1)
            left = idx[side == 0]
            right = idx[side == 1]
        rec(left, q1, base)
        rec(right, q - q1, base + q1)

    rec(np.arange(nelv, dtype=np.int64), nparts, 0)
    return part


# ---------------------------------------------------------------------------
# Backend: metis (multilevel k-way, same weighted graph)
# ---------------------------------------------------------------------------
def metis_partition(A, nelv, nparts):
    try:
        import pymetis
    except ImportError:
        sys.exit('Error: --backend metis needs pymetis (pip install pymetis)')
    if nparts == 1:
        return np.zeros(nelv, dtype=np.int64)
    Ai = A.tocsr()
    if Ai.nnz >= 2**31:
        sys.exit('Error: dual graph too large for METIS 32-bit indices '
                 '(%d edges); use --backend geometric' % Ai.nnz)
    xadj = Ai.indptr.astype(np.int32)
    adjncy = Ai.indices.astype(np.int32)
    eweights = np.rint(Ai.data).astype(np.int32)

    def _part(recursive=False):
        try:                                  # pymetis >= 2025.2
            adj = pymetis.CSRAdjacency(adj_starts=xadj, adjacent=adjncy)
            return pymetis.part_graph(nparts, adjacency=adj,
                                      eweights=eweights, recursive=recursive)
        except (AttributeError, TypeError):   # older pymetis
            return pymetis.part_graph(nparts, xadj=xadj, adjncy=adjncy,
                                      eweights=eweights, recursive=recursive)

    _, membership = _part()
    part = np.asarray(membership, dtype=np.int64)
    if np.unique(part).size < nparts:
        # k-way METIS can collapse parts when the mean part size gets tiny
        # (nparts a large fraction of nelv); recursive bisection handles that.
        log('  note: k-way METIS returned %d non-empty parts; '
            'retrying with recursive bisection' % np.unique(part).size)
        _, membership = _part(recursive=True)
        part = np.asarray(membership, dtype=np.int64)
    if np.unique(part).size < nparts:
        sys.exit('Error: METIS produced only %d non-empty parts of %d '
                 'requested; use --backend spectral or geometric for this '
                 'nparts/nelv ratio.' % (np.unique(part).size, nparts))
    return part


# ---------------------------------------------------------------------------
# Backend: geometric (recursive coordinate bisection on centroids)
# ---------------------------------------------------------------------------
def geometric_partition(elems, nelv, nparts):
    cent = elems['v']['xyz'].mean(axis=1)               # (nelv, 3)
    part = np.zeros(nelv, dtype=np.int64)

    def rec(idx, q, base):
        n = idx.size
        if q <= 1 or n == 0:
            part[idx] = base
            return
        q1 = q // 2
        n1 = max(1, min(n - 1, nint(n * q1 / q)))
        c = cent[idx]
        axis = int(np.argmax(c.max(axis=0) - c.min(axis=0)))
        order = np.argsort(c[:, axis], kind='stable')
        rec(idx[order[:n1]], q1, base)
        rec(idx[order[n1:]], q - q1, base + q1)

    rec(np.arange(nelv, dtype=np.int64), nparts, 0)
    return part


# ---------------------------------------------------------------------------
# Reorder + write (the prepart contract)
# ---------------------------------------------------------------------------
def reorder_and_write(path, nelv, elems, zones, curves, part, pos_of_elid):
    # stable sort by partition: contiguous blocks, original order within a part
    order = np.argsort(part, kind='stable')             # new position -> old pos
    newid_of_pos = np.empty(nelv, dtype=np.int64)       # old pos -> new 1-based id
    newid_of_pos[order] = np.arange(1, nelv + 1)

    elems_out = elems[order].copy()
    elems_out['id'] = np.arange(1, nelv + 1, dtype=np.int32)

    def remap_el(ids):
        ids = ids.astype(np.int64)
        ok = (ids >= 1) & (ids <= nelv)
        out = ids.copy()
        out[ok] = newid_of_pos[pos_of_elid[ids[ok]]]
        return out.astype(np.int32)

    z5 = zones[zones['t'] == 5].copy()
    z7 = zones[zones['t'] == 7].copy()
    zx = zones[(zones['t'] != 5) & (zones['t'] != 7)].copy()
    if z5.size:
        z5['e'] = remap_el(z5['e'])
        z5['p_e'] = remap_el(z5['p_e'])
    if z7.size:
        z7['e'] = remap_el(z7['e'])
        z7['p_e'] = 0                                   # Neko leaves these unset;
        z7['g'] = 0                                     # write deterministic zeros
    if zx.size:
        # Legacy zone types (e.g. type 1/2 in older meshes: poisson, nekbone).
        # Neko's current reader ignores them, but they are part of the file --
        # carry them through verbatim with only the element id renumbered.
        log('  note: carrying %d zone records of other types %s through '
            '(element ids renumbered, payload verbatim)'
            % (zx.shape[0], sorted(set(int(t) for t in zx['t']))))
        zx['e'] = remap_el(zx['e'])
    curves_out = curves.copy()
    if curves_out.size:
        curves_out['e'] = remap_el(curves_out['e'])

    write_nmsh(path, elems_out, z5, z7, zx, curves_out)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        prog='prepart_light.py',
        description='Partition a Neko .nmsh and write a reordered .nmsh whose '
                    'linear read reproduces the partition (like prepart).')
    ap.add_argument('mesh', help='input .nmsh')
    ap.add_argument('nparts', type=int, help='number of partitions')
    ap.add_argument('out', nargs='?', default=None,
                    help='output .nmsh (default <base>_<nparts>.nmsh)')
    ap.add_argument('--backend', choices=('spectral', 'metis', 'geometric'),
                    default='spectral')
    ap.add_argument('--metis', dest='backend', action='store_const',
                    const='metis', help='shorthand for --backend metis')
    ap.add_argument('--geometric', dest='backend', action='store_const',
                    const='geometric', help='shorthand for --backend geometric')
    ap.add_argument('--no-stats', action='store_true',
                    help='skip the edge-cut report (avoids building the dual '
                         'graph for --backend geometric on very large meshes)')
    args = ap.parse_args()

    if args.nparts < 1:
        sys.exit('Error: nparts must be a positive integer')
    if args.out is None:
        base = args.mesh
        i = base.rfind('.')
        j = base.rfind('/')
        base = base[:i] if i > j else base
        args.out = '%s_%d.nmsh' % (base, args.nparts)

    log(BANNER)
    log('  input     : %s' % args.mesh)
    log('  output    : %s' % args.out)
    log('  nparts    : %d' % args.nparts)
    log('  backend   : %s' % args.backend)

    t0 = time.time()
    log('  [1/4] reading mesh ...')
    nelv, elems, zones, curves = read_nmsh(args.mesh)
    log('        %d hex elements, %d zones, %d curved elements'
        % (nelv, zones.shape[0], curves.shape[0]))
    if args.nparts > nelv:
        sys.exit('Error: nparts (%d) > number of elements (%d)'
                 % (args.nparts, nelv))

    cell, pos_of_elid = merged_cell(nelv, elems, zones)

    need_graph = (args.backend in ('spectral', 'metis')) or not args.no_stats
    A = None
    if need_graph:
        log('  [2/4] building element dual graph (periodic-merged) ...')
        A = dual_graph(cell)
    else:
        log('  [2/4] skipping dual graph (--no-stats, geometric backend)')

    log('  [3/4] partitioning (%s) ...' % args.backend)
    if args.backend == 'spectral':
        part = spectral_partition(A, nelv, args.nparts)
    elif args.backend == 'metis':
        part = metis_partition(A, nelv, args.nparts)
    else:
        part = geometric_partition(elems, nelv, args.nparts)

    sizes = np.bincount(part, minlength=args.nparts)
    log('        part sizes: min %d / max %d' % (sizes.min(), sizes.max()))
    if A is not None and not args.no_stats:
        cut, total = weighted_cut(A, part)
        log('        edge cut: %d / %d shared-vertex links cut (%7.3f%%)'
            % (cut, total, 100.0 * cut / max(total, 1.0)))

    log('  [4/4] renumbering and writing reordered mesh ...')
    reorder_and_write(args.out, nelv, elems, zones, curves, part, pos_of_elid)
    log('  done -> %s   (%.1f s)' % (args.out, time.time() - t0))


if __name__ == '__main__':
    main()
