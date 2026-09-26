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
"""Mesh partitioning: three backends over the same weighted element dual
graph, all producing block sizes that EXACTLY match how Neko distributes
elements when it reads the file.

The block model (the whole point of the reordering contract): Neko's reader
splits nelv elements over P ranks with ``linear_dist_t``
(src/common/datadist.f90): ``L = nelv // P``, ``R = nelv % P``, and the
first R ranks get L+1 elements.  The written element order must therefore be
composed of blocks of exactly those sizes, one block per rank, or a linear
read does not reproduce the partition.  The recursive backends split on
exact rank-share sums (no rounding); METIS output is repaired to the exact
sizes by moving a few boundary elements.

Backends:
  spectral   recursive spectral bisection -- Nek5000 genmap's algorithm with
             a robust eigensolver (dense / shift-invert Lanczos /
             pyamg-preconditioned LOBPCG by sub-problem size).
  metis      multilevel k-way via pymetis on the same weighted graph,
             followed by the exact-size repair pass.
  geometric  recursive coordinate bisection on element centroids: no graph,
             no eigensolve, numpy only, near-instant at any size.
"""

import sys

import numpy as np


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# Neko's linear distribution (the only block model in this package)
# ---------------------------------------------------------------------------
def neko_linear_sizes(nelv, nparts):
    """Per-rank element counts of Neko's linear_dist_t: the first
    ``nelv % nparts`` ranks get one extra element."""
    L, R = divmod(nelv, nparts)
    sizes = np.full(nparts, L, dtype=np.int64)
    sizes[:R] += 1
    return sizes


def _share(nelv, nparts, base, q1):
    """Sum of the linear shares of ranks base .. base+q1-1 (closed form)."""
    L, R = divmod(nelv, nparts)
    return q1 * L + max(0, min(base + q1, R) - base)


def weighted_cut(A, part):
    C = A.tocoo()
    mask = part[C.row] != part[C.col]
    return float(C.data[mask].sum()) / 2.0, float(C.data.sum()) / 2.0


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


def _sides_connected(subA, side):
    from scipy.sparse.csgraph import connected_components
    for s in (0, 1):
        idx = np.flatnonzero(side == s)
        if idx.size == 0:
            continue
        ncomp, _ = connected_components(subA[idx][:, idx], directed=False)
        if ncomp > 1:
            return False
    return True


def _region_grow(subA, n, n1):
    """Deterministic balanced split by BFS region growing.  Tries a BFS from
    the first node (connected side 0) and, if the complement comes out
    disconnected, the symmetric grow from the last-reached node (connected
    side 1).  On graphs where no balanced connected split exists (e.g. a
    star), returns the first split with a warning -- the reordering contract
    is unaffected, only cut quality."""
    from scipy.sparse.csgraph import breadth_first_order

    def grow(start, take_n, claimed_side):
        order = breadth_first_order(subA, start, directed=False,
                                    return_predecessors=False)
        side = np.full(n, 1 - claimed_side, dtype=np.int8)
        take = order[:min(take_n, order.size)]
        side[take] = claimed_side
        short = take_n - take.size
        if short > 0:      # disconnected input: top up deterministically
            rest = np.flatnonzero(side != claimed_side)[:short]
            side[rest] = claimed_side
        return side, order

    side, order = grow(0, n1, 0)
    if _sides_connected(subA, side):
        return side
    alt, _ = grow(int(order[-1]), n - n1, 1)
    if _sides_connected(subA, alt):
        return alt
    log('  note: no balanced connected bisection exists for a %d-element '
        'sub-graph; one side stays disconnected (cut quality only)' % n)
    return side


def spectral_partition(A, nelv, nparts, kmodes=4):
    part = np.zeros(nelv, dtype=np.int64)

    def rec(idx, q, base):
        n = idx.size
        if q <= 1 or n == 0:
            part[idx] = base
            return
        q1 = q // 2
        n1 = _share(nelv, nparts, base, q1)   # exact linear share of the left ranks
        if n <= 2 or n1 == 0 or n1 == n:
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
                conn = _sides_connected(subA, side)
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
# Backend: metis (multilevel k-way + exact-size repair)
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
        # k-way METIS can collapse parts when the mean part size gets tiny;
        # recursive bisection handles that.
        log('  note: k-way METIS returned %d non-empty parts; '
            'retrying with recursive bisection' % np.unique(part).size)
        _, membership = _part(recursive=True)
        part = np.asarray(membership, dtype=np.int64)
    if np.unique(part).size < nparts:
        sys.exit('Error: METIS produced only %d non-empty parts of %d '
                 'requested; use --backend spectral or geometric for this '
                 'nparts/nelv ratio.' % (np.unique(part).size, nparts))
    return repair_sizes(A, part, neko_linear_sizes(nelv, nparts))


def repair_sizes(A, part, target):
    """Move elements between parts until the part sizes exactly equal
    ``target`` (Neko's linear shares).  METIS is near-balanced, so this
    typically moves only a handful of elements; each move prefers a boundary
    element with the strongest ties to the destination part, keeping the
    extra cut small.  Relabelling first (largest parts onto the L+1 ranks)
    minimises the number of moves."""
    nparts = target.size
    sizes = np.bincount(part, minlength=nparts)
    # relabel: sort parts by size, ranks by target share (both descending)
    order_p = np.argsort(-sizes, kind='stable')
    order_r = np.argsort(-target, kind='stable')
    relabel = np.empty(nparts, dtype=np.int64)
    relabel[order_p] = order_r
    part = relabel[part]
    sizes = np.bincount(part, minlength=nparts)

    Ac = A.tocsr()
    moved = 0
    while True:
        over = np.flatnonzero(sizes > target)
        if over.size == 0:
            break
        src = int(over[0])
        under = np.flatnonzero(sizes < target)
        cand = np.flatnonzero(part == src)
        Asub = Ac[cand]
        # connectivity of every candidate to every under-full part: one
        # sparse matvec per destination (vectorised over the candidates)
        gains = np.stack([Asub @ (part == dst).astype(np.float64)
                          for dst in under])            # (n_under, n_cand)
        k = int(np.argmax(gains))
        dst = int(under[k // cand.size])
        e = int(cand[k % cand.size])
        part[e] = dst
        sizes[src] -= 1
        sizes[dst] += 1
        moved += 1
        if moved > len(part):
            sys.exit('Error: size-repair pass failed to converge (bug)')
    if moved:
        log('  note: moved %d element(s) to match Neko\'s exact linear '
            'block sizes' % moved)
    return part


# ---------------------------------------------------------------------------
# Backend: geometric (recursive coordinate bisection; numpy only)
# ---------------------------------------------------------------------------
def geometric_partition(cent, nelv, nparts):
    """Recursive coordinate bisection on element centroids (``cent`` is
    (nelv, 3)).  Ignores periodic wrap-around (a known quality trade-off)."""
    part = np.zeros(nelv, dtype=np.int64)

    def rec(idx, q, base):
        n = idx.size
        if q <= 1 or n == 0:
            part[idx] = base
            return
        q1 = q // 2
        n1 = _share(nelv, nparts, base, q1)
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
def reorder_and_write(path, mesh, part, pos_of_elid, inputs=()):
    """Write the reordered .nmsh: contiguous per-partition blocks (stable
    within a block), element ids renumbered 1..nelv, zone/curve element
    references remapped, all payloads verbatim."""
    from .formats import write_nmsh
    nelv = mesh.nelv
    order = np.argsort(part, kind='stable')             # new position -> old pos
    newid_of_pos = np.empty(nelv, dtype=np.int64)       # old pos -> new 1-based id
    newid_of_pos[order] = np.arange(1, nelv + 1)

    elems_out = mesh.elems[order].copy()
    elems_out['id'] = np.arange(1, nelv + 1, dtype=np.int32)

    def remap_el(ids):
        return newid_of_pos[pos_of_elid[ids.astype(np.int64)]].astype(np.int32)

    zones, curves = mesh.zones, mesh.curves
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
        # Legacy zone types (e.g. type 1/2 in older meshes).  Neko's current
        # reader ignores them, but they are part of the file -- carry them
        # through verbatim with only the element id renumbered.
        log('  note: carrying %d zone records of other types %s through '
            '(element ids renumbered, payload verbatim)'
            % (zx.shape[0], sorted(set(int(t) for t in zx['t']))))
        zx['e'] = remap_el(zx['e'])
    curves_out = curves.copy()
    if curves_out.size:
        curves_out['e'] = remap_el(curves_out['e'])

    write_nmsh(path, elems_out, (z5, z7, zx), curves_out, inputs=inputs)
