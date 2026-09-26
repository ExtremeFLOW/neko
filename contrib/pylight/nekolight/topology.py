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
"""Mesh topology: point de-duplication, periodic merges, face/edge tables,
external surface (skin) extraction and the element dual graph.

Everything is a plain function over numpy arrays.  All indices along the
8*nelv corner axis are int64 (8 * 3e8 corners overflows int32); packed sort
keys are uint64 (defined wrap-around).
"""

import sys

import numpy as np

from .formats import FACE_RE2, EDGE_RE2


# ---------------------------------------------------------------------------
# Element id maps
# ---------------------------------------------------------------------------
def pos_of_elid_map(nelv, elems):
    """0-based record position of each global element id (1-based index).

    A valid .nmsh may store its element records in any order; the ids must be
    a permutation of 1..nelv.  Anything else is a hard error.
    """
    elids = elems['id'].astype(np.int64)
    if elids.size and (elids.min() < 1 or elids.max() > nelv):
        sys.exit('Error: element id out of [1,%d]' % nelv)
    pos = np.full(nelv + 1, -1, dtype=np.int64)
    pos[elids] = np.arange(nelv)
    if (pos[1:] < 0).any():
        sys.exit('Error: element ids are not a permutation of 1..nelv '
                 '(duplicate or missing id)')
    return pos


# ---------------------------------------------------------------------------
# Point de-duplication (bit-exact, first-appearance order)
# ---------------------------------------------------------------------------
def dedup_points(xyz):
    """Assign each distinct corner coordinate a 1-based id in order of first
    appearance, comparing float64 triples bit-exactly (matching Neko's own
    point table, which hashes the raw bits).

    ``xyz`` is (nelv, 8, 3) f64.  Returns (vid (nelv, 8) int32, n_unique).

    Memory note: this views the corners as 24-byte keys and sorts them, so
    the transient cost is roughly 45 bytes per corner (int64 indices
    included).  At 3e8 elements that is a fat-node job; below ~1e7 elements
    it is instant.
    """
    n8 = xyz.shape[0] * 8
    flat = np.ascontiguousarray(xyz).reshape(n8, 3)
    keys = flat.view([('', 'V24')]).ravel()
    _, first, inverse = np.unique(keys, return_index=True,
                                  return_inverse=True)
    # renumber the (byte-sorted) unique keys to first-appearance order
    rank = np.empty(first.size, dtype=np.int64)
    rank[np.argsort(first, kind='stable')] = np.arange(first.size)
    vid = (rank[inverse] + 1).reshape(-1, 8)
    if first.size > np.iinfo(np.int32).max:
        sys.exit('Error: more than 2^31 unique points -- the .nmsh format '
                 'stores 32-bit point ids')
    return vid.astype(np.int32), int(first.size)


# ---------------------------------------------------------------------------
# Periodic merges
# ---------------------------------------------------------------------------
def periodic_min_merge(nelv, vidx, zones, pos_of_elid):
    """The vectorised min-id periodic merge used by the partitioners: every
    corner named by a periodic zone record takes ``min(own id, stored
    glb_pt_id)``, applied over the whole set at once.

    ``vidx`` is (nelv, 8) int64 vertex ids.  Returns merged ids (same shape).
    Zone records must already be validated (see formats.validate_zones).
    """
    z5 = zones[zones['t'] == 5]
    if not z5.size:
        return vidx
    uniq = np.unique(vidx)
    merge = uniq.copy()
    pos = pos_of_elid[z5['e'].astype(np.int64)]
    slots = FACE_RE2[z5['f'].astype(np.int64) - 1]      # (nz5, 4)
    raw = vidx[pos[:, None], slots]
    ridx = np.searchsorted(uniq, raw.ravel())
    np.minimum.at(merge, ridx, z5['g'].astype(np.int64).ravel())
    return merge[np.searchsorted(uniq, vidx.ravel())].reshape(nelv, 8)


def periodic_replace_merge(nelv, vidx, zones, pos_of_elid):
    """The checker's merge: corner k of a periodic facet takes glb_pt_ids(k)
    directly (last record wins), exactly like Neko's apply_periodic_facet on
    read.  Boundary-sized, so a record-order loop is the faithful (and cheap)
    implementation."""
    z5 = zones[zones['t'] == 5]
    if not z5.size:
        return vidx
    remap = {}
    for rec in z5:
        pos = pos_of_elid[int(rec['e'])]
        slots = FACE_RE2[int(rec['f']) - 1]
        for k in range(4):
            remap[int(vidx[pos, slots[k]])] = int(rec['g'][k])
    keys = np.fromiter(remap.keys(), dtype=np.int64, count=len(remap))
    vals = np.fromiter(remap.values(), dtype=np.int64, count=len(remap))
    order = np.argsort(keys)
    keys, vals = keys[order], vals[order]
    idx = np.searchsorted(keys, vidx.ravel())
    idx[idx >= keys.size] = 0
    hit = keys[idx] == vidx.ravel()
    out = vidx.ravel().copy()
    out[hit] = vals[idx[hit]]
    return out.reshape(nelv, 8)


def compress_ids(vidx):
    """Compress arbitrary vertex ids to dense 0..n-1 (sorted-id order)."""
    _, cell = np.unique(vidx.ravel(), return_inverse=True)
    return cell.reshape(vidx.shape).astype(np.int64)


# ---------------------------------------------------------------------------
# Face / edge tables (packed-key sort; shared by the checker and the skin)
# ---------------------------------------------------------------------------
def _pack_faces(vidx):
    """Canonical 16-byte key per (element, facet): the sorted 4 corner ids of
    each of the 6 facets, packed into two uint64.  Returns (6*nelv, 2)."""
    f = vidx[:, FACE_RE2].reshape(-1, 4).astype(np.uint64)   # (6n, 4)
    f.sort(axis=1)
    keys = np.empty((f.shape[0], 2), dtype=np.uint64)
    keys[:, 0] = (f[:, 0] << np.uint64(32)) | f[:, 1]
    keys[:, 1] = (f[:, 2] << np.uint64(32)) | f[:, 3]
    return keys


def face_multiplicity(vidx):
    """For every (element, facet): how many times its canonical face occurs
    in the whole mesh.  Returns (counts (nelv, 6), n_unique_faces)."""
    keys = _pack_faces(vidx).view([('', 'V16')]).ravel()
    _, inverse, counts = np.unique(keys, return_inverse=True,
                                   return_counts=True)
    return counts[inverse].reshape(-1, 6), int(counts.size)


def count_edges(vidx):
    """Number of unique edges (12 per element, sorted 2-tuples)."""
    e = vidx[:, EDGE_RE2].reshape(-1, 2).astype(np.uint64)
    lo = np.minimum(e[:, 0], e[:, 1])
    hi = np.maximum(e[:, 0], e[:, 1])
    return int(np.unique((lo << np.uint64(32)) | hi).size)


def skin(vidx):
    """External surface: the (element, facet) pairs whose face occurs exactly
    once.  Returns (elem_pos (m,), facet0 (m,)) int64 arrays.

    Use RAW (unmerged) vertex ids for viewing -- periodic boundaries are then
    part of the skin, which is what you want to look at; use merged ids to
    reproduce the checker's topological externality.
    """
    mult, _ = face_multiplicity(vidx)
    epos, fct = np.nonzero(mult == 1)
    return epos.astype(np.int64), fct.astype(np.int64)


# ---------------------------------------------------------------------------
# Element dual graph
# ---------------------------------------------------------------------------
def dual_graph(cell):
    """Weighted element dual graph A = E.E^T with the diagonal removed:
    A[i,j] = number of shared (merged) vertices between elements i and j.
    ``cell`` is (nelv, 8) dense 0-based ids.  Needs scipy."""
    import scipy.sparse as sp
    nelv = cell.shape[0]
    npts = int(cell.max()) + 1
    rows = np.repeat(np.arange(nelv, dtype=np.int64), 8)
    E = sp.csr_matrix((np.ones(nelv * 8), (rows, cell.ravel())),
                      shape=(nelv, npts))
    A = (E @ E.T).tocsr()
    A.setdiag(0)
    A.eliminate_zeros()
    return A
