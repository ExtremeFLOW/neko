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
"""prepart -- partition a Neko .nmsh and write a reordered .nmsh whose linear
read reproduces the partition exactly (the contrib/prepart contract).

The written element order is grouped into contiguous blocks whose sizes
exactly match Neko's linear_dist_t shares (L = nelv // P, first nelv % P
ranks get L+1), so rank i of a P-rank run receives exactly partition i.
Curves and zones are carried through with only their element references
renumbered; point ids and coordinates are never changed.

Backends (--backend):
  spectral   (default) recursive spectral bisection of the periodic-merged,
             shared-vertex-weighted element dual graph -- Nek5000 genmap's
             algorithm with a robust scipy eigensolver.  Deterministic.
  metis      multilevel k-way via pymetis (pip install pymetis) on the same
             graph, repaired to the exact linear block sizes; usually the
             best cut.
  geometric  recursive coordinate bisection on element centroids: numpy
             only, no graph, near-instant on very large meshes.  Ignores
             periodic wrap-around.

Usage: prepart.py mesh.nmsh nparts [out.nmsh] [--backend ...] [--no-stats]
Default output is <base>_<nparts>.nmsh.
"""

import argparse
import sys
import time

import numpy as np

from nekolight import (banner, read_nmsh, validate_zones, validate_curves,
                       pos_of_elid_map, periodic_min_merge, compress_ids,
                       dual_graph, neko_linear_sizes, spectral_partition,
                       metis_partition, geometric_partition, weighted_cut,
                       reorder_and_write)


def log(msg):
    print(msg, flush=True)


def main():
    ap = argparse.ArgumentParser(
        prog='prepart.py',
        description='Partition a Neko .nmsh and write a reordered .nmsh '
                    'whose linear read reproduces the partition.')
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
                    help='skip the edge-cut report (with --backend geometric '
                         'this also skips the vertex merge and dual graph '
                         'entirely -- the numpy-only fast path)')
    args = ap.parse_args()

    if args.nparts < 1:
        sys.exit('Error: nparts must be a positive integer')
    if args.out is None:
        base = args.mesh
        i, j = base.rfind('.'), base.rfind('/')
        base = base[:i] if i > j else base
        args.out = '%s_%d.nmsh' % (base, args.nparts)

    log(banner('prepart  (mesh partitioner)'))
    log('  input     : %s' % args.mesh)
    log('  output    : %s' % args.out)
    log('  nparts    : %d' % args.nparts)
    log('  backend   : %s' % args.backend)

    t0 = time.time()
    log('  [1/4] reading mesh ...')
    mesh = read_nmsh(args.mesh)
    validate_zones(mesh.nelv, mesh.zones, args.mesh)
    validate_curves(mesh.nelv, mesh.curves, args.mesh)
    log('        %d hex elements, %d zones, %d curved elements'
        % (mesh.nelv, mesh.zones.shape[0], mesh.curves.shape[0]))
    if args.nparts > mesh.nelv:
        sys.exit('Error: nparts (%d) > number of elements (%d)'
                 % (args.nparts, mesh.nelv))
    pos_of_elid = pos_of_elid_map(mesh.nelv, mesh.elems)

    need_graph = (args.backend in ('spectral', 'metis')) or not args.no_stats
    A = None
    if need_graph:
        log('  [2/4] building element dual graph (periodic-merged) ...')
        vidx = mesh.elems['v']['idx'].astype(np.int64)
        cell = compress_ids(periodic_min_merge(mesh.nelv, vidx, mesh.zones,
                                               pos_of_elid))
        A = dual_graph(cell)
    else:
        log('  [2/4] skipping dual graph (--no-stats, geometric backend)')

    log('  [3/4] partitioning (%s) ...' % args.backend)
    if args.backend == 'spectral':
        part = spectral_partition(A, mesh.nelv, args.nparts)
    elif args.backend == 'metis':
        part = metis_partition(A, mesh.nelv, args.nparts)
    else:
        cent = mesh.elems['v']['xyz'].mean(axis=1)
        part = geometric_partition(cent, mesh.nelv, args.nparts)

    sizes = np.bincount(part, minlength=args.nparts)
    want = neko_linear_sizes(mesh.nelv, args.nparts)
    if not np.array_equal(np.sort(sizes), np.sort(want)) or \
       not np.array_equal(sizes, want):
        sys.exit('Error: internal error -- block sizes %s do not match '
                 'Neko\'s linear distribution %s' % (sizes, want))
    log('        part sizes: min %d / max %d (exact linear_dist blocks)'
        % (sizes.min(), sizes.max()))
    if A is not None and not args.no_stats:
        cut, total = weighted_cut(A, part)
        log('        edge cut: %d / %d shared-vertex links cut (%7.3f%%)'
            % (cut, total, 100.0 * cut / max(total, 1.0)))

    log('  [4/4] renumbering and writing reordered mesh ...')
    reorder_and_write(args.out, mesh, part, pos_of_elid,
                      inputs=(args.mesh,))
    log('  done -> %s   (%.1f s)' % (args.out, time.time() - t0))


if __name__ == '__main__':
    main()
