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
"""mesh_checker -- validate and describe a Neko .nmsh.

Reads and validates the ENTIRE file -- elements, zones, curve records and the
exact end-of-file position -- and reports sizes (with the periodic merge
applied, matching Neko's glb_mfcs/glb_meds), the bounding box, boundary
zones and unlabelled external faces.  Any malformed record (out-of-range
reference, bad label, truncated section) is a hard error: this tool never
blesses a file it could not fully parse.

Options:
  --jacobian            also check for negative/zero Jacobians (exact for
                        straight-sided elements)
  --write-zone-indices  also write <mesh>_zone_indices.fld marking labelled
                        boundary faces by zone index (implies --jacobian)

Exit status 0 means: the file parsed completely and no check failed.
"""

import argparse
import os
import sys

import numpy as np

from nekolight import (banner, read_nmsh, validate_zones, validate_curves,
                       pos_of_elid_map, periodic_replace_merge,
                       face_multiplicity, count_edges, gll_xyz, min_jacobian,
                       facet_gll_mask, write_zone_indices_fld, FACE_RE2,
                       MAX_ZLBLS)


def log(msg):
    print(msg, flush=True)


def main():
    ap = argparse.ArgumentParser(
        prog='mesh_checker.py',
        description='Validate a Neko .nmsh and report its sizes, zones and '
                    'quality.')
    ap.add_argument('mesh', help='input .nmsh')
    ap.add_argument('--jacobian', action='store_true',
                    help='also check for negative/zero Jacobians')
    ap.add_argument('--write-zone-indices', action='store_true',
                    help='also write <mesh>_zone_indices.fld (implies '
                         '--jacobian)')
    args = ap.parse_args()
    do_jac = args.jacobian or args.write_zone_indices

    log(banner('mesh_checker'))
    failed = False

    log('  [1/3] reading and validating the file ...')
    mesh = read_nmsh(args.mesh)                  # errors on any truncation
    validate_zones(mesh.nelv, mesh.zones, args.mesh)
    validate_curves(mesh.nelv, mesh.curves, args.mesh)
    pos_of_elid = pos_of_elid_map(mesh.nelv, mesh.elems)
    if mesh.trailing:
        log('        note: %d trailing bytes past the curve section '
            '(MPI-IO no-truncate artifact; Neko ignores them)'
            % mesh.trailing)

    vidx = mesh.elems['v']['idx'].astype(np.int64)
    xyz = mesh.elems['v']['xyz']
    lo, hi = xyz.reshape(-1, 3).min(axis=0), xyz.reshape(-1, 3).max(axis=0)
    max_vidx = int(vidx.max())

    # ---- zones: facet marking + the checker's replace-merge ----
    log('  [2/3] applying the periodic merge and counting faces/edges ...')
    zones = mesh.zones
    z5 = zones[zones['t'] == 5]
    z7 = zones[zones['t'] == 7]
    zleg = zones[(zones['t'] >= 1) & (zones['t'] <= 4)]
    ftype = np.zeros((mesh.nelv, 6), dtype=np.int8)      # 0 none, 1 lbl,
    flabel = np.zeros((mesh.nelv, 6), dtype=np.int8)     # 2 per, 3 legacy
    if zleg.size:
        # legacy zone types (1..4, pre-labelled-zone Neko): the current Neko
        # reader ignores these records, but they document the boundary --
        # count their facets as accounted for, not as unlabelled
        ftype[pos_of_elid[zleg['e'].astype(np.int64)],
              zleg['f'].astype(np.int64) - 1] = 3
    if z5.size:
        ftype[pos_of_elid[z5['e'].astype(np.int64)],
              z5['f'].astype(np.int64) - 1] = 2
    labeled_cnt = np.zeros(MAX_ZLBLS + 1, dtype=np.int64)
    if z7.size:
        lbl = z7['p_f'].astype(np.int64)                 # label lives in p_f
        labeled_cnt = np.bincount(lbl, minlength=MAX_ZLBLS + 1)
        p = pos_of_elid[z7['e'].astype(np.int64)]
        f0 = z7['f'].astype(np.int64) - 1
        ftype[p, f0] = 1
        flabel[p, f0] = lbl.astype(np.int8)

    merged = periodic_replace_merge(mesh.nelv, vidx, zones, pos_of_elid)
    mult, n_faces = face_multiplicity(merged)
    n_edges = count_edges(merged)
    n_unlabeled = int(((mult == 1) & (ftype == 0)).sum())

    # ---- report (mirrors Neko's mesh_checker) ----
    log('')
    log(' --------------Size-------------')
    log(' Number of elements: %d' % mesh.nelv)
    log(' Number of points:   %d' % max_vidx)
    log(' Number of faces:    %d' % n_faces)
    log(' Number of edges:    %d' % n_edges)
    log(' Bounding box:')
    log('    x %14.6g %14.6g' % (lo[0], hi[0]))
    log('    y %14.6g %14.6g' % (lo[1], hi[1]))
    log('    z %14.6g %14.6g' % (lo[2], hi[2]))
    log('')
    log(' --------------Zones------------')
    log(' Number of periodic faces: %d' % z5.shape[0])
    if zleg.size:
        log(' Legacy zone records (types 1-4): %d (ignored by Neko\'s '
            'current reader; their facets are treated as documented '
            'boundaries here)' % zleg.shape[0])
    log('')
    log(' Labelled zones:')
    for i in range(1, MAX_ZLBLS + 1):
        if labeled_cnt[i] > 0:
            log('    Zone %2d: %d faces' % (i, labeled_cnt[i]))

    jac_min = None
    if do_jac:
        log('')
        log(' ------------Jacobian----------')
        n_bad, first_bad, jac_min = 0, 0, np.inf
        chunk = 1 << 20
        for s in range(0, mesh.nelv, chunk):
            jm = min_jacobian(xyz[s:s + chunk])
            jac_min = min(jac_min, float(jm.min()))
            bad = np.flatnonzero(jm <= 0.0)
            if bad.size:
                if n_bad == 0:
                    first_bad = s + int(bad[0]) + 1
                n_bad += int(bad.size)
        log(' Min Jacobian (straight-sided): %14.6g' % jac_min)
        if n_bad > 0:
            failed = True
            log(' Error: Found %d element(s) with a negative/zero Jacobian '
                '(first at record %d).' % (n_bad, first_bad))
        else:
            log(' No negative/zero Jacobians.')

    if n_unlabeled > 0:
        failed = True
        log(' Error: Found %d unlabelled external faces.' % n_unlabeled)

    if args.write_zone_indices:
        log('')
        log('  [3/3] writing zone-index field ...')
        base = os.path.splitext(args.mesh)[0] + '_zone_indices'
        gxyz = gll_xyz(xyz)
        mask = facet_gll_mask()                          # (6, 27)
        sc = np.zeros((mesh.nelv, 27), dtype=np.float32)
        for f0 in range(6):
            v = flabel[:, f0].astype(np.float32)[:, None]  # (nelv, 1)
            on = mask[f0][None, :] & (v > 0)
            np.maximum(sc, np.where(on, v, 0.0), out=sc)   # highest label wins
        # the ACTUAL element ids in record order -- a valid .nmsh may store
        # its records shuffled, and the id list is what maps a block to an
        # element
        write_zone_indices_fld(base, mesh.elems['id'], gxyz, sc)
        log('  wrote %s.fld (+ .nek5000 companion)' % base)
    else:
        log('')
        log('  [3/3] no field output requested')

    log(' Done')
    if failed:
        print('Mesh check failed with one or several errors.',
              file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
