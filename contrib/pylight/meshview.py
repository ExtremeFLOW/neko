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
"""meshview -- look at a Neko .nmsh interactively.

Extracts the external surface (skin) of the mesh in numpy -- O(n^(2/3))
quads, so even huge meshes reduce to a few million faces -- and shows it
with PyVista (VTK: ParaView-grade camera, picking and clipping widgets,
`pip install pyvista`).  Colouring:

  --color zone        labelled-zone index of each boundary face (0 = none)
  --color partition   the rank that would own each element when the file is
                      read by --nparts P ranks (Neko's linear distribution)
                      -- shows partition blocks of a prepart-ed mesh
  --color jacobian    owner element's minimum Jacobian (straight-sided)

Other outputs: --export skin.vtu writes the skin for ParaView (no pyvista
needed); --screenshot out.png renders off-screen (works headless);
--matplotlib uses matplotlib instead of pyvista (small meshes only).
"""

import argparse
import sys

import numpy as np

from nekolight import (banner, read_nmsh, validate_zones, validate_curves,
                       pos_of_elid_map, min_jacobian, neko_linear_sizes)
from nekolight.view import (extract_skin, show_pyvista, show_matplotlib,
                            export_vtu)


def log(msg):
    print(msg, flush=True)


def main():
    ap = argparse.ArgumentParser(
        prog='meshview.py',
        description='View the external surface of a Neko .nmsh.')
    ap.add_argument('mesh', help='input .nmsh')
    ap.add_argument('--color', choices=('zone', 'partition', 'jacobian'),
                    default='zone', help='per-face colouring (default zone)')
    ap.add_argument('--nparts', type=int, default=0,
                    help='number of ranks for --color partition')
    ap.add_argument('--export', metavar='SKIN.VTU', default=None,
                    help='write the skin as a .vtu for ParaView and exit '
                         '(unless a viewer option is also given)')
    ap.add_argument('--screenshot', metavar='OUT.PNG', default=None,
                    help='render off-screen to an image (headless-safe)')
    ap.add_argument('--matplotlib', action='store_true',
                    help='use matplotlib instead of pyvista (small meshes)')
    args = ap.parse_args()

    log(banner('meshview'))
    log('  reading %s ...' % args.mesh)
    mesh = read_nmsh(args.mesh)
    validate_zones(mesh.nelv, mesh.zones, args.mesh)
    validate_curves(mesh.nelv, mesh.curves, args.mesh)
    pos_of_elid = pos_of_elid_map(mesh.nelv, mesh.elems)
    log('  %d hex elements, %d zones, %d curved elements'
        % (mesh.nelv, mesh.zones.shape[0], mesh.curves.shape[0]))

    # per-element data to paint on the skin
    data = {}
    z7 = mesh.zones[mesh.zones['t'] == 7]
    zone_of_facet = np.zeros((mesh.nelv, 6), dtype=np.int32)
    if z7.size:
        zone_of_facet[pos_of_elid[z7['e'].astype(np.int64)],
                      z7['f'].astype(np.int64) - 1] = z7['p_f']
    if args.color == 'partition':
        if args.nparts < 1:
            sys.exit('Error: --color partition needs --nparts P')
        # rank owning each RECORD POSITION under Neko's linear read
        bounds = np.cumsum(neko_linear_sizes(mesh.nelv, args.nparts))
        data['partition'] = np.searchsorted(bounds, np.arange(mesh.nelv),
                                            side='right').astype(np.int32)
    if args.color == 'jacobian':
        data['jacobian'] = min_jacobian(mesh.elems['v']['xyz'])

    log('  extracting the external surface ...')
    sk = extract_skin(mesh, data)
    # zone colouring is per (element, facet), painted after extraction
    sk.celldata['zone'] = zone_of_facet[sk.elem_pos, sk.facet]
    log('  skin: %d quads, %d points' % (sk.quads.shape[0],
                                         sk.points.shape[0]))

    if args.export:
        export_vtu(args.export, sk)
        log('  wrote %s (open with ParaView)' % args.export)
        if not (args.screenshot or args.matplotlib):
            return
    if args.matplotlib:
        show_matplotlib(sk, color_by=args.color, screenshot=args.screenshot)
    else:
        show_pyvista(sk, color_by=args.color, screenshot=args.screenshot,
                     title='meshview: %s' % args.mesh)


if __name__ == '__main__':
    main()
