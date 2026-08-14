#!/usr/bin/env python3
"""Build the retargeted O-grid diffusor: pipe (with trip hump) -> diffusor ->
hot-zone entry section.

Uses cvd_diffusor_mesh/ogrid_smooth.py's 12-block morphing butterfly (4 core +
8 BL ring sectors), retargeted in ogrid_cvd.py. That topology is the reason
this reaches minSJ ~0.6 where the rounded-rectangle FRAME blocking inherited
from gen_exit_train.py managed only 0.148 -- a frame has corner pads that get
crushed as the corner radius approaches the full radius, i.e. in a circle.

Streamwise distribution is clustered on the hump and on the diffusor throat,
which are the two places the section changes fastest.
"""
import numpy as np
import ogrid_cvd as og
import ogrid_smooth_cvd as ogs

NH, NV, NR = 12, 5, 3            # 444 elements/section (build_sem.py "33k" settings)
BLFRAC, RATIO = 0.2, 1.15

def xdist(n_pipe=34, n_diff=34, n_duct=10):
    """Stations: fine over the hump, fine at the throat, relaxing downstream."""
    # pipe 0 -> X_PE, clustered on the hump
    t = np.linspace(0, 1, n_pipe + 1)
    xp = og.X_PE * t
    w = 1.0 + 2.0*np.exp(-((xp - og.HUMP_X)/(1.6*og.HUMP_L))**2)   # density bump
    s = np.cumsum(1.0/w); s = (s - s[0])/(s[-1] - s[0])
    pipe = og.X_PE * s
    # diffusor: geometric, fine at the throat (X_PE end)
    g = 1.02 ** np.arange(n_diff)
    pos = np.concatenate([[0.0], np.cumsum(g)]); pos /= pos[-1]
    diff = og.X_PE + pos[1:] * (og.X_CS - og.X_PE)
    # straight duct
    duct = og.X_CS + np.linspace(0, 1, n_duct + 1)[1:] * (og.X_OUT - og.X_CS)
    return np.concatenate([pipe, diff, duct])


xs = xdist()
print(f"stations {len(xs)}: x {xs[0]:.1f} .. {xs[-1]:.1f} mm")
d = np.abs(np.diff(xs))
print(f"  dx: pipe {d[:34].min():.2f}-{d[:34].max():.2f} | "
      f"diffusor {d[34:68].min():.2f}-{d[34:68].max():.2f} | duct {d[68:].max():.2f} mm")

r = ogs.build(xs, NH, NV, NR, "diffusor.msh", blfrac=BLFRAC, ratio=RATIO)
print("nodes=%d hexes=%d inlet=%d outlet=%d pipe_w=%d diff_w=%d rect_w=%d" % r)

import gmsh
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.open("diffusor.msh")
gmsh.option.setNumber("Geometry.Tolerance", 1e-5)
gmsh.model.mesh.removeDuplicateNodes()
tags = gmsh.model.mesh.getElementsByType(5)[0]
for m in ("minSJ", "minSICN"):
    q = np.array(gmsh.model.mesh.getElementQualities(tags, m))
    print(f"  {m}: min {q.min():.4f}  p0.1% {np.percentile(q,0.1):.4f}  "
          f"median {np.median(q):.4f}  negative {(q<0).sum()}")
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("diffusor.msh")
gmsh.finalize()
print("wrote diffusor.msh")
