#!/usr/bin/env python3
"""ONE-PIECE O-grid sweep of the full CVD train:

    inlet pipe (trip hump) -> diffusor -> hot zone -> outlet extension

Why a single sweep instead of "joining" to the existing hotzone_trip mesh:
a butterfly cross-section carries 8 singular edges (3- and 5-element edges)
running streamwise; the hot-zone tensor section has none. Singular edges
cannot terminate at a conformal interface plane, so a node-for-node join of
the two topologies is impossible. The conformal answer is to continue the
butterfly through the SAME DOWNSTREAM GEOMETRY (build_hotzone_trip's bp
table, heated zones, wafer, plus a 200 mm outlet extension) -- there is then
no joint at all. Reference for the approach: cvd_sem_33k.msh sweeps this
topology through a 250x50 channel at minSJ 0.59.

Costs vs the old downstream mesh (stated, not hidden):
  * substrate: the wafer disc cannot be a union of element faces in a swept
    section (that needed hotzone_trip's plan-view butterfly). Zone 5 is the
    set of floor faces whose centroid lies in the disc -> staircase rim at
    ~element scale. Area error is reported at build time.
  * duct corners are superellipse-rounded (~2.5 mm at n=12), not sharp.

Zones (same numbering as hotzone_trip so case files port):
  1 inlet, 2 outlet, 3 wallCold, 4 wallHot, 5 substrate, 6 fluid.
"""
import numpy as np
import gmsh
import ogrid_cvd as og
import ogrid_smooth_cvd as ogs

NH, NV, NR = 12, 5, 3            # 444 elements/section (SEM "33k" class)
# blfrac 0.15 and the morph-rate station distribution below came from the
# tuning sweeps (tune_diffusor.py): base minSJ 0.552 -> 0.684. The binding
# mechanism was STREAMWISE SHEAR through the circle->rect morph, not in-plane
# corner shape -- corner softening/refinement all made it worse.
BLFRAC, RATIO = 0.15, 1.15
OUT = "train.msh"

WAFER_X, WAFER_R = -771.4, 68.5          # wafer center/radius in train x (mm)
HOT_X0, HOT_X1 = -950.9, -664.4          # heated section (z_hz 0.285..-0.0015)


def xdist():
    """Stations pinned to every breakpoint: geometry kinks (XT), heated-zone
    limits, wafer rim. Pipe clustered on the hump; diffusor geometric, fine at
    the throat; extension coarsening outward."""
    n = 24
    t = np.linspace(0, 1, n + 1)
    xp = og.X_PE * t
    w = 1.0 + 2.0*np.exp(-((xp - og.HUMP_X)/(1.6*og.HUMP_L))**2)
    s = np.cumsum(1.0/w); s = (s - s[0])/(s[-1] - s[0])
    pipe = og.X_PE * s

    # diffusor: stations spaced uniformly in (x, morph-parameter) arc length,
    # i.e. clustered where the section changes fastest (mid-morph). 56 stations:
    # the sweep showed monotone minSJ gains 28->56 (0.552->0.628 before the
    # other knobs); shear ~ morph-per-station.
    xd = np.linspace(og.X_PE, og.X_CS, 400)
    sv = np.array([og._sstep((xx - og.X_PE)/(og.X_CS - og.X_PE)) for xx in xd])
    wv = np.hypot(np.gradient(xd), 120.0*np.gradient(sv))
    cum = np.cumsum(wv); cum -= cum[0]; cum /= cum[-1]
    diff = np.interp(np.linspace(0, 1, 57), cum, xd)[1:]

    # constriction taper (-572.4..-583.4, height 30.8->19.4 over 11 mm) refined
    # 2->5: after the diffusor fix the global minSJ moved there (0.6733 at
    # x=-581); same streamwise-shear mechanism, same cure.
    seg = [(-410.0, -572.4, 16), (-572.4, -583.4, 5), (-583.4, -607.9, 4),
           (-607.9, -664.4, 6), (-664.4, -702.9, 4), (-702.9, -839.9, 14),
           (-839.9, -950.9, 11), (-950.9, -1017.9, 6)]
    hz = np.concatenate([np.linspace(a, b, m + 1)[1:] for a, b, m in seg])

    g = 1.18 ** np.arange(8)
    pos = np.concatenate([[0.0], np.cumsum(g)]); pos /= pos[-1]
    ext = -1017.9 + pos[1:]*(og.X_OUT + 1017.9)
    return np.concatenate([pipe, diff, hz, ext])


def section_full(x, nh, nv, nr, blfrac=BLFRAC, ratio=RATIO):
    """Smooth butterfly section, rigidly shifted to the local vertical center."""
    blocks = ogs.section_smooth(x, nh, nv, nr, blfrac, ratio)
    c = og.CEN(x)
    if c:
        off = np.array([0.0, c])
        blocks = [b + off for b in blocks]
    return blocks


# ---- build ------------------------------------------------------------------
og.section8 = section_full           # og.build sweeps whatever section8 is
xs = xdist()
print(f"stations: {len(xs)}  x {xs[0]:.1f} .. {xs[-1]:.1f} mm")
r = og.build(xs, NH, NV, NR, "train_raw.msh", blfrac=BLFRAC, ratio=RATIO)
print("nodes=%d hexes=%d inlet=%d outlet=%d (wall quads %d/%d/%d by old x-split)" % r)

# ---- retag zones ------------------------------------------------------------
# og.build tags walls 3/4/5 by pipe/diffusor/channel x-position. Retag by
# physics: substrate = floor faces with centroid inside the wafer disc;
# wallHot = other wall faces in the heated x-range; wallCold = the rest.
lines = open("train_raw.msh").read().splitlines()
i_n = lines.index("$Nodes"); nn = int(lines[i_n + 1])
P = np.empty((nn + 1, 3))
for ln in lines[i_n + 2:i_n + 2 + nn]:
    p = ln.split(); P[int(p[0])] = [float(p[1]), float(p[2]), float(p[3])]

i_e = lines.index("$Elements"); ne = int(lines[i_e + 1])
counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
sub_area = 0.0
out = []
for ln in lines[i_e + 2:i_e + 2 + ne]:
    p = ln.split()
    etype, tag = int(p[1]), int(p[3])
    if etype == 3 and tag in (3, 4, 5):                    # wall quad
        nd = [int(q) for q in p[5:9]]
        c = P[nd].mean(axis=0)                             # (x, y, z)
        if c[2] < 0 and (c[0] - WAFER_X)**2 + c[1]**2 <= WAFER_R**2:
            tag = 5
            a, b_, cc, dd = P[nd]
            sub_area += 0.5*np.linalg.norm(np.cross(b_ - a, cc - a)) \
                      + 0.5*np.linalg.norm(np.cross(cc - a, dd - a))
        elif HOT_X0 <= c[0] <= HOT_X1:
            tag = 4
        else:
            tag = 3
        p[3] = p[4] = str(tag)
    counts[tag] = counts.get(tag, 0) + 1 if etype == 3 else counts.get(tag, 0)
    out.append(" ".join(p))
lines[i_e + 2:i_e + 2 + ne] = out

i_p = lines.index("$PhysicalNames")
n_p = int(lines[i_p + 1])
lines[i_p:i_p + n_p + 3] = [
    "$PhysicalNames", "6",
    '2 1 "inlet"', '2 2 "outlet"', '2 3 "wallCold"',
    '2 4 "wallHot"', '2 5 "substrate"', '3 6 "fluid"',
    "$EndPhysicalNames"]
open(OUT, "w").write("\n".join(lines) + "\n")

disc = np.pi * WAFER_R**2
print(f"zones: inlet {counts[1]}  outlet {counts[2]}  wallCold {counts[3]}  "
      f"wallHot {counts[4]}  substrate {counts[5]}")
print(f"substrate staircase: {sub_area:.0f} mm^2 vs disc {disc:.0f} mm^2  "
      f"({100*(sub_area/disc - 1):+.1f}% net area error)")

# ---- heal + quality ---------------------------------------------------------
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.open(OUT)
gmsh.option.setNumber("Geometry.Tolerance", 1e-5)
gmsh.model.mesh.removeDuplicateNodes()
tags = gmsh.model.mesh.getElementsByType(5)[0]
q = np.array(gmsh.model.mesh.getElementQualities(tags, "minSJ"))
# locate per-region minima
_, nc, _ = gmsh.model.mesh.getNodes()
PP = nc.reshape(-1, 3)
_, _, en = gmsh.model.mesh.getElements(3)
conn = np.concatenate([n for n in en]).reshape(len(tags), 8)
nt_, _, _ = gmsh.model.mesh.getNodes()
nid = {int(t): i for i, t in enumerate(nt_)}
cen = PP[np.vectorize(nid.get)(conn)].mean(axis=1)
print(f"quality: minSJ {q.min():.4f}  p0.1% {np.percentile(q, 0.1):.4f}  "
      f"median {np.median(q):.4f}  negative {(q < 0).sum()}")
for lab, m in (("pipe", cen[:, 0] >= og.X_PE),
               ("diffusor", (cen[:, 0] < og.X_PE) & (cen[:, 0] >= og.X_CS)),
               ("hot zone", (cen[:, 0] < og.X_CS) & (cen[:, 0] >= -1017.9)),
               ("outlet ext", cen[:, 0] < -1017.9)):
    print(f"  {lab:10s} n={m.sum():6d}  minSJ {q[m].min():.4f}")
i = int(np.argmin(q))
print(f"  worst cell at x={cen[i,0]:.1f}, y={cen[i,1]:.1f}, z={cen[i,2]:.1f}")
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(OUT)
gmsh.finalize()
print(f"wrote {OUT} (healed, v2.2)")
