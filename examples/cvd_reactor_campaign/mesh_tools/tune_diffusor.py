#!/usr/bin/env python3
"""Sweep the diffusor mesh knobs against minSJ on a short pipe+diffusor domain.

Knobs (all mesh-only unless stated):
  N_INSET  corner exponent of the interface curve alone (wall geometry fixed)
  NV       azimuthal cells on the short (corner-side) sectors
  CCH      azimuthal clustering toward the corners (smaller = tighter)
  SMID     diffusor stations uniform in the MORPH parameter s instead of
           geometric from the throat (clusters where the section changes fastest)

Module attributes are patched in-process: wall_inset/section_smooth read
N_INSET/CCH from module globals at CALL time, so no reload is needed (the
function-object rebinding trap only bites `from x import f` names, and NV is
a plain argument).
"""
import numpy as np
import gmsh
import ogrid_cvd as og
import ogrid_smooth_cvd as ogs

BLFRAC, RATIO = 0.2, 1.15
X_END = -520.0


def xdist_short(smid=False, ndiff=28):
    t = np.linspace(0, 1, 17)
    pipe = og.X_PE * t
    if smid:
        xd = np.linspace(og.X_PE, og.X_CS, 400)
        sv = np.array([og._sstep((x - og.X_PE)/(og.X_CS - og.X_PE)) for x in xd])
        # arc-length in (x, s)-space with s scaled to the corner migration size
        w = np.hypot(np.gradient(xd), 120.0*np.gradient(sv))
        cum = np.cumsum(w); cum -= cum[0]; cum /= cum[-1]
        diff = np.interp(np.linspace(0, 1, ndiff + 1), cum, xd)[1:]
    else:
        g = 1.02 ** np.arange(ndiff)
        pos = np.concatenate([[0.0], np.cumsum(g)]); pos /= pos[-1]
        diff = og.X_PE + pos[1:]*(og.X_CS - og.X_PE)
    duct = og.X_CS + np.linspace(0, 1, 9)[1:]*(X_END - og.X_CS)
    return np.concatenate([pipe, diff, duct])


def run(label, n_inset=12.0, nv=5, cch=0.7, smid=False, ndiff=28,
        ratio=RATIO, blfrac=BLFRAC, nr=3, tie=False, psi=None):
    og.N_INSET = float(n_inset)
    ogs.CCH = float(cch)
    ogs.TIE_APEX = tie
    if psi is not None:
        ogs.PSI_S, ogs.PSI_L = np.deg2rad(psi[0]), np.deg2rad(psi[1])
    else:
        ogs.PSI_S, ogs.PSI_L = np.deg2rad(15.0), np.deg2rad(30.0)

    def section_full(x, nh, nv_, nr_, blfrac=blfrac, ratio=ratio):
        blocks = ogs.section_smooth(x, nh, nv_, nr_, blfrac, ratio)
        c = og.CEN(x)
        return [b + np.array([0.0, c]) for b in blocks] if c else blocks

    og.section8 = section_full
    xs = xdist_short(smid, ndiff)
    r = og.build(xs, 12, nv, nr, "_tune.msh", blfrac=blfrac, ratio=ratio)
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open("_tune.msh")
    gmsh.model.mesh.removeDuplicateNodes()
    tags = gmsh.model.mesh.getElementsByType(5)[0]
    q = np.array(gmsh.model.mesh.getElementQualities(tags, "minSJ"))
    _, nc, _ = gmsh.model.mesh.getNodes()
    P = nc.reshape(-1, 3)
    nt_, _, _ = gmsh.model.mesh.getNodes()
    nid = {int(t): i for i, t in enumerate(nt_)}
    _, _, en = gmsh.model.mesh.getElements(3)
    conn = np.concatenate([n for n in en]).reshape(len(tags), 8)
    cen = P[np.vectorize(nid.get)(conn)].mean(axis=1)
    gmsh.finalize()
    i = int(np.argmin(q))
    print(f"{label:26s} hexes {r[1]:6d}  minSJ {q.min():.4f}  p0.1% "
          f"{np.percentile(q, 0.1):.4f}  worst @ x={cen[i,0]:7.1f} "
          f"y={cen[i,1]:6.1f} z={cen[i,2]:6.1f}")
    return q.min()


print("sweep 3: winner combos")
if True:
    run("SMID56+PSI25/40", smid=True, ndiff=56, psi=(25, 40))
    run("SMID56+bl0.15", smid=True, ndiff=56, blfrac=0.15)
    run("SMID56+PSI+bl0.15", smid=True, ndiff=56, psi=(25, 40), blfrac=0.15)
    run("SMID70+PSI+bl0.15", smid=True, ndiff=70, psi=(25, 40), blfrac=0.15)
    run("SMID56+PSI20/35+bl0.15", smid=True, ndiff=56, psi=(20, 35), blfrac=0.15)
    run("SMID56+PSI30/45+bl0.15", smid=True, ndiff=56, psi=(30, 45), blfrac=0.15)
    raise SystemExit
print("sweep 2: streamwise-shear hypothesis")
run("base 28 geo")
run("ST42", ndiff=42)
run("ST56", ndiff=56)
run("SMID42", smid=True, ndiff=42)
run("SMID56", smid=True, ndiff=56)
run("TIE_APEX", tie=True)
run("RATIO 1.30", ratio=1.30)
run("blfrac 0.15", blfrac=0.15)
run("blfrac 0.30", blfrac=0.30)
run("NR4", nr=4)
run("PSI 25/40", psi=(25, 40))
print("-- combos --")
run("SMID56 + TIE", smid=True, ndiff=56, tie=True)
run("SMID56 + I8", smid=True, ndiff=56, n_inset=8)
run("SMID56 + TIE + I8", smid=True, ndiff=56, tie=True, n_inset=8)
