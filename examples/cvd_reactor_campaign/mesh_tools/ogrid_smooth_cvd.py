"""Smoothed morphing-butterfly O-grid for the CVD inlet.
Same 12-block topology as ogrid.py (4 core + 8 BL ring sectors) and same build/IO,
but the section generator is built for smooth, clean cells:
  1. core/ring interface: C1 bulged sides (cubic Hermite, tangent = wall tangent at the
     E/N/W/S rays -> the octagon's 4 axis kinks are gone); a real corner only at the
     4 diagonal block corners, with controlled interior angle 180-2*PSI; blends into
     the true wall-offset curve through the diffusor (channel behaviour unchanged).
  2. azimuthal spacing continuous across sector boundaries (no 3:1 jump between nh- and
     nv-sectors): in the pipe the corner rays sit at the uniform-index angle
     2*pi*nv/(4nh+4nv) and morph to the geometric corner of the channel; within-sector
     geometric grading matches cell sizes at shared boundaries.
  3. every streamwise blend is C2 smootherstep -> no kinked longitudinal lines at
     x=-40/-300.
Usage:
  python ogrid_smooth.py sections    # 2x4 figure: current (top) vs smooth (bottom)
  python ogrid_smooth.py zoom       # close-ups of the junction region
  python ogrid_smooth.py build      # coarse 3D build of both + minSJ comparison
"""
import numpy as np, sys
import ogrid_cvd as ogrid
from ogrid_cvd import (wall, wall_inset, X_PE, X_CS, Y_CH, Z_CH, _sstep, grade_wall, lin,
                   coons, coons_p)

# Interface tangent rotation at the diagonal corner. Swept-train tuning showed
# the optimum DIFFERS by region: the pipe Hermite wants the original 15/30
# (larger values fold it toward the wall), while the flat channel corner wants
# 20/35 (smaller core corner angle -> better corner cells). PSI is therefore
# blended with the same morph parameter s as the geometry.
PSI_S_PIPE, PSI_S_CH = np.deg2rad(15.0), np.deg2rad(20.0)
PSI_L_PIPE, PSI_L_CH = np.deg2rad(30.0), np.deg2rad(35.0)
CCH = 0.7                 # channel: azimuthal clustering factor at the 4 corner boundaries
TIE_APEX = False          # if True, the interface corner ray uses the WALL apex angle
                          # (tcI = tcW): wall and interface node patterns then migrate
                          # together through the morph -- streamwise-shear test knob.


def _wall_dense(x, a0, a1, n=241):
    return np.array([wall(x, t) for t in np.linspace(a0, a1, n)])


def _rot(v, phi):
    c, s = np.cos(phi), np.sin(phi)
    return np.array([c*v[0] - s*v[1], s*v[0] + c*v[1]])


def _apex_angle(pts, th):
    """curvature-weighted mean angle: polar angle of the corner apex of a quarter curve."""
    d1 = np.gradient(pts, axis=0); d2 = np.gradient(d1, axis=0)
    k = np.abs(d1[:, 0]*d2[:, 1] - d1[:, 1]*d2[:, 0]) / (np.linalg.norm(d1, axis=1)**3 + 1e-30)
    w = k**6 * np.linalg.norm(d1, axis=1)
    return float((th*w).sum() / w.sum())


def _hermite(P0, P1, T0, T1, n=201):
    u = np.linspace(0, 1, n)[:, None]
    h00 = 2*u**3 - 3*u**2 + 1; h10 = u**3 - 2*u**2 + u
    h01 = -2*u**3 + 3*u**2;    h11 = u**3 - u**2
    return h00*P0[None] + h10*T0[None] + h01*P1[None] + h11*T1[None]


def _wall_tangent(x, th, h=1e-4):
    t = wall(x, th + h) - wall(x, th - h)
    return t / np.linalg.norm(t)


def _wall_normal_in(x, th, h=1e-4):
    t = wall(x, th + h) - wall(x, th - h); t /= np.linalg.norm(t)
    n = np.array([t[1], -t[0]])
    return n if n @ wall(x, th) < 0 else -n


def _resample_by_sizes(dense, sizes):
    """resample dense polyline to len(sizes)+1 nodes, cell lengths proportional to sizes."""
    sl = np.r_[0.0, np.cumsum(np.hypot(*np.diff(dense, axis=0).T))]
    tg = np.r_[0.0, np.cumsum(sizes)]
    tg = tg / tg[-1] * sl[-1]
    return np.c_[np.interp(tg, sl, dense[:, 0]), np.interp(tg, sl, dense[:, 1])]


def _geo_sizes(n, d0, d1):
    """n relative cell sizes grading geometrically d0 -> d1."""
    if n == 1:
        return np.array([0.5 * (d0 + d1)])
    r = (d1 / d0) ** (1.0 / (n - 1))
    return d0 * r ** np.arange(n)


def section_smooth(x, nh, nv, nr, blfrac=0.12, ratio=1.25):
    """drop-in replacement for ogrid.section8: same block list, smooth junctions."""
    P = 4 * (nh + nv)
    tc0 = 2 * np.pi * nv / P                    # pipe corner ray: uniform-index angle
    if x >= X_PE:   s = 0.0
    elif x <= X_CS: s = 1.0
    else:           s = _sstep((X_PE - x) / (X_PE - X_CS))
    d = blfrac * min(wall(x, 0.0)[0], wall(x, np.pi/2)[1])
    tv = grade_wall(nr, ratio)
    psi_s = (1 - s) * PSI_S_PIPE + s * PSI_S_CH
    psi_l = (1 - s) * PSI_L_PIPE + s * PSI_L_CH

    # --- diagonal ray angles: pipe uniform-index -> true curvature apex of the corner,
    #     computed separately for the wall and for the interface (inset) curve ---
    if s > 0.0:
        thq = np.linspace(1e-3, np.pi/2 - 1e-3, 721)
        apW = _apex_angle(np.array([wall(x, t) for t in thq]), thq)
        apI = _apex_angle(np.array([wall_inset(x, t, d) for t in thq]), thq)
    else:
        apW = apI = tc0
    tcW = (1 - s) * tc0 + s * apW               # wall sector-boundary ray
    tcI = tcW if TIE_APEX else (1 - s) * tc0 + s * apI   # interface corner-vertex ray
    A = [0, tcW, np.pi/2, np.pi - tcW, np.pi, np.pi + tcW, 3*np.pi/2, 2*np.pi - tcW, 2*np.pi]
    B = [0, tcI, np.pi/2, np.pi - tcI, np.pi, np.pi + tcI, 3*np.pi/2, 2*np.pi - tcI, 2*np.pi]
    narc = [nv, nh, nh, nv, nv, nh, nh, nv]

    # --- wall arcs: azimuthal sizes continuous across sector boundaries ---
    dense = [_wall_dense(x, A[k], A[k+1]) for k in range(8)]
    L = [float(np.hypot(*np.diff(dn, axis=0).T).sum()) for dn in dense]
    u = [L[k] / narc[k] for k in range(8)]
    cl = (1 - s) + s * CCH                               # corner clustering (channel only)
    db = [np.sqrt(u[k-1] * u[k]) * (cl if k % 2 else 1.0) for k in range(8)]
    warc, frac = [], []
    for k in range(8):
        sz = _geo_sizes(narc[k], db[k], db[(k+1) % 8])
        warc.append(_resample_by_sizes(dense[k], sz))
        frac.append(np.r_[0.0, np.cumsum(sz)] / sz.sum())  # node fractions (for ie + sv)

    # --- block-corner vertices: blend pipe offset with channel inset (at its apex) ---
    Vp = [wall(x, B[j]) + d * _wall_normal_in(x, B[j]) for j in range(8)]
    Vi = [wall_inset(x, B[j], d) for j in range(8)]
    V = [(1 - s) * Vp[j] + s * Vi[j] for j in range(8)]

    # --- interface curve ---
    # pipe: Hermite sides, C1 at the axis rays (tangent = wall tangent), controlled
    #       corner at the 4 diagonals (rotation split short/long side, no S-bend).
    # channel: wall_inset (rounded, fold-free), ends snapped to V.
    ie = []
    for k in range(8):
        j0, j1 = k, (k + 1) % 8
        t0 = _wall_tangent(x, B[k]); t1 = _wall_tangent(x, B[k+1])
        psi = psi_s if narc[k] == nv else psi_l
        if j0 % 2 == 1: t0 = _rot(t0, +psi)              # leaving a diagonal corner
        if j1 % 2 == 1: t1 = _rot(t1, -psi)              # arriving at a diagonal corner
        ch = np.linalg.norm(V[j1] - V[j0])
        herm = _hermite(V[j0], V[j1], ch*t0, ch*t1)
        ins = np.array([wall_inset(x, t, d)
                        for t in np.linspace(B[k], B[k+1], 201)])
        w01 = np.linspace(0, 1, len(ins))[:, None]
        ins = ins + (1-w01)*(V[j0]-ins[0]) + w01*(V[j1]-ins[-1])   # snap ends to V
        cur = (1 - s) * _resample_by_sizes(herm, np.diff(frac[k])) \
              + s * _resample_by_sizes(ins, np.diff(frac[k]))
        ie.append(cur)

    # --- 4 butterfly core quads (identical structure to ogrid.section8) ---
    C = np.zeros(2)
    QA = coons(lin(C, V[0], nh+1), ie[1][::-1], lin(C, V[2], nv+1), ie[0])
    QB = coons(lin(C, V[2], nv+1), ie[3][::-1], lin(C, V[4], nh+1), ie[2])
    QC = coons(lin(C, V[4], nh+1), ie[5][::-1], lin(C, V[6], nv+1), ie[4])
    QD = coons(lin(C, V[6], nv+1), ie[7][::-1], lin(C, V[0], nh+1), ie[6])

    # --- 8 BL ring sectors: radial rays, graded to the wall ---
    sectors = []
    for k in range(8):
        s0 = V[k]           + tv[:, None] * (warc[k][0]  - V[k])
        s1 = V[(k+1) % 8]   + tv[:, None] * (warc[k][-1] - V[(k+1) % 8])
        sectors.append(coons_p(ie[k], warc[k], s0, s1, tv, frac[k]))
    return [QA, QB, QC, QD] + sectors


def build(xs, nh, nv, nr, fname, blfrac=0.12, ratio=1.25):
    """ogrid.build with the smooth section generator swapped in."""
    old = ogrid.section8
    ogrid.section8 = section_smooth
    try:
        return ogrid.build(xs, nh, nv, nr, fname, blfrac=blfrac, ratio=ratio)
    finally:
        ogrid.section8 = old


# ---------- diagnostics ----------
def _plot_blocks(ax, blocks, lw=0.4):
    for B in blocks:
        for j in range(B.shape[0]): ax.plot(B[j, :, 0], B[j, :, 1], 'b-', lw=lw)
        for i in range(B.shape[1]): ax.plot(B[:, i, 0], B[:, i, 1], 'b-', lw=lw)
    ax.set_aspect('equal'); ax.grid(alpha=0.3)


def _quad_minangle(blocks):
    """min/max interior corner-angle stats over all 2D cells of a section."""
    worst = 180.0; best = 0.0
    for B in blocks:
        p00 = B[:-1, :-1]; p10 = B[1:, :-1]; p11 = B[1:, 1:]; p01 = B[:-1, 1:]
        for a, b, c in ((p00, p01, p10), (p01, p11, p00), (p11, p10, p01), (p10, p00, p11)):
            v1 = b - a; v2 = c - a
            cs = np.einsum('...i,...i', v1, v2) / (
                np.linalg.norm(v1, axis=-1) * np.linalg.norm(v2, axis=-1) + 1e-30)
            ang = np.degrees(np.arccos(np.clip(cs, -1, 1)))
            worst = min(worst, float(ang.min())); best = max(best, float(ang.max()))
    return worst, best


NH, NV, NR = ogrid.NH, ogrid.NV, ogrid.NR

if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    mode = sys.argv[1] if len(sys.argv) > 1 else "sections"
    XS = [0, -150, -280, -700]
    if mode == "sections":
        fig, axs = plt.subplots(2, 4, figsize=(20, 10))
        for col, x in enumerate(XS):
            _plot_blocks(axs[0, col], ogrid.section8(x, NH, NV, NR))
            axs[0, col].set_title(f"current  x={x}")
            _plot_blocks(axs[1, col], section_smooth(x, NH, NV, NR))
            axs[1, col].set_title(f"smooth  x={x}")
        fig.tight_layout(); fig.savefig("sections_compare.png", dpi=110)
        print("wrote sections_compare.png")
        for x in XS:
            a0 = _quad_minangle(ogrid.section8(x, NH, NV, NR))
            a1 = _quad_minangle(section_smooth(x, NH, NV, NR))
            print(f"x={x:6.0f}  min/max cell angle  current {a0[0]:5.1f}/{a0[1]:6.1f}"
                  f"   smooth {a1[0]:5.1f}/{a1[1]:6.1f}")
    elif mode == "zoom":
        views = [(0, (0, 20), (0, 20)), (0, (-20, 0), (0, 20)),
                 (-150, (30, 60), (5, 25)), (-220, (60, 105), (0, 30))]
        fig, axs = plt.subplots(2, len(views), figsize=(5*len(views), 10))
        for col, (x, xl, yl) in enumerate(views):
            for row, gen in enumerate((ogrid.section8, section_smooth)):
                _plot_blocks(axs[row, col], gen(x, NH, NV, NR), lw=0.7)
                axs[row, col].set_xlim(*xl); axs[row, col].set_ylim(*yl)
                axs[row, col].set_title(f"{'current' if row == 0 else 'smooth'} x={x}")
        fig.tight_layout(); fig.savefig("sections_zoom.png", dpi=110)
        print("wrote sections_zoom.png")
    elif mode == "build":
        import gmsh
        xs = ogrid.xdist(7, 40, 40, chan_grade=1.02)
        for name, fn in (("cur", ogrid.build), ("smooth", build)):
            r = fn(xs, NH, NV, NR, f"cmp_{name}.msh")
            gmsh.initialize(); gmsh.open(f"cmp_{name}.msh")
            q = gmsh.model.mesh.getElementQualities(
                gmsh.model.mesh.getElementsByType(5)[0], "minSJ")
            g = gmsh.model.mesh.getElementQualities(
                gmsh.model.mesh.getElementsByType(5)[0], "gamma")
            print(f"{name:7s} hexes={r[1]}  minSJ min {q.min():.3f} mean {q.mean():.3f} "
                  f"neg% {100*(q < 0).mean():.3f}   gamma min {g.min():.3f}")
            gmsh.finalize()
