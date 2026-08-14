"""Retargeted O-grid for the SiC CVD train (from cvd_diffusor_mesh/ogrid.py).
CAD sampler removed: the superellipse fit it produced was never read, and
wall() is analytic.  Geometry retargeted to the hot-zone entry section, plus
an axisymmetric trip hump in the pipe.
"""

import gmsh, numpy as np, sys

import os
# ---- RETARGETED to the SiC CVD hot-zone train (mm) -------------------------
# Pipe radius from the exit train (rif = 15.1 mm). Channel half-dimensions are
# the hot-zone ENTRY section, 160.4 x 30.823 mm -> 80.2 x 15.4115 half.
# Flow runs -x here (ogrid convention); the mesh is rotated into the hot-zone
# frame (flow -z) on export.
R_PIPE = 15.1; Y_CH = 80.2; Z_CH = 15.4115
X_PE = -150.0; X_CS = -410.0; X_OUT = -1217.9  # pipe 150 | diffusor 260 | hot zone 607.9 | ext 200

# ---- hot-zone + outlet-extension section tables (train x, mm) ---------------
# Converted from build_hotzone_trip.py's bp table (hot-zone z, metres) via
#     x = X_CS - (539.4 - 1000*z_hz)
# so the DOWNSTREAM GEOMETRY is identical to hotzone_trip: entry z=0.5394 ->
# x=-410, heated walls z=0.285..-0.0015 -> x=-664.4..-950.9, exit z=-0.0685 ->
# x=-1017.9, plus a 200 mm outlet extension. The duct is asymmetric about the
# pipe axis (floor -14.523, ceiling +16.3/16.5), so each section is built as a
# SYMMETRIC superellipse of half-height HT and shifted vertically by CEN(x) --
# a rigid per-station translation that cannot degrade cell quality.
XT  = np.array([-1017.9, -950.9, -664.4, -607.9, -583.4, -572.4, -410.0])
YT  = np.array([  75.5,    75.5,   74.5,   75.2,   75.2,   80.2,   80.2])
FT_ = np.array([-14.523, -14.523,-14.523,  -8.4,   -8.4, -14.523,-14.523])
CT_ = np.array([  16.5,    16.5,   16.5,   11.0,   11.0,   16.3,   16.3])
HT  = 0.5*(CT_ - FT_)          # half-height
CCT = 0.5*(CT_ + FT_)          # vertical center offset
C_ENTRY = float(CCT[-1])       # 0.8885 at the hot-zone entry = diffusor exit

def _hz(x):
    return float(np.interp(x, XT, YT)), float(np.interp(x, XT, HT))

def CEN(x):
    """Vertical center of the section: 0 in the pipe, table in the hot zone,
    smootherstep blend across the diffusor."""
    if x >= X_PE: return 0.0
    if x <= X_CS: return float(np.interp(x, XT, CCT))
    return _sstep((x - X_PE)/(X_CS - X_PE))*C_ENTRY

# ---- axisymmetric trip hump in the pipe -----------------------------------
HUMP_X = -100.0          # crest, 100 mm downstream of the inlet
HUMP_L = 20.0            # half-length
HUMP_H = 3.0             # crest radius reduction: 20% of R_PIPE = 36% AREA blockage.
                         # Mesh quality is insensitive to this -- minSJ is set by
                         # the diffusor at x=-346, r=67, not by the pipe. Measured
                         # identical minSJ 0.5767 for HUMP_H = 2,3,4,5 mm, so this
                         # can go to 5 mm (55% area) without costing anything.

def _hump(x):
    t = (x - HUMP_X)/HUMP_L
    return 0.5*(1.0 + np.cos(np.pi*t)) if abs(t) < 1.0 else 0.0

def _rect_pt(th, Y, Z):
    c, s = np.cos(th), np.sin(th)
    t = min(Y/abs(c) if abs(c) > 1e-9 else 1e18, Z/abs(s) if abs(s) > 1e-9 else 1e18)
    return np.array([t*c, t*s])

def _superr(Y, Z, n, th):
    return (np.abs(np.cos(th)/Y)**n + np.abs(np.sin(th)/Z)**n)**(-1.0/n)

N_CAP = 12.0   # superellipse exponent at the channel (rounded-rect corners ~3.9mm; smooth morph)
N_INSET = 12.0 # corner exponent of the INTERFACE (inset) curve ONLY. The wall keeps
               # N_CAP, so lowering this changes NO geometry -- it just rounds the
               # butterfly interface corner, spreading its turn over more ring cells.
               # (The inset built from (Y-d, H-d) at the wall exponent is SHARPER
               # than the wall corner, which is where the worst cells came from.)

def _sstep(t): return t*t*t*(t*(t*6 - 15) + 10)   # smootherstep 6t^5-15t^4+10t^3: C2 (zero 1st&2nd deriv at ends)

def wall(x, th):
    """(y,z) on the analytic wall: circle (pipe), smoothstep superellipse morph (diffusor), rounded rect (channel)."""
    if x >= X_PE:
        rp = R_PIPE - HUMP_H*_hump(x)
        return np.array([rp*np.cos(th), rp*np.sin(th)])
    if x <= X_CS:
        Y, H = _hz(x)
        r = _superr(Y, H, N_CAP, th); return np.array([r*np.cos(th), r*np.sin(th)])
    s = _sstep((x - X_PE)/(X_CS - X_PE))
    Y = R_PIPE + s*(Y_CH - R_PIPE); Z = R_PIPE + s*(Z_CH - R_PIPE); n = 2.0 + s*(N_CAP - 2.0)
    r = _superr(Y, Z, n, th)
    return np.array([r*np.cos(th), r*np.sin(th)])

def wall_inset(x, th, d):
    """wall shape shrunk inward by BL thickness d (radially aligned with the wall -> no corner fold)."""
    if x >= X_PE:
        r = R_PIPE - HUMP_H*_hump(x) - d
        return np.array([r*np.cos(th), r*np.sin(th)])
    if x <= X_CS:
        Y, H = _hz(x)
        r = _superr(Y - d, H - d, N_INSET, th); return np.array([r*np.cos(th), r*np.sin(th)])
    s = _sstep((x - X_PE)/(X_CS - X_PE))
    Y = R_PIPE + s*(Y_CH - R_PIPE) - d; Z = R_PIPE + s*(Z_CH - R_PIPE) - d
    n = 2.0 + s*(N_INSET - 2.0)
    r = _superr(Y, Z, n, th); return np.array([r*np.cos(th), r*np.sin(th)])

def corner_angle(x):
    tc = np.arctan2(Z_CH, Y_CH)
    if x >= X_PE: return np.pi/4
    if x <= X_CS: return tc
    f = (X_PE - x)/(X_PE - X_CS)
    return (1-f)*np.pi/4 + f*tc

# ---------- helpers ----------
def lin(a, b, n):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return a[None,:] + np.linspace(0,1,n)[:,None]*(b-a)[None,:]

def arclen(dense, n):
    """resample a dense polyline to n nodes equally spaced by ARC LENGTH (even cells)."""
    sl = np.r_[0.0, np.cumsum(np.hypot(*np.diff(dense, axis=0).T))]
    tg = np.linspace(0, sl[-1], n)
    return np.c_[np.interp(tg, sl, dense[:,0]), np.interp(tg, sl, dense[:,1])]

def offset_inward(warc, d, I0, I1):
    """true parallel offset of the wall arc inward by d (concentric, correct reduced curvature),
    with endpoints snapped to the erosion vertices I0, I1."""
    t = np.gradient(warc, axis=0); t /= np.linalg.norm(t, axis=1, keepdims=True)
    n = np.c_[t[:, 1], -t[:, 0]]                              # normal
    n[np.einsum('ij,ij->i', n, warc) < 0] *= -1              # make outward
    off = warc - d*n                                          # parallel curve
    s = np.linspace(0, 1, len(warc))[:, None]
    return off + (1-s)*(I0 - off[0]) + s*(I1 - off[-1])      # snap ends to erosion vertices

def coons(e_b, e_t, e_l, e_r):
    e_b,e_t,e_l,e_r = map(lambda e: np.asarray(e,float), (e_b,e_t,e_l,e_r))
    ns, nt = len(e_b), len(e_l)
    S = np.linspace(0,1,ns)[None,:,None]; T = np.linspace(0,1,nt)[:,None,None]
    G = (1-T)*e_b[None] + T*e_t[None] + (1-S)*e_l[:,None] + S*e_r[:,None]
    G -= (1-S)*(1-T)*e_b[0] + S*(1-T)*e_b[-1] + (1-S)*T*e_t[0] + S*T*e_t[-1]
    return G  # [t, s, 2]

def join(p, q, r, n):
    """polyline p->q->r with n+1 nodes per leg (shared midpoint q)."""
    return np.vstack([lin(p,q,n+1), lin(q,r,n+1)[1:]])

def coons_p(e_b, e_t, e_l, e_r, tv, sv):
    """Coons patch with explicit param distributions tv (radial t) and sv (arc s)."""
    e_b,e_t,e_l,e_r = map(lambda e: np.asarray(e,float), (e_b,e_t,e_l,e_r))
    S = sv[None,:,None]; T = tv[:,None,None]
    G = (1-T)*e_b[None] + T*e_t[None] + (1-S)*e_l[:,None] + S*e_r[:,None]
    G -= (1-S)*(1-T)*e_b[0] + S*(1-T)*e_b[-1] + (1-S)*T*e_t[0] + S*T*e_t[-1]
    return G

def grade_wall(nr, ratio):
    """nr+1 params in [0,1], clustered toward 1 (wall): geometric, first wall-cell smallest."""
    sizes = ratio**np.arange(nr-1, -1, -1)          # large near inner -> small near wall
    pos = np.concatenate([[0.0], np.cumsum(sizes)])
    return pos/pos[-1]

# ---------- one cross-section: 8-sector octagon O-grid (graded BL ring) ----------
# nh = arc cells on the top/bottom (y-long) sectors; nv = arc cells on the side (z) sectors.
# core = (2*nv) x (2*nh) cells -> for the 5:1 channel use nh>nv to keep cells ~square.
def section8(x, nh, nv, nr, blfrac=0.12, ratio=1.25):
    """morphing octagonal butterfly (morph_butterfly.geo structure):
    4-quad butterfly core + 8-sector graded BL ring; inner octagon = wall inset by the BL thickness."""
    tc = corner_angle(x)
    A = [0, tc, np.pi/2, np.pi-tc, np.pi, np.pi+tc, 3*np.pi/2, 2*np.pi-tc]   # E,NE,N,NW,W,SW,S,SE
    narc = [nv, nh, nh, nv, nv, nh, nh, nv]        # per-sector arc cells (k=0..7)
    d = blfrac*min(wall(x, 0.0)[0], wall(x, np.pi/2)[1])           # BL ring thickness
    tv = grade_wall(nr, ratio)
    warc = []
    for k in range(8):
        a0, a1 = A[k], (A[k+1] if k < 7 else 2*np.pi)
        warc.append(arclen(np.array([wall(x, a) for a in np.linspace(a0, a1, 120)]), narc[k]+1))
    def enorm(p0, p1, ref):                                        # outward unit normal of segment p0->p1
        t = p1 - p0; t = t/np.linalg.norm(t); n = np.array([t[1], -t[0]])
        return n if n @ ref > 0 else -n
    I = []                                                         # inner octagon vertices = wall offset inward by d
    for k in range(8):
        P = warc[k][0]
        nb = enorm(warc[(k-1) % 8][-2], warc[(k-1) % 8][-1], P)    # normal of edge before vertex
        na = enorm(warc[k][0], warc[k][1], P)                      # normal of edge after vertex
        I.append(P - d*(nb + na)/(1.0 + nb @ na))                  # Minkowski erosion vertex (concentric)
    # straight octagon sides in pipe & channel; concentric splines through the diffusor transition
    if x >= X_PE:   w = 0.0
    elif x <= X_CS: w = 1.0
    else:           w = (X_PE - x)/(X_PE - X_CS)
    ie = [(1-w)*lin(I[k], I[(k+1) % 8], narc[k]+1) + w*offset_inward(warc[k], d, I[k], I[(k+1) % 8])
          for k in range(8)]
    # 4 butterfly core quads: center -> axis points, bounded by straight inner-octagon sides
    C = np.zeros(2)
    QA = coons(lin(C,I[0],nh+1), ie[1][::-1], lin(C,I[2],nv+1), ie[0])   # E-NE-N
    QB = coons(lin(C,I[2],nv+1), ie[3][::-1], lin(C,I[4],nh+1), ie[2])   # N-NW-W
    QC = coons(lin(C,I[4],nh+1), ie[5][::-1], lin(C,I[6],nv+1), ie[4])   # W-SW-S
    QD = coons(lin(C,I[6],nv+1), ie[7][::-1], lin(C,I[0],nh+1), ie[6])   # S-SE-E
    sectors = []
    for k in range(8):
        sv = np.linspace(0, 1, narc[k]+1)
        s0 = I[k]          + tv[:,None]*(warc[k][0]  - I[k])      # radial BL (graded, ~wall-normal)
        s1 = I[(k+1) % 8]  + tv[:,None]*(warc[k][-1] - I[(k+1) % 8])
        sectors.append(coons_p(ie[k], warc[k], s0, s1, tv, sv))
    return [QA, QB, QC, QD] + sectors  # 4 core + 8 ring (ring = last 8, wall at t=-1)

# ---------- x distribution ----------
def xdist(n_pipe, n_diff, n_chan, chan_grade=1.0):
    pipe = np.linspace(0, X_PE, n_pipe+1)
    diff = np.linspace(X_PE, X_CS, n_diff+1)[1:]
    if abs(chan_grade - 1.0) < 1e-9:
        chan = np.linspace(X_CS, X_OUT, n_chan+1)[1:]
    else:                                                   # geometric: fine at the diffusor exit -> coarse downstream
        sizes = chan_grade**np.arange(n_chan)
        pos = np.concatenate([[0.0], np.cumsum(sizes)]); pos /= pos[-1]
        chan = (X_CS + pos*(X_OUT - X_CS))[1:]
    return np.concatenate([pipe, diff, chan])

# ---------- full 3D hex assembly ----------
def build(xs, nh, nv, nr, fname, blfrac=0.12, ratio=1.25):
    secs = [section8(x, nh, nv, nr, blfrac=blfrac, ratio=ratio) for x in xs]
    node_id = {}; coords = []
    def gid(x, y, z):
        k = (round(x,4), round(y,4), round(z,4)); i = node_id.get(k)
        if i is None: i = len(coords); node_id[k] = i; coords.append((x,y,z))
        return i
    IDS = []
    for ix, x in enumerate(xs):
        IDS.append([np.array([[gid(x, b[it,js,0], b[it,js,1]) for js in range(b.shape[1])]
                              for it in range(b.shape[0])]) for b in secs[ix]])
    nx = len(xs); nblk = len(secs[0]); hexes = []
    for ix in range(nx-1):
        for blk in range(nblk):
            g0, g1 = IDS[ix][blk], IDS[ix+1][blk]; nt, ns = g0.shape
            for it in range(nt-1):
                for js in range(ns-1):
                    hexes.append([g0[it,js],g0[it+1,js],g0[it+1,js+1],g0[it,js+1],
                                  g1[it,js],g1[it+1,js],g1[it+1,js+1],g1[it,js+1]])
    CA = np.array(coords)
    for k, h in enumerate(hexes):
        P = CA[h]
        if np.dot(np.cross(P[1]-P[0], P[3]-P[0]), P[4]-P[0]) < 0:
            hexes[k] = [h[0],h[3],h[2],h[1],h[4],h[7],h[6],h[5]]
    def caps(ix):
        q = []
        for blk in range(nblk):
            g = IDS[ix][blk]; nt, ns = g.shape
            for it in range(nt-1):
                for js in range(ns-1):
                    q.append([g[it,js],g[it+1,js],g[it+1,js+1],g[it,js+1]])
        return q
    inlet, outlet = caps(0), caps(nx-1)
    pipe_w, diff_w, rect_w = [], [], []           # wall split by streamwise zone
    for blk in range(nblk-8, nblk):               # last 8 = ring sectors, outer edge it=nt-1
        for ix in range(nx-1):
            g0, g1 = IDS[ix][blk], IDS[ix+1][blk]; it = g0.shape[0]-1
            xmid = 0.5*(xs[ix] + xs[ix+1])
            tgt = pipe_w if xmid >= X_PE else (diff_w if xmid >= X_CS else rect_w)
            for js in range(g0.shape[1]-1):
                tgt.append([g0[it,js], g0[it,js+1], g1[it,js+1], g1[it,js]])
    with open(fname, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write('$PhysicalNames\n6\n'
                '2 1 "inlet"\n2 2 "outlet"\n2 3 "pipe_wall"\n'
                '2 4 "diffusor_wall"\n2 5 "rectangular_wall"\n3 6 "fluid"\n$EndPhysicalNames\n')
        f.write(f"$Nodes\n{len(coords)}\n")
        for i,(x,y,z) in enumerate(coords): f.write(f"{i+1} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("$EndNodes\n")
        E = ([(3,1,q) for q in inlet]+[(3,2,q) for q in outlet]
             +[(3,3,q) for q in pipe_w]+[(3,4,q) for q in diff_w]+[(3,5,q) for q in rect_w]
             +[(5,6,h) for h in hexes])
        f.write(f"$Elements\n{len(E)}\n")
        for i,(t,p,nd) in enumerate(E):
            f.write(f"{i+1} {t} 2 {p} {p} " + " ".join(str(n+1) for n in nd) + "\n")
        f.write("$EndElements\n")
    return len(coords), len(hexes), len(inlet), len(outlet), len(pipe_w), len(diff_w), len(rect_w)

# ---------- curve to true CAD profile (2nd-order, wall nodes snapped to CAD) ----------
def curve_to_cad(infile, outfile):
    """promote to 2nd order (hex27) and snap all wall nodes (corners + mid-edge/face) onto the
    analytic wall -> curved boundary elements. Wall = physical groups 3,4,5 (pipe/diffusor/rect)."""
    gmsh.initialize(); gmsh.open(infile)
    gmsh.model.mesh.setOrder(2)
    for tag in (3, 4, 5):                                          # pipe_wall, diffusor_wall, rectangular_wall
        tags, coords = gmsh.model.mesh.getNodesForPhysicalGroup(2, tag)
        coords = coords.reshape(-1, 3)
        for t, c in zip(tags, coords):
            x, y, z = c
            cz = CEN(x)                       # sections are shifted by CEN(x)
            yz = wall(x, np.arctan2(z - cz, y))
            gmsh.model.mesh.setNode(int(t), [x, yz[0], yz[1] + cz], [])
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.SaveAll", 0)                      # only physical-group elements (gmsh2nek)
    gmsh.write(outfile); gmsh.finalize()

NH, NV, NR = 30, 10, 8         # arc cells per sector (NH>NV: even cell SIZE on the 5:1 channel); radial BL layers
if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    mode = sys.argv[1] if len(sys.argv) > 1 else "build"
    if mode == "sections":
        fig, axs = plt.subplots(1, 4, figsize=(20, 5))
        for ax, x in zip(axs, [0, -150, -280, -700]):
            for B in section8(x, NH, NV, NR):
                for j in range(B.shape[0]): ax.plot(B[j,:,0], B[j,:,1], 'b-', lw=0.4)
                for i in range(B.shape[1]): ax.plot(B[:,i,0], B[:,i,1], 'b-', lw=0.4)
            ax.set_aspect('equal'); ax.set_title(f"x={x}"); ax.grid(alpha=0.3)
        fig.tight_layout(); fig.savefig("ogrid_sections.png", dpi=110); print("wrote ogrid_sections.png")
    elif mode == "curve":
        curve_to_cad("cvd.msh", "cvd_curved.msh"); print("wrote cvd_curved.msh (2nd-order, wall on CAD)")
    else:
        xs = xdist(7, 40, 95, chan_grade=1.01)   # pipe, diffusor, graded long channel (fine->coarse)
        print("nodes=%d hexes=%d inlet=%d outlet=%d pipe_wall=%d diffusor_wall=%d rect_wall=%d stations=%d"
              % (build(xs, NH, NV, NR, "cvd.msh") + (len(xs),)))
        gmsh.initialize(); gmsh.open("cvd.msh")          # heal coincident nodes -> watertight
        gmsh.option.setNumber("Geometry.Tolerance", 1e-5)
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)   # msh2 for gmsh2nek
        gmsh.write("cvd.msh"); gmsh.finalize()
        print("healed -> watertight cvd.msh (v2.2)")
