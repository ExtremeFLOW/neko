#!/usr/bin/env python3
"""Full CVD train as a STANDALONE mesh: inlet pipe (with trip hump) -> CAD
diffusor -> hot zone -> extended outlet.

Adapted from CVD_showerhead/gen_exit_train.py, which emitted an .inc to be
Included by sh_param.geo. Changes:
  * SELF-CONTAINED: station 0 builds its own ring instead of borrowing
    pI0..pI11 from the showerhead geo, and we emit a complete .geo with
    physical groups rather than an include fragment.
  * INLET PIPE EXTENDED 40 -> L_PIPE mm, carrying an axisymmetric trip hump
    (radius reduction) so the flow is disturbed before the diffusor.
  * OUTLET EXTENDED by L_OUT so hot gas has room to leave without the
    outflow BC sitting on the thermal plume.
  * Walls split into wallHot / wallCold by station, plus the wafer substrate,
    so the case file can drive them separately.

Original docstring follows.
---
Emit exit_train.inc: exit pipe -> CAD diffusor -> hot-zone duct as ONE
O-grid/frame hex blocking, copying the reference methodology:
  - ogrid_smooth.py (cvd_diffusor_mesh): O-grid ring in the round pipe that
    morphs into a Cartesian core with BL bands on all four walls
  - build_ogrid_v2.py (sic_cvd_hot_zone): side strips, wafer butterfly tied
    to the core corners, bowed box edges
Cross-section topology everywhere: 3x3 frame = core + 4 wall BL bands + 4
corner pads (12 outer pts + 4 core corners). The diffusor blends the CAD
rings (cad_rings.npz) into the EXACT hot-zone entry rectangle over the last
30%% of the loft. The wafer box replaces core+floor/ceiling bands with a
4-level butterfly stack; pads and side bands sweep straight through.

Expects geo vars pI0..pI11 / aI0..aI11 (12-split interface ring). Appends
to sExit[], sWafDisk[], sOut[], sIface[], vAll[]. Regenerate after changing
geometry parameters. Counts: core 12x12, bands 3 layers ratio 1.2 ->
~324 cells/section (their SEM 328/section class).
"""
import numpy as np

ax, ay = 0.0, 0.0013
z0, Lexit, Ldiff, rif = 1.0040, 0.040, 0.260, 0.0151
z2, z3 = z0 - Lexit, z0 - Lexit - Ldiff
hzS = [0.0, 0.1624, 0.1734, 0.1979, 0.2544, 0.2884, 0.4344, 0.5409, 0.6079]
hzX = [0.0802, 0.0802, 0.0752, 0.0752, 0.0745, 0.074619, 0.075128,
       0.0755, 0.0755]
hzF = [-0.014523]*2 + [-0.0084]*2 + [-0.014523]*5
hzC = [0.0163, 0.0163, 0.0110, 0.0110, 0.0165, 0.0165, 0.0165, 0.0165,
       0.0165]
BOXA, BOXB = 5, 6
rWaf, sWaf, bow = 0.0685, 0.3614, 0.010
QC, QM = 0.60, 0.48

# ---------------------------------------------------------------- NEW GEOMETRY
L_PIPE  = 0.150          # inlet pipe length (was Lexit = 0.040)
HUMP_S  = 0.100          # hump centre, measured upstream from the pipe exit z2
HUMP_L  = 0.020          # hump half-length
HUMP_H  = 0.0020         # crest radius reduction (13% of rif) -- the trip
NPIPE_ST = 15            # stations along the pipe (must resolve the hump)
L_OUT   = 0.200          # straight outlet extension past the hot zone
NOUT_ST = 5              # stations along it
NPIPE_SEG, NOUT_SEG = 3, 5   # pts per segment in pipe / outlet

def pipe_r(z):
    """Pipe radius with an axisymmetric raised-cosine hump (trip)."""
    t = (z - (z2 + HUMP_S))/HUMP_L
    b = 0.5*(1.0 + np.cos(np.pi*t)) if abs(t) < 1.0 else 0.0
    return rif - HUMP_H*b

DELTA = 0.0025           # BL band width (= hot-zone strip)
N = 7
NCW = 2*N - 1            # pts on core edges (12 cells)
NRB = 4                  # pts across BL band (3 layers)
RATIO = 1.2              # BL growth ratio, fine at the wall
NRAD, NDG, DPR = 4, 5, 1.2
NPIPE, NCAD = 6, 2
sc = (N - 1)/14
hzNp = [round(20*sc)+1, round(6*sc)+1, round(10*sc)+1, round(14*sc)+1,
        round(10*sc)+1, NCW, round(20*sc)+1, round(12*sc)+1]

L, cnt = [], {"p": 0, "c": 0, "f": 0, "v": 0, "ll": 0, "sl": 0}
def nm(k):
    cnt[k] += 1
    return f"{k}{cnt[k]}"

def P(x, y, z):
    p = nm("p")
    L.append(f"{p} = newp; Point({p}) = {{{x:.9g}, {y:.9g}, {z:.9g}, lc}};")
    return p

def LINE(a, b):
    c = nm("c")
    L.append(f"{c} = newl; Line({c}) = {{{a}, {b}}};")
    return c

def SPL(pts):
    if len(pts) == 2:
        return LINE(pts[0], pts[1])
    c = nm("c")
    L.append(f"{c} = newl; Spline({c}) = {{{', '.join(pts)}}};")
    return c

def neg(x):
    return x[1:] if x.startswith("-") else "-" + x

def TC(curves, n, grad=""):
    if curves:
        L.append(f"Transfinite Curve {{{', '.join(curves)}}} = {n}{grad};")

def FACE(loop, corners=None, plane=False):
    ll, f = nm("ll"), nm("f")
    L.append(f"{ll} = newll; Line Loop({ll}) = {{{', '.join(loop)}}};")
    L.append(f"{f} = news; {'Plane Surface' if plane else 'Surface'}({f}) = "
             f"{{{ll}}}; Recombine Surface {{{f}}};")
    if corners:
        L.append(f"Transfinite Surface {{{f}}} = {{{', '.join(corners)}}};")
    else:
        L.append(f"Transfinite Surface {{{f}}};")
    return f

def VOL(faces, corners):
    sl, v = nm("sl"), nm("v")
    L.append(f"{sl} = newsl; Surface Loop({sl}) = {{{', '.join(faces)}}};")
    L.append(f"{v} = newv; Volume({v}) = {{{sl}}};")
    L.append(f"Transfinite Volume {{{v}}} = {{{', '.join(corners)}}};")
    L.append(f"vAll[] += {{{v}}};")

def quarter_walk(a, ct, fb, rc, m, s):
    """point at arclength s from corner-mid m along quarter m (CCW) of the
    rounded rect: half-width a, top ct, bottom fb, corner radius rc.
    Quarter m runs corner-arc-mid m -> corner-arc-mid m+1."""
    cc = [np.array([a - rc, ct - rc]), np.array([-(a - rc), ct - rc]),
          np.array([-(a - rc), -(fb - rc)]), np.array([a - rc, -(fb - rc)])]
    amid = np.deg2rad(45 + 90*m)
    Larc = rc*np.pi/4                       # half-corner arc
    if m in (0, 2):                         # top / bottom straight
        Lst = 2*(a - rc)
    else:
        Lst = (ct - rc) + (fb - rc)
    mp = (m + 1) % 4
    if s <= Larc:                           # on corner-arc m, past its mid
        th = amid + (s/rc if rc > 0 else 0)
        return cc[m] + rc*np.array([np.cos(th), np.sin(th)])
    if s <= Larc + Lst:                     # on the straight edge
        u = s - Larc
        p0 = cc[m] + rc*np.array([np.cos(amid + np.pi/4),
                                  np.sin(amid + np.pi/4)])
        dv = np.array([[-1, 0], [0, -1], [1, 0], [0, 1]][m])
        return p0 + u*dv
    u = s - Larc - Lst                      # on corner-arc m+1, before mid
    th0 = np.deg2rad(90*mp)
    th = th0 + (u/rc if rc > 0 else 0)
    return cc[mp] + rc*np.array([np.cos(th), np.sin(th)])

def quarter_len(a, ct, fb, rc, m):
    return rc*np.pi/2 + (2*(a - rc) if m in (0, 2) else (ct - rc) + (fb - rc))

# ------------------------------------------------------------- stations
# each station: analytic rounded-rect (a, ct, fb, rc) at height z (+bow)
stations = []
def add_station(z, a, ct, fb, rc, bw=0, rot=0.0):
    stations.append(dict(z=z, a=a, ct=ct, fb=fb, rc=rc, bow=bw, rot=rot))

# ROT was -pi/8 in the showerhead build, to line the pipe frame up with the 8
# bores. Standalone there are no bores, and a non-zero ROT twists the diffusor
# (the intermediate rounded-rectangles rotate, pushing the y-extent from
# -0.013..0.018 out to -0.021..0.025). Straight duct -> no rotation.
ROT = 0.0
# --- inlet pipe, z2+L_PIPE -> z2, circular, radius humped by pipe_r ---
for i in range(NPIPE_ST):
    zp = z2 + L_PIPE*(1.0 - i/(NPIPE_ST - 1.0))
    r = pipe_r(zp)
    add_station(zp, r, r, r, r, rot=ROT)
XHe, Fe, Ce = hzX[0], hzF[0], hzC[0]
NDIF = 12
for i in range(NDIF):
    s = (i + 1.0)/(NDIF + 1)
    h = s*s*s*(10 + s*(-15 + 6*s))           # smootherstep: clean junctions
    add_station(z2 - s*Ldiff,
                rif + h*(XHe - rif), rif + h*(Ce - rif),
                rif + h*(-Fe - rif), rif*(1 - h), rot=ROT*(1 - h))
for i in range(9):
    bw = bow if i == BOXA else (-bow if i == BOXB else 0)
    add_station(z3 - hzS[i], hzX[i], hzC[i], -hzF[i], 0.0, bw)
IHZ_END = len(stations) - 1              # last hot-zone station
# --- straight outlet extension at the hot-zone exit section ---
for i in range(1, NOUT_ST + 1):
    add_station(z3 - hzS[-1] - L_OUT*i/NOUT_ST,
                hzX[-1], hzC[-1], -hzF[-1], 0.0)
NST = len(stations)
IHZ0 = IHZ_END - 8                       # first hot-zone station
IBA = IHZ0 + BOXA
IPIPE_END = NPIPE_ST - 1                 # station at z2
seg_n = ([NPIPE_SEG]*IPIPE_END + [NCAD]*(IHZ0 - IPIPE_END)
         + hzNp + [NOUT_SEG]*NOUT_ST)
assert len(seg_n) == NST - 1, (len(seg_n), NST - 1)

# -------------------------------------------------------- emit sections
S = []
for k, st in enumerate(stations):
    z, a, ct, fb, rc, bw = (st["z"], st["a"], st["ct"], st["fb"], st["rc"],
                            st["bow"])
    rot = st["rot"]
    cR, sR = np.cos(rot), np.sin(rot)
    def RT(xy):
        return (xy[0]*cR - xy[1]*sR, xy[0]*sR + xy[1]*cR)
    ai, cti, fbi, rci = a - DELTA, ct - DELTA, fb - DELTA, max(rc - DELTA, 0)
    o, e = [], []
    if True:                     # every station builds its own ring now
        for m in range(4):
            Lq = quarter_len(a, ct, fb, rc, m)
            for s_ in [0.0, DELTA, Lq - DELTA]:
                x, y = RT(quarter_walk(a, ct, fb, rc, m, s_))
                o.append(P(ax + x, ay + y, z))
        for j in range(12):
            m, piece = j//3, j % 3
            Lq = quarter_len(a, ct, fb, rc, m)
            segs = [(0.0, DELTA, 3), (DELTA, Lq - DELTA, 9),
                    (Lq - DELTA, Lq, 3)][piece]
            s0_, s1_, npts = segs
            ss = np.linspace(s0_, s1_, npts + 2)[1:-1]
            main = piece == 1
            pts = [o[j]]
            for si_, u in zip(ss, np.linspace(0, 1, npts + 2)[1:-1]):
                x, y = RT(quarter_walk(a, ct, fb, rc, m, si_))
                w = bw*np.sin(np.pi*u) if (main and bw and m in (0, 2)) else 0
                pts.append(P(ax + x, ay + y, z + w))
            pts.append(o[(j + 1) % 12])
            e.append(SPL(pts))
    q = []
    for m in range(4):
        x, y = RT(quarter_walk(ai, cti, fbi, rci, m, 0.0))
        q.append(P(ax + x, ay + y, z))
    iq = []
    for m in range(4):
        Lq = quarter_len(ai, cti, fbi, rci, m)
        pts = [q[m]]
        for si_, u in zip(np.linspace(0, Lq, 11)[1:-1],
                          np.linspace(0, 1, 11)[1:-1]):
            x, y = RT(quarter_walk(ai, cti, fbi, rci, m, si_))
            w = bw*np.sin(np.pi*u) if (bw and m in (0, 2)) else 0
            pts.append(P(ax + x, ay + y, z + w))
        pts.append(q[(m + 1) % 4])
        iq.append(SPL(pts))
    ka = [LINE(o[3*m + 1], q[m]) for m in range(4)]
    kb = [LINE(o[(3*m + 11) % 12], q[m]) for m in range(4)]
    TC([e[3*m + 1] for m in range(4)] + iq, NCW)
    TC([e[3*m] for m in range(4)], NRB, f" Using Progression {RATIO:.6g}")
    TC([e[3*m + 2] for m in range(4)], NRB,
       f" Using Progression {1/RATIO:.6g}")
    TC(ka + kb, NRB, f" Using Progression {RATIO:.6g}")
    pads = [FACE([e[3*m], ka[m], f"-{kb[m]}", e[(3*m + 11) % 12]])
            for m in range(4)]
    bands = [FACE([e[3*m + 1], kb[(m+1) % 4], f"-{iq[m]}", f"-{ka[m]}"])
             for m in range(4)]
    core = FACE(iq, corners=q)
    S.append(dict(o=o, e=e, q=q, iq=iq, ka=ka, kb=kb, pads=pads,
                  bands=bands, core=core))
L.append("sIface[] = {" + ", ".join(S[0]["pads"] + S[0]["bands"]
                                    + [S[0]["core"]]) + "};")
# export interface-section entities for the chamber generator
L.append("ifQ[] = {" + ", ".join(S[0]["q"]) + "};")
L.append("ifIQ[] = {" + ", ".join(S[0]["iq"]) + "};")
L.append("ifKA[] = {" + ", ".join(S[0]["ka"]) + "};")
L.append("ifKB[] = {" + ", ".join(S[0]["kb"]) + "};")
L.append("ifPads[] = {" + ", ".join(S[0]["pads"]) + "};")
L.append("ifBands[] = {" + ", ".join(S[0]["bands"]) + "};")
L.append("ifCore = " + S[0]["core"] + ";")

# ------------------------------------------------------------- segments
walls, walls_hot = [], []
# Heated section: hot-zone stations 4..7 (z_train 0.4496..0.1631), which is
# z_hz = 0.285..-0.0015 in the hot-zone frame -- i.e. REACTOR_K = {1,2,3} of
# build_ogrid_v2.py. Segments IHZ0+4 .. IHZ0+6 lie inside it.
HOT_SEG = set(range(IHZ0 + 4, IHZ0 + 7))
for k in range(NST - 1):
    A, B = S[k], S[k+1]
    box = (k == IBA)
    ro = [LINE(A["o"][j], B["o"][j]) for j in range(12)]
    rq = [LINE(A["q"][m], B["q"][m]) for m in range(4)]
    gr = ""
    if k == IBA - 1:
        gr = f" Using Progression {1/1.12:.6g}"
    if k == IBA + 1:
        gr = " Using Progression 1.12"
    TC(ro + rq, seg_n[k], gr)
    fw = {}
    for j in range(12):
        if box and j in (1, 7):
            continue                       # main floor/ceiling -> butterfly
        fw[j] = FACE([A["e"][j], ro[(j+1) % 12], f"-{B['e'][j]}", f"-{ro[j]}"])
        (walls_hot if k in HOT_SEG else walls).append(fw[j])
    fiq, fka, fkb = {}, {}, {}
    for m in range(4):
        if not (box and m in (0, 2)):
            fiq[m] = FACE([A["iq"][m], rq[(m+1) % 4], f"-{B['iq'][m]}",
                           f"-{rq[m]}"])
        fka[m] = FACE([A["ka"][m], rq[m], f"-{B['ka'][m]}", f"-{ro[3*m + 1]}"])
        fkb[m] = FACE([A["kb"][m], rq[m], f"-{B['kb'][m]}",
                       f"-{ro[(3*m + 11) % 12]}"])
    for m in range(4):
        VOL([A["pads"][m], B["pads"][m], fw[3*m], fw[(3*m + 11) % 12],
             fka[m], fkb[m]],
            [A["o"][3*m], A["o"][3*m + 1], A["q"][m], A["o"][(3*m + 11) % 12],
             B["o"][3*m], B["o"][3*m + 1], B["q"][m], B["o"][(3*m + 11) % 12]])
    for m in range(4):
        if box and m in (0, 2):
            continue
        mp = (m+1) % 4
        VOL([A["bands"][m], B["bands"][m], fw[3*m + 1], fiq[m], fka[m],
             fkb[mp]],
            [A["o"][3*m + 1], A["o"][3*m + 2], A["q"][mp], A["q"][m],
             B["o"][3*m + 1], B["o"][3*m + 2], B["q"][mp], B["q"][m]])
    if not box:
        VOL([A["core"], B["core"]] + [fiq[m] for m in range(4)],
            A["q"] + B["q"])
    else:
        BX = dict(A=A, B=B, fiq=fiq, fka=fka, fkb=fkb, ro=ro, rq=rq)

# ------------------------------------------------------ wafer butterfly
A, B = BX["A"], BX["B"]
F5, C5 = hzF[BOXA], hzC[BOXA]
zwc = z3 - sWaf
qw, qs, qm = rWaf/np.sqrt(2), QC*rWaf/np.sqrt(2), QM*rWaf
qco = (2*qs*qs - qm*qm)/(2*(qs - qm))
levels = [F5, F5 + DELTA, C5 - DELTA, C5]
lvc = [[A["o"][8], A["o"][7], B["o"][7], B["o"][8]],      # UR UL DL DR
       [A["q"][3], A["q"][2], B["q"][2], B["q"][3]],
       [A["q"][0], A["q"][1], B["q"][1], B["q"][0]],
       [A["o"][1], A["o"][2], B["o"][2], B["o"][1]]]
lvU = [f"-{A['e'][7]}", f"-{A['iq'][2]}", A["iq"][0], A["e"][1]]  # UR->UL
lvD = [B["e"][7], B["iq"][2], f"-{B['iq'][0]}", f"-{B['e'][1]}"]  # DL->DR
lvR = [BX["ro"][8], BX["rq"][3], BX["rq"][0], BX["ro"][1]]        # UR->DR
lvL = [BX["ro"][7], BX["rq"][2], BX["rq"][1], BX["ro"][2]]        # UL->DL
lvV = [[A["kb"][3], A["ka"][2], B["ka"][2], B["kb"][3]],          # lower->up
       [A["iq"][3], f"-{A['iq'][1]}", f"-{B['iq'][1]}", B["iq"][3]],
       [f"-{A['ka'][0]}", f"-{A['kb'][1]}", f"-{B['kb'][1]}",
        f"-{B['ka'][0]}"]]
# outer faces bordering the box stack per level gap [N, W, S, E]:
lvSide = [[A["bands"][2], BX["fka"][2], B["bands"][2], BX["fkb"][3]],
          [A["core"], BX["fiq"][1], B["core"], BX["fiq"][3]],
          [A["bands"][0], BX["fkb"][1], B["bands"][0], BX["fka"][0]]]
dxs, dzs = [1, -1, -1, 1], [1, 1, -1, -1]
bfP, bfQ, bfA, bfQA, bfSP, bfDG, bfF = [], [], [], [], [], [], []
for lv in range(4):
    y = ay + levels[lv]
    O = P(ax, y, zwc)
    Pm = [P(ax + dxs[m]*qw, y, zwc + dzs[m]*qw) for m in range(4)]
    Qm = [P(ax + dxs[m]*qs, y, zwc + dzs[m]*qs) for m in range(4)]
    Bc = [P(ax + [0, -1, 0, 1][m]*qco, y, zwc + [1, 0, -1, 0][m]*qco)
          for m in range(4)]
    arc, qa, sp = [], [], []
    for m in range(4):
        mp = (m+1) % 4
        c = nm("c")
        L.append(f"{c} = newl; Circle({c}) = {{{Pm[m]}, {O}, {Pm[mp]}}};")
        arc.append(c)
        c = nm("c")
        L.append(f"{c} = newl; Circle({c}) = {{{Qm[m]}, {Bc[m]}, {Qm[mp]}}};")
        qa.append(c)
        sp.append(LINE(Qm[m], Pm[m]))
    dg = [LINE(lvc[lv][m], Pm[m]) for m in range(4)]
    TC(arc + qa, NCW)
    TC(sp, NRAD)
    TC(dg, NDG, f" Using Progression {DPR}")
    fN = FACE([lvU[lv], dg[1], f"-{arc[0]}", f"-{dg[0]}"],
              corners=[lvc[lv][0], lvc[lv][1], Pm[1], Pm[0]], plane=True)
    fW = FACE([lvL[lv], dg[2], f"-{arc[1]}", f"-{dg[1]}"],
              corners=[lvc[lv][1], lvc[lv][2], Pm[2], Pm[1]], plane=True)
    fS = FACE([lvD[lv], dg[3], f"-{arc[2]}", f"-{dg[2]}"],
              corners=[lvc[lv][2], lvc[lv][3], Pm[3], Pm[2]], plane=True)
    fE = FACE([f"-{lvR[lv]}", dg[0], f"-{arc[3]}", f"-{dg[3]}"],
              corners=[lvc[lv][3], lvc[lv][0], Pm[0], Pm[3]], plane=True)
    ring = [FACE([qa[m], sp[(m+1) % 4], f"-{arc[m]}", f"-{sp[m]}"],
                 plane=True) for m in range(4)]
    corefc = FACE(qa, corners=Qm, plane=True)
    bfP.append(Pm); bfQ.append(Qm); bfA.append(arc); bfQA.append(qa)
    bfSP.append(sp); bfDG.append(dg)
    bfF.append(dict(N=fN, W=fW, S=fS, E=fE, ring=ring, core=corefc))
L.append("sWafDisk[] = {" + ", ".join([bfF[0]["core"]] + bfF[0]["ring"])
         + "};")
L.append("sExit[] += {" + ", ".join(
    [bfF[0][kk] for kk in "NWSE"] + [bfF[3][kk] for kk in "NWSE"]
    + bfF[3]["ring"] + [bfF[3]["core"]]) + "};")
for g in range(3):
    a, b = g, g + 1
    ncv = NCW if g == 1 else NRB
    grd = ("" if g == 1 else
           f" Using Progression {RATIO if g == 0 else 1/RATIO:.6g}")
    vP = [LINE(bfP[a][m], bfP[b][m]) for m in range(4)]
    vQ = [LINE(bfQ[a][m], bfQ[b][m]) for m in range(4)]
    TC(vP + vQ, ncv, grd)
    wC = [FACE([bfA[a][m], vP[(m+1) % 4], f"-{bfA[b][m]}", f"-{vP[m]}"])
          for m in range(4)]
    qC = [FACE([bfQA[a][m], vQ[(m+1) % 4], f"-{bfQA[b][m]}", f"-{vQ[m]}"])
          for m in range(4)]
    pC = [FACE([bfSP[a][m], vP[m], f"-{bfSP[b][m]}", f"-{vQ[m]}"])
          for m in range(4)]
    dC = [FACE([bfDG[a][m], vP[m], f"-{bfDG[b][m]}", neg(lvV[g][m])])
          for m in range(4)]
    VOL([bfF[a]["core"], bfF[b]["core"]] + qC, bfQ[a] + bfQ[b])
    for m in range(4):
        mp = (m+1) % 4
        VOL([bfF[a]["ring"][m], bfF[b]["ring"][m], qC[m], wC[m], pC[m],
             pC[mp]],
            [bfQ[a][m], bfQ[a][mp], bfP[a][mp], bfP[a][m],
             bfQ[b][m], bfQ[b][mp], bfP[b][mp], bfP[b][m]])
    for key, m in [("N", 0), ("W", 1), ("S", 2), ("E", 3)]:
        VOL([bfF[a][key], bfF[b][key], lvSide[g][m], wC[m], dC[m],
             dC[(m+1) % 4]],
            [lvc[a][m], lvc[a][(m+1) % 4], bfP[a][(m+1) % 4], bfP[a][m],
             lvc[b][m], lvc[b][(m+1) % 4], bfP[b][(m+1) % 4], bfP[b][m]])

L.append("sOut[] = {" + ", ".join(S[-1]["pads"] + S[-1]["bands"]
                                  + [S[-1]["core"]]) + "};")
L.append("sExit[] += {" + ", ".join(walls) + "};")
L.append("sHot[] = {" + ", ".join(walls_hot) + "};")

HEAD = ['// Full CVD train: inlet pipe (trip hump) -> diffusor -> hot zone',
        '//                 -> extended outlet.  Generated by build_full_train.py',
        'SetFactory("Built-in");',
        'lc = 0.002;',
        'sExit[] = {}; vAll[] = {}; sWafDisk[] = {}; sOut[] = {}; sIface[] = {};']
TAIL = ['',
        '// zone ids match hotzone_trip so case files port across:',
        '//   1 inlet  2 outlet  3 wallCold  4 wallHot  5 substrate',
        'Physical Surface("inlet", 1)     = {sIface[]};',
        'Physical Surface("outlet", 2)    = {sOut[]};',
        'Physical Surface("wallCold", 3)  = {sExit[]};',
        'Physical Surface("wallHot", 4)   = {sHot[]};',
        'Physical Surface("substrate", 5) = {sWafDisk[]};',
        'Physical Volume("fluid", 6)      = {vAll[]};',
        'Mesh.RecombineAll = 1;',
        'Mesh.Recombine3DAll = 1;']
open("full_train.geo", "w").write("\n".join(HEAD + L + TAIL) + "\n")
ncells = ((NCW-1)**2 + 4*(NCW-1)*(NRB-1) + 4*(NRB-1)**2)
nax = sum(n - 1 for n in seg_n)
print(f"wrote full_train.geo: {NST} stations, {cnt['v']} volumes, "
      f"{cnt['f']} faces, ~{ncells} cells/section x {nax} axial "
      f"~= {ncells*nax} hexes")
print(f"  z: {stations[-1]['z']:+.4f} (outlet) .. {stations[0]['z']:+.4f} (inlet)"
      f"   total length {stations[0]['z']-stations[-1]['z']:.4f} m")
print(f"  pipe  {L_PIPE*1e3:.0f} mm, hump crest r = {pipe_r(z2+HUMP_S)*1e3:.2f} mm "
      f"vs {rif*1e3:.2f} mm  ({HUMP_H/rif*100:.0f}% blockage)")
print(f"  diffusor {Ldiff*1e3:.0f} mm | hot zone {hzS[-1]*1e3:.0f} mm | "
      f"outlet ext {L_OUT*1e3:.0f} mm")
print(f"  wall faces: cold {len(walls)}, hot {len(walls_hot)}")
