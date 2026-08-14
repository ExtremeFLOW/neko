#!/usr/bin/env python3
"""Preview of train.msh: profile, sections with block edges, quality, substrate."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import gmsh
import ogrid_cvd as og
import ogrid_smooth_cvd as ogs

NH, NV, NR, BLFRAC, RATIO = 12, 5, 3, 0.2, 1.15
WAFER_X, WAFER_R = -771.4, 68.5
HOT_X0, HOT_X1 = -950.9, -664.4


def section_full(x):
    blocks = ogs.section_smooth(x, NH, NV, NR, BLFRAC, RATIO)
    c = og.CEN(x)
    return [b + np.array([0.0, c]) for b in blocks] if c else blocks


# ---- load mesh --------------------------------------------------------------
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.open("train.msh")
nt, nc, _ = gmsh.model.mesh.getNodes()
P = nc.reshape(-1, 3)
nid = {int(t): i for i, t in enumerate(nt)}
tags = gmsh.model.mesh.getElementsByType(5)[0]
q = np.array(gmsh.model.mesh.getElementQualities(tags, "minSJ"))
_, _, en = gmsh.model.mesh.getElements(3)
conn = np.concatenate([n for n in en]).reshape(len(tags), 8)
cen = P[np.vectorize(nid.get)(conn)].mean(axis=1)
# substrate quads
sub = []
for d, t in gmsh.model.getPhysicalGroups(2):
    if gmsh.model.getPhysicalName(d, t) != "substrate":
        continue
    for e in gmsh.model.getEntitiesForPhysicalGroup(d, t):
        _, _, qn = gmsh.model.mesh.getElements(2, e)
        for arr in qn:
            sub.append(np.array(arr).reshape(-1, 4))
sub = np.vstack(sub) if sub else np.empty((0, 4), int)
gmsh.finalize()

fig = plt.figure(figsize=(16, 11))

# ---- (1) wall profile -------------------------------------------------------
ax = fig.add_subplot(3, 1, 1)
xd = np.linspace(0, og.X_OUT, 1200)
Y = np.array([og.wall(x, 0.0)[0] for x in xd])
H = np.array([og.wall(x, np.pi/2)[1] for x in xd])
C = np.array([og.CEN(x) for x in xd])
ax.plot(xd, C + H, "k-", lw=1.5, label="ceiling")
ax.plot(xd, C - H, "k-", lw=1.5, label="floor")
ax.plot(xd, Y, "b--", lw=1.2, label="half-width (plan)")
ax.plot(xd, -Y, "b--", lw=1.2)
ax.axvspan(HOT_X0, HOT_X1, color="red", alpha=.10, lw=0, label="heated section")
ax.plot([WAFER_X - WAFER_R, WAFER_X + WAFER_R],
        [np.interp(WAFER_X, xd[::-1], (C - H)[::-1])]*2, lw=5, color="#7a3d00",
        label="wafer (substrate)", solid_capstyle="butt")
for xx, lab in ((og.X_PE, "pipe end"), (og.X_CS, "diffusor exit = hot-zone entry"),
                (-1017.9, "hot-zone exit")):
    ax.axvline(xx, color="green", ls=":", lw=1.2)
    ax.text(xx, 86, lab, rotation=90, fontsize=7, va="top", ha="right", color="green")
ax.set(xlabel="train x (mm)  —  inlet right (x=0), outlet left",
       ylabel="mm", title="wall profile: pipe + hump -> diffusor -> hot zone -> extension")
ax.grid(alpha=.3); ax.legend(fontsize=8, loc="lower left", ncol=3)

# ---- (2) cross-sections with block edges -----------------------------------
XSEC = [(-100.0, "hump crest"), (-280.0, "mid-diffusor"), (-410.0, "entry"),
        (-596.0, "constriction"), (-771.4, "wafer mid")]
for k, (x, lab) in enumerate(XSEC):
    ax = fig.add_subplot(3, 5, 6 + k)
    for B in section_full(x):
        for j in range(B.shape[0]):
            ax.plot(B[j, :, 0], B[j, :, 1], "b-", lw=0.35)
        for i in range(B.shape[1]):
            ax.plot(B[:, i, 0], B[:, i, 1], "b-", lw=0.35)
    ax.set_aspect("equal")
    ax.set_title(f"{lab}\nx={x:.0f}", fontsize=8)
    ax.tick_params(labelsize=6)

# ---- (3) quality + substrate footprint -------------------------------------
ax = fig.add_subplot(3, 2, 5)
s = ax.scatter(cen[:, 0], cen[:, 2], c=q, s=1.0, cmap="RdYlGn", vmin=0.3, vmax=1.0, lw=0)
fig.colorbar(s, ax=ax, pad=.01, aspect=12, label="min scaled Jacobian")
bad = q < 0.6
ax.scatter(cen[bad, 0], cen[bad, 2], s=5, facecolors="none", edgecolors="k",
           lw=.3, label=f"minSJ < 0.6  ({bad.sum()})")
ax.set(xlabel="x (mm)", ylabel="z (mm)", title=f"element quality (min {q.min():.3f})")
ax.legend(fontsize=7); ax.grid(alpha=.3)

ax = fig.add_subplot(3, 2, 6)
for quad in sub:
    p = P[[nid[int(n)] for n in quad]]
    ax.fill(p[:, 0], p[:, 1], facecolor="#c9862b", edgecolor="k", lw=.25, alpha=.85)
th = np.linspace(0, 2*np.pi, 400)
ax.plot(WAFER_X + WAFER_R*np.cos(th), WAFER_R*np.sin(th), "r-", lw=1.6,
        label="exact wafer disc r=68.5")
ax.set_aspect("equal")
ax.set(xlabel="x (mm)", ylabel="y (mm)",
       title=f"substrate zone: {len(sub)} floor faces, staircase rim (+1.3% area)")
ax.legend(fontsize=8); ax.grid(alpha=.3)

fig.suptitle(f"train.msh — one-piece O-grid sweep, {len(tags):,} hexes, "
             f"minSJ {q.min():.3f}  (frame attempt was 0.148)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.965])
fig.savefig("train_preview.png", dpi=110)
print("wrote train_preview.png")
