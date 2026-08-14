#!/usr/bin/env python3
"""z-normal (top view) slice of the fulltrain run at mid-height z=0.

Frame: x streamwise (inlet x=0, flow -x), y spanwise, z vertical.
Unlike the y=0 mid-plane (which coincides with butterfly sector faces),
z=0 need not be an exact element face everywhere, so for every element
whose z-range straddles the plane we take the GLL layer closest to it
(accepted if within a quarter of the element's z-extent -- visually
indistinguishable, no interpolation needed).

usage: slice_train_z.py [field0.fNNNNN] [z0]   (default: newest, z0=0)
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek

comm = MPI.COMM_WORLD
X_PE, X_CS, X_HZEND = -0.150, -0.410, -1.0179
HOT0, HOT1 = -0.9509, -0.6644
HUMP = -0.100

_latest = lambda: sorted(f for f in glob.glob("field0.f0*")
                         if "lock" not in f and not glob.glob(f + "-*.lock"))[-1]
fn = sys.argv[1] if len(sys.argv) > 1 else _latest()
z0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
geo = preadnek("field0.f00000", comm)
fd = preadnek(fn, comm)
t = getattr(fd, "time", float("nan"))

Xs, Ys, TRI, plan = [], [], [], []
base = 0
nskip = 0
for e in range(geo.nel):
    ze = geo.elem[e].pos[2]
    zmin, zmax = ze.min(), ze.max()
    if not (zmin - 1e-9 <= z0 <= zmax + 1e-9):
        continue
    best = None
    for ax_ in range(3):
        for i in range(ze.shape[ax_]):
            dev = np.abs(np.take(ze, i, axis=ax_) - z0).max()
            if best is None or dev < best[0]:
                best = (dev, ax_, i)
    dev, ax_, i = best
    if dev > 0.25*max(zmax - zmin, 1e-12):
        nskip += 1
        continue
    xe = np.take(geo.elem[e].pos[0], i, axis=ax_)
    ye = np.take(geo.elem[e].pos[1], i, axis=ax_)
    n1, n2 = xe.shape
    plan.append((e, ax_, i))
    Xs.append(xe.ravel()); Ys.append(ye.ravel())
    idx = base + np.arange(n1*n2).reshape(n1, n2)
    a = idx[:-1, :-1].ravel(); b = idx[1:, :-1].ravel()
    c = idx[1:, 1:].ravel();   d = idx[:-1, 1:].ravel()
    TRI.append(np.column_stack([a, b, c]))
    TRI.append(np.column_stack([a, c, d]))
    base += n1*n2
X = np.concatenate(Xs); Y = np.concatenate(Ys)
T = mtri.Triangulation(X, Y, np.concatenate(TRI))
print(f"{fn}: t = {t:.4f}; z={z0:g} slice {len(plan)} elements "
      f"({nskip} skipped), {len(X):,} points")


def grab(getter):
    return np.concatenate([np.take(getter(fd.elem[e]), i, axis=ax_).ravel()
                           for (e, ax_, i) in plan])


U = grab(lambda el: el.vel[0])       # streamwise (negative = downstream)
V = grab(lambda el: el.vel[1])       # spanwise -- in-plane turbulence indicator
Tt = grab(lambda el: el.temp[0])

fig, ax = plt.subplots(3, 1, figsize=(17, 7.4))
panels = [(-U, "streamwise  -u  (positive = downstream)", "viridis", None),
          (V, "spanwise  v   (zero in laminar flow -- turbulence indicator)",
           "RdBu_r", "sym"),
          (Tt, "temperature", "inferno", (1.0, 6.0))]
for a, (q, lab, cm, lim) in zip(ax, panels):
    if lim == "sym":
        m = max(np.percentile(np.abs(q), 99.8), 1e-12)
        lo, hi = -m, m
    elif lim:
        lo, hi = lim
    else:
        lo, hi = np.percentile(q, 0.3), np.percentile(q, 99.7)
    s = a.tripcolor(T, np.clip(q, lo, hi), shading="gouraud", cmap=cm,
                    vmin=lo, vmax=hi)
    fig.colorbar(s, ax=a, pad=.006, aspect=9)
    for xx, c, ls in ((HUMP, "darkorange", "-"), (X_PE, "green", ":"),
                      (X_CS, "green", "--"), (HOT1, "red", ":"),
                      (HOT0, "red", ":"), (X_HZEND, "purple", ":")):
        a.axvline(xx, color=c, ls=ls, lw=1.1, alpha=.8)
    a.set(ylabel="y (m)", title=lab)
    a.set_xlim(X.min(), X.max())
ax[-1].set_xlabel("x (m)  —  inlet right (x=0), flow right-to-left; "
                  "orange=hump, green=pipe end/diffusor exit, red=heated section")
fig.suptitle(f"fulltrain (Re_pipe 4000)  —  z={z0:g} top view,  t = {t:.3f}",
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = f"slicez_{os.path.basename(fn).replace('.', '_')}.png"
fig.savefig(out, dpi=105)
print("wrote", out)
