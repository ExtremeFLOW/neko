#!/usr/bin/env python3
"""Element-aware y=0 mid-plane slice of the fulltrain run (side view).

Frame: x streamwise (inlet x=0, flow -x), y spanwise, z vertical.
The y=0 plane passes through the butterfly's N/S sector boundaries, so it is
covered by element faces -- same per-element slicing as render_slice.py, with
the sliced coordinate being y (index 1) and the plot plane (x, z).

usage: slice_train.py [field0.fNNNNN]   (default: newest)
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
geo = preadnek("field0.f00000", comm)
fd = preadnek(fn, comm)
t = getattr(fd, "time", float("nan"))

Xs, Zs, TRI, plan = [], [], [], []
base = 0
for e in range(geo.nel):
    ye = geo.elem[e].pos[1]
    best = None
    for ax_ in range(3):
        others = tuple(k for k in range(3) if k != ax_)
        err = np.abs(ye).max(axis=others)              # plane y = 0
        i = int(np.argmin(err))
        if best is None or err[i] < best[0]:
            best = (err[i], ax_, i)
    emin, ax_, i = best
    if emin > 1e-7:
        continue
    xe = np.take(geo.elem[e].pos[0], i, axis=ax_)
    ze = np.take(geo.elem[e].pos[2], i, axis=ax_)
    n1, n2 = xe.shape
    plan.append((e, ax_, i))
    Xs.append(xe.ravel()); Zs.append(ze.ravel())
    idx = base + np.arange(n1*n2).reshape(n1, n2)
    a = idx[:-1, :-1].ravel(); b = idx[1:, :-1].ravel()
    c = idx[1:, 1:].ravel();   d = idx[:-1, 1:].ravel()
    TRI.append(np.column_stack([a, b, c]))
    TRI.append(np.column_stack([a, c, d]))
    base += n1*n2
X = np.concatenate(Xs); Z = np.concatenate(Zs)
T = mtri.Triangulation(X, Z, np.concatenate(TRI))
print(f"{fn}: t = {t:.4f}; slice {len(plan)} elements, {len(X):,} points")


def grab(getter):
    return np.concatenate([np.take(getter(fd.elem[e]), i, axis=ax_).ravel()
                           for (e, ax_, i) in plan])


U = grab(lambda el: el.vel[0])       # streamwise (negative = downstream)
W = grab(lambda el: el.vel[2])       # vertical -- ~0 wherever laminar
Tt = grab(lambda el: el.temp[0])

fig, ax = plt.subplots(3, 1, figsize=(17, 8.2))
panels = [(-U, "streamwise  -u  (positive = downstream)", "viridis", None),
          (W, "vertical  w   (zero in laminar flow -- turbulence indicator)",
           "RdBu_r", "sym"),
          (Tt, "temperature  (ramping walls)", "inferno", (1.0, 6.0))]
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
    a.set(ylabel="z (m)", title=lab, aspect=6.0)
    a.set_xlim(X.min(), X.max())
ax[-1].set_xlabel("x (m)  —  inlet right (x=0), flow right-to-left; orange=hump, "
                  "green=pipe end/diffusor exit, red=heated section")
fig.suptitle(f"fulltrain (Re_pipe 4000)  —  y=0 mid-plane,  t = {t:.3f}", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = f"slice_{os.path.basename(fn).replace('.', '_')}.png"
fig.savefig(out, dpi=105)
print("wrote", out)
