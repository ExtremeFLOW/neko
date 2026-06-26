# Visual validation: Neko low-Mach lowMach_test vs the analytic (Nek5000) solution.
#   u = T = 0.5(3 + tanh(x/d)),  QTL = 0.5/d (1 - tanh^2(x/d)),  d = 0.2.
# Neko fld geometry is written only in f00000 (var[0]=3); later snapshots carry
# fields only -> read coords from f00000, fields from the latest file (same mesh
# ordering). Overlays every GLL point of Neko (dots) on the closed-form curves.
# Run: /home/hochi/projects/matrixenv/bin/python render_validation.py
import numpy as np, glob, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysemtools.io.ppymech.neksuite import preadnek
from mpi4py import MPI
comm = MPI.COMM_WORLD
d = 0.2

def cat(fd, a, k):
    return np.concatenate([getattr(fd.elem[e], a)[k].ravel() for e in range(fd.nel)])

# geometry from f00000
geo = preadnek('field0.f00000', comm)
x = cat(geo, 'pos', 0); y = cat(geo, 'pos', 1); z = cat(geo, 'pos', 2)
# fields from the latest snapshots
fl = preadnek(sorted(glob.glob('field0.f0*'))[-1], comm)
u = cat(fl, 'vel', 0); T = cat(fl, 'temp', 0)
ql = preadnek(sorted(glob.glob('qt_out0.f0*'))[-1], comm)
# Q_T sits in qt_out's temp or pres slot (field_writer packed [temperature, Q_T])
qtl_at = 0.5/d*(1 - np.tanh(x/d)**2)
QTp, QTt = cat(ql, 'pres', 0), cat(ql, 'temp', 0)
QT = QTt if np.max(np.abs(QTt-qtl_at)) < np.max(np.abs(QTp-qtl_at)) else QTp

ua = lambda xx: 0.5*(3 + np.tanh(xx/d))
qa = lambda xx: 0.5/d*(1 - np.tanh(xx/d)**2)
eu = np.max(np.abs(u-ua(x))); eT = np.max(np.abs(T-ua(x))); eq = np.max(np.abs(QT-qa(x)))
xa = np.linspace(-1, 1, 600)

# --- profile overlay ---
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
panels = [(u, ua, 'u (velocity)', eu), (T, ua, 'T (temperature)', eT),
          (QT, qa, r'$Q_T$ (thermal divergence)', eq)]
for a, (val, fa, lab, err) in zip(ax, panels):
    a.plot(xa, fa(xa), 'k-', lw=2.2, label='analytic (Nek5000)', zorder=1)
    a.scatter(x, val, s=7, c='crimson', alpha=0.4, label='Neko (all GLL pts)', zorder=2)
    a.set_title('%s\nmax err = %.2e' % (lab, err), fontsize=11)
    a.set_xlabel('x'); a.grid(alpha=0.3); a.legend(fontsize=9, loc='best')
fig.suptitle('Neko low-Mach vs analytic lowMach_test   t=2.0 (2000 steps), %s'
             % sorted(glob.glob('field0.f0*'))[-1], fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('viz_profiles.png', dpi=110)
print('wrote viz_profiles.png   max err: u=%.2e T=%.2e Q_T=%.2e' % (eu, eT, eq))

# --- 2D contours on one z-plane (solution is x-only -> vertical bands) ---
m = np.isclose(z, z.min())
fig2, ax2 = plt.subplots(1, 3, figsize=(15, 3.6))
for a, (val, lab, cm) in zip(ax2, [(u[m], 'u', 'coolwarm'), (T[m], 'T', 'coolwarm'),
                                    (QT[m], r'$Q_T$', 'viridis')]):
    sc = a.tricontourf(x[m], y[m], val, levels=30, cmap=cm)
    a.set_title('Neko ' + lab); a.set_xlabel('x'); a.set_ylabel('y')
    fig2.colorbar(sc, ax=a)
fig2.suptitle('Neko fields, z-plane (x-only solution -> vertical bands)')
fig2.tight_layout(rect=[0, 0, 1, 0.93])
fig2.savefig('viz_contour.png', dpi=110)
print('wrote viz_contour.png')
