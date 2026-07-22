# Where are the temperature under/overshoots? Latest turbulent hump frame:
# scan the FULL 3D field for T < 1 (undershoot) and T > 2 (overshoot), project
# the violating points onto (x,y) over a mid-z T context map + a wake zoom.
# Run: /home/hochi/projects/matrixenv/bin/python render_gibbs_map.py
import glob
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from pysemtools.io.ppymech.neksuite import preadnek
from mpi4py import MPI
comm = MPI.COMM_WORLD
RUN = '../heated_hump_res/Old_53K_turbulent/'

def cat(fd, a, k):
    return np.concatenate([getattr(fd.elem[e], a)[k].ravel() for e in range(fd.nel)])

geo = preadnek(RUN + 'field0.f00000', comm)
x, y, z = (cat(geo, 'pos', k) for k in range(3))
last = sorted(glob.glob(RUN + 'field0.f0*'))[-1]
fd = preadnek(last, comm)
T = cat(fd, 'temp', 0)

under = np.maximum(1.0 - T, 0.0)
over = np.maximum(T - 2.0, 0.0)
mu_ = under > 0.02
mo_ = over > 0.02
iu, io = np.argmax(under), np.argmax(over)
print(f'points: {T.size};  T range [{T.min():.3f}, {T.max():.3f}]')
print(f'undershoot >0.02: {mu_.sum()} pts ({100*mu_.sum()/T.size:.3f}%), worst {under[iu]:.3f} at (x,y,z)=({x[iu]:.2f},{y[iu]:.2f},{z[iu]:.2f})')
print(f'overshoot  >0.02: {mo_.sum()} pts ({100*mo_.sum()/T.size:.3f}%), worst {over[io]:.3f} at (x,y,z)=({x[io]:.2f},{y[io]:.2f},{z[io]:.2f})')
print(f'clamp-active (T<0.8): {np.sum(T < 0.8)} pts')

# mid-z context slice with hump masked
zlev = np.unique(np.round(z, 8))
mz = np.isclose(z, zlev[np.argmin(np.abs(zlev - zlev.mean()))])
xi, yi = x[mz], y[mz]
def ywall_fn(xq):
    yw = np.zeros_like(xq)
    m = (xq > 4.0) & (xq < 8.0)
    yw[m] = 0.5 * (1.0 + np.cos(np.pi * (xq[m] - 6.0) / 2.0))
    return yw
tri = mtri.Triangulation(xi, yi)
xc = xi[tri.triangles].mean(axis=1); yc = yi[tri.triangles].mean(axis=1)
tri.set_mask(yc < ywall_fn(xc) - 1e-3)
xu = np.linspace(xi.min(), xi.max(), 800)

fig, ax = plt.subplots(2, 1, figsize=(14, 8))
for a, (x0, x1, y0, y1, ttl) in zip(ax, [(0, 25, 0, 3, 'full domain'),
                                          (8, 20, 0, 1.2, 'wake / heated-wall zoom')]):
    a.tricontourf(tri, np.clip(T[mz], 1, 2), levels=np.linspace(1, 2, 31),
                  cmap='Greys', alpha=0.75)
    a.fill_between(xu, 0, ywall_fn(xu), color='0.3', zorder=3)
    su = a.scatter(x[mu_], y[mu_], s=4 + 300*under[mu_], c='#1f77b4', alpha=0.5,
                   label=f'T < 1 (undershoot), n={mu_.sum()}')
    so = a.scatter(x[mo_], y[mo_], s=4 + 300*over[mo_], c='#d62728', alpha=0.5,
                   label=f'T > 2 (overshoot), n={mo_.sum()}')
    a.set_xlim(x0, x1); a.set_ylim(y0, y1)
    a.set_xlabel('x'); a.set_ylabel('y'); a.set_title(ttl, fontsize=11)
    a.set_aspect('equal')
ax[0].legend(fontsize=9, loc='upper left')
ax[0].axvline(11, color='crimson', lw=0.8, ls='--')
ax[0].annotate('heated wall starts', xy=(11.1, 2.6), fontsize=9, color='crimson')
fig.suptitle(f'Gibbs violations, all z-planes projected onto (x,y) — {last.split("/")[-1]} (t = 44)\n'
             'grey = mid-z temperature; dot size ∝ violation magnitude', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(RUN + 'viz_gibbs_map.png', dpi=110)
print('wrote', RUN + 'viz_gibbs_map.png')
