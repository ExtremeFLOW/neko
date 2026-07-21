# Instantaneous velocity contours from the latest muT snapshot: three cuts of u
# (wall-parallel near the hot wall, streamwise x-y, cross-stream z-y) plus v in
# the cross-plane to show the buoyant up/downdrafts driving it.
# Run: /home/hochi/projects/matrixenv/bin/python render_vel.py
import glob
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysemtools.io.ppymech.neksuite import preadnek
from mpi4py import MPI
comm = MPI.COMM_WORLD
RUN = 'overnight_gravity_muT/'
DT_OUT = 5.0

def cat(fd, a, k):
    return np.concatenate([getattr(fd.elem[e], a)[k].ravel() for e in range(fd.nel)])

geo = preadnek(RUN + 'field0.f00000', comm)
x, y, z = (cat(geo, 'pos', k) for k in range(3))
last = sorted(glob.glob(RUN + 'field0.f0*'))[-2]
fd = preadnek(last, comm)
u, v = cat(fd, 'vel', 0), cat(fd, 'vel', 1)

def lev(coord, target):
    lv = np.unique(np.round(coord, 8))
    return lv[np.argmin(np.abs(lv - target))]

my = np.isclose(y, lev(y, 0.10))       # wall-parallel, ~y+ 18 off hot wall
mz = np.isclose(z, lev(z, np.pi/2))    # streamwise x-y cut
mx = np.isclose(x, lev(x, np.pi))      # cross-stream z-y cut

fig, ax = plt.subplots(2, 2, figsize=(15.5, 8))
t_now = int(last[-5:]) * DT_OUT

s = ax[0,0].tricontourf(x[my], z[my], u[my], levels=40, cmap='viridis')
ax[0,0].set_title('u — wall-parallel plane y = 0.10 (near hot wall)')
ax[0,0].set_xlabel('x'); ax[0,0].set_ylabel('z'); ax[0,0].set_aspect('equal')
fig.colorbar(s, ax=ax[0,0], label='u')

s = ax[0,1].tricontourf(x[mz], y[mz], u[mz], levels=40, cmap='viridis')
ax[0,1].set_title('u — streamwise plane z = π/2')
ax[0,1].set_xlabel('x'); ax[0,1].set_ylabel('y'); ax[0,1].set_aspect('equal')
fig.colorbar(s, ax=ax[0,1], label='u')

s = ax[1,0].tricontourf(z[mx], y[mx], u[mx], levels=40, cmap='viridis')
ax[1,0].set_title('u — cross-stream plane x = π')
ax[1,0].set_xlabel('z'); ax[1,0].set_ylabel('y'); ax[1,0].set_aspect('equal')
fig.colorbar(s, ax=ax[1,0], label='u')

vmax = np.percentile(np.abs(v[mx]), 99)
s = ax[1,1].tricontourf(z[mx], y[mx], v[mx], levels=40, cmap='RdBu_r',
                        vmin=-vmax, vmax=vmax)
ax[1,1].set_title('v — same cross-plane (red up / blue down)')
ax[1,1].set_xlabel('z'); ax[1,1].set_ylabel('y'); ax[1,1].set_aspect('equal')
fig.colorbar(s, ax=ax[1,1], label='v')

fig.suptitle(f'Instantaneous velocity, variable-viscosity soak — {last}  (t ≈ {t_now:.0f})',
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(RUN + 'viz_velocity.png', dpi=110)
print('wrote', RUN + 'viz_velocity.png',
      ' u [%.2f, %.2f]  v [%.2f, %.2f]' % (u.min(), u.max(), v.min(), v.max()))
