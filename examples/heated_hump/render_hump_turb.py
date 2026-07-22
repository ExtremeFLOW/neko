# Heated hump (53K, order 3) — Dong-restart run frames scp'd from Sigma.
# Renders the mid-z slice: streamwise u, temperature, wall-normal v, and an
# outlet-region zoom at the latest frame (t past the Re 3000, order 6).
# Run: /home/hochi/projects/matrixenv/bin/python render_hump_res.py
import glob
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysemtools.io.ppymech.neksuite import preadnek
from mpi4py import MPI
comm = MPI.COMM_WORLD
RUN = '../heated_hump_res/Old_53K_turbulent/'

def cat(fd, a, k):
    return np.concatenate([getattr(fd.elem[e], a)[k].ravel() for e in range(fd.nel)])

geo = preadnek(RUN + 'field0.f00000', comm)
x, y, z = (cat(geo, 'pos', k) for k in range(3))
files = sorted(glob.glob(RUN + 'field0.f0*'))
last = files[-1]
fd = preadnek(last, comm)
t_now = getattr(fd, 'time', None)
u, v, T = cat(fd, 'vel', 0), cat(fd, 'vel', 1), cat(fd, 'temp', 0)

zlev = np.unique(np.round(z, 8))
mz = np.isclose(z, zlev[np.argmin(np.abs(zlev - zlev.mean()))])

fig, ax = plt.subplots(3, 1, figsize=(16, 10))
s = ax[0].tricontourf(x[mz], y[mz], u[mz], levels=48, cmap='viridis')
ax[0].set_title('streamwise u — mid-z slice')
fig.colorbar(s, ax=ax[0], label='u')
s = ax[1].tricontourf(x[mz], y[mz], T[mz], levels=48, cmap='RdBu_r', vmin=1, vmax=2)
ax[1].set_title('temperature (hot wall x > 11 at T = 2)')
fig.colorbar(s, ax=ax[1], label='T')
vmax = np.percentile(np.abs(v[mz]), 99.5)
s = ax[2].tricontourf(x[mz], y[mz], v[mz], levels=48, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
ax[2].set_title('wall-normal v — wake structures')
fig.colorbar(s, ax=ax[2], label='v')
for a in ax:
    a.set_xlabel('x'); a.set_ylabel('y'); a.set_aspect('equal')
ttl = f't = {t_now:.2f}' if t_now is not None else last.split('/')[-1]
fig.suptitle(f'Heated hump, Re 3000 turbulent, Dong outflow — {ttl}  (past the Re 3000, order 6)', fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(RUN + 'viz_hump_state.png', dpi=110)
print('wrote', RUN + 'viz_hump_state.png',
      ' u [%.2f, %.2f] v [%.2f, %.2f] T [%.3f, %.3f]' % (u.min(), u.max(), v.min(), v.max(), T.min(), T.max()))

# outlet-region zoom: u and v near the exit (x > 21) — the old crime scene
mo = mz & (x > 21.0)
fig2, ax2 = plt.subplots(1, 2, figsize=(13, 4.2))
s = ax2[0].tricontourf(x[mo], y[mo], u[mo], levels=40, cmap='viridis')
ax2[0].set_title('u near outlet (x > 21)')
fig2.colorbar(s, ax=ax2[0], label='u')
vm = np.percentile(np.abs(v[mo]), 99.5)
s = ax2[1].tricontourf(x[mo], y[mo], v[mo], levels=40, cmap='RdBu_r', vmin=-vm, vmax=vm)
ax2[1].set_title('v near outlet — Dong keeping backflow tame')
fig2.colorbar(s, ax=ax2[1], label='v')
for a in ax2:
    a.set_xlabel('x'); a.set_ylabel('y'); a.set_aspect('equal')
fig2.suptitle(f'Outlet region, {ttl}')
fig2.tight_layout(rect=[0, 0, 1, 0.93])
fig2.savefig(RUN + 'viz_hump_outlet.png', dpi=110)
print('wrote', RUN + 'viz_hump_outlet.png')
