# Visual check of the overnight gravity soak (Rayleigh-Benard-Poiseuille channel):
#   viz_overnight_state.png     — latest snapshot: hot-wall streaks, cross-plane
#                                 T plumes + v rolls, plane-averaged profiles
#   viz_overnight_evolution.png — T cross-plane at t ~ 5 / 25 / 50 / latest
# Run: /home/hochi/projects/matrixenv/bin/python render_overnight.py
import glob
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysemtools.io.ppymech.neksuite import preadnek
from mpi4py import MPI
comm = MPI.COMM_WORLD
RUN = 'overnight_gravity/'
DT_OUT = 5.0

def cat(fd, a, k):
    return np.concatenate([getattr(fd.elem[e], a)[k].ravel() for e in range(fd.nel)])

geo = preadnek(RUN + 'field0.f00000', comm)
x, y, z = (cat(geo, 'pos', k) for k in range(3))

files = sorted(glob.glob(RUN + 'field0.f0*'))
last = files[-2]                       # -2: newest may still be mid-write
ilast = int(last[-5:])
fl = preadnek(last, comm)
u, v, T = cat(fl, 'vel', 0), cat(fl, 'vel', 1), cat(fl, 'temp', 0)

# nearest existing GLL levels to the requested slice positions
def nearest_level(coord, target):
    lv = np.unique(np.round(coord, 8))
    return lv[np.argmin(np.abs(lv - target))]

y_streak = nearest_level(y, 0.10)                    # ~y+ 18 off the hot wall
x_cut = nearest_level(x, np.pi)
my = np.isclose(y, y_streak)
mx = np.isclose(x, x_cut)

# plane-averaged profiles over (x, z)
yr = np.round(y, 8)
ylev, inv = np.unique(yr, return_inverse=True)
cnt = np.bincount(inv)
u_mean = np.bincount(inv, u) / cnt
T_mean = np.bincount(inv, T) / cnt

fig = plt.figure(figsize=(16, 8.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1])

ax1 = fig.add_subplot(gs[0, :2])
s1 = ax1.tricontourf(x[my], z[my], u[my], levels=40, cmap='viridis')
ax1.set_title(f'streamwise u, wall-parallel plane y = {y_streak:.3f} (near hot wall) — streaks')
ax1.set_xlabel('x'); ax1.set_ylabel('z'); ax1.set_aspect('equal')
fig.colorbar(s1, ax=ax1, label='u')

ax2 = fig.add_subplot(gs[0, 2])
s2 = ax2.tricontourf(z[mx], y[mx], T[mx], levels=40, cmap='RdBu_r',
                     vmin=1.0, vmax=4.0)
ax2.set_title(f'T, cross-plane x = {x_cut:.2f} — plumes')
ax2.set_xlabel('z'); ax2.set_ylabel('y'); ax2.set_aspect('equal')
fig.colorbar(s2, ax=ax2, label='T')

ax3 = fig.add_subplot(gs[1, 0])
vmax = np.percentile(np.abs(v[mx]), 99)
s3 = ax3.tricontourf(z[mx], y[mx], v[mx], levels=40, cmap='RdBu_r',
                     vmin=-vmax, vmax=vmax)
ax3.set_title('wall-normal v, same cross-plane — rolls')
ax3.set_xlabel('z'); ax3.set_ylabel('y'); ax3.set_aspect('equal')
fig.colorbar(s3, ax=ax3, label='v')

ax4 = fig.add_subplot(gs[1, 1])
ax4.plot(u_mean, ylev, lw=2, color='#1f77b4')
ax4.set_title('plane-averaged  ⟨u⟩(y)')
ax4.set_xlabel('⟨u⟩'); ax4.set_ylabel('y'); ax4.grid(alpha=0.3)

ax5 = fig.add_subplot(gs[1, 2])
ax5.plot(T_mean, ylev, lw=2, color='#d62728', label='⟨T⟩, now')
ax5.plot([4, 1], [0, 2], 'k--', lw=1.2, label='initial linear')
ax5.set_title('plane-averaged  ⟨T⟩(y)')
ax5.set_xlabel('⟨T⟩'); ax5.set_ylabel('y'); ax5.grid(alpha=0.3)
ax5.legend(fontsize=9)

fig.suptitle(f'Overnight gravity soak (g_nd = 0.1, unstable): snapshot {last}  t ≈ {ilast*DT_OUT:.0f}',
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(RUN + 'viz_overnight_state.png', dpi=110)
print('wrote', RUN + 'viz_overnight_state.png',
      ' T range [%.2f, %.2f]  u range [%.2f, %.2f]' % (T.min(), T.max(), u.min(), u.max()))

# --- evolution of the cross-plane temperature ---
picks = [i for i in [1, 5, 10] if i < ilast] + [ilast]
fig2, axs = plt.subplots(1, len(picks), figsize=(4.2*len(picks), 4.2), sharey=True)
for a, i in zip(np.atleast_1d(axs), picks):
    fd = preadnek(RUN + 'field0.f%05d' % i, comm)
    Ti = cat(fd, 'temp', 0)
    sc = a.tricontourf(z[mx], y[mx], Ti[mx], levels=40, cmap='RdBu_r', vmin=1, vmax=4)
    a.set_title(f't ≈ {i*DT_OUT:.0f}')
    a.set_xlabel('z'); a.set_aspect('equal')
np.atleast_1d(axs)[0].set_ylabel('y')
fig2.colorbar(sc, ax=list(np.atleast_1d(axs)), label='T', shrink=0.85)
fig2.suptitle('T cross-plane: development of buoyant convection (hot wall below)')
fig2.savefig(RUN + 'viz_overnight_evolution.png', dpi=110, bbox_inches='tight')
print('wrote', RUN + 'viz_overnight_evolution.png')
