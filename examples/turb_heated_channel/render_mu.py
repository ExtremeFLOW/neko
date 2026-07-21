# mu(T) field evolution for the variable-viscosity soak (overnight_gravity_muT/).
#   Top row: mu in the (z,y) cross-plane at successive times (sequential viridis).
#   Bottom: mu-vs-T scatter of ALL GLL points vs the analytic power law
#           mu = mu_ref (T/T_ref)^0.7  — collapse onto the curve verifies the
#           built-in property model; plus plane-averaged <mu>(y) profiles.
# props0 field_writer slot packing (known gotcha): pres=rho, vel0=mu_tot,
# vel1=lambda, vel2=T.
# Run: /home/hochi/projects/matrixenv/bin/python render_mu.py
import glob
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysemtools.io.ppymech.neksuite import preadnek
from mpi4py import MPI
comm = MPI.COMM_WORLD
RUN = 'overnight_gravity_muT/'
DT_OUT = 5.0
MU_REF, T_REF, EXPO = 3.571428571e-4, 1.0, 0.7

def cat(fd, a, k):
    return np.concatenate([getattr(fd.elem[e], a)[k].ravel() for e in range(fd.nel)])

geo = preadnek(RUN + 'field0.f00000', comm)
x, y, z = (cat(geo, 'pos', k) for k in range(3))
files = sorted(glob.glob(RUN + 'props0.f0*'))
ilast = int(files[-2][-5:])                       # -2: newest may be mid-write
picks = sorted(set([0, 1, max(1, ilast//2), ilast]))

xlev = np.unique(np.round(x, 8))
mx = np.isclose(x, xlev[np.argmin(np.abs(xlev - np.pi))])

fig = plt.figure(figsize=(4.3*len(picks), 7.6))
gs = fig.add_gridspec(2, 2*len(picks), height_ratios=[1, 1.15], hspace=0.35)
axs_top = [fig.add_subplot(gs[0, 2*k:2*k+2]) for k in range(len(picks))]
mu_lo, mu_hi = MU_REF*0.95, MU_REF*(4.0/T_REF)**EXPO*1.02
snap = {}
for a, i in zip(axs_top, picks):
    fd = preadnek(RUN + 'props0.f%05d' % i, comm)
    mu, T = cat(fd, 'vel', 0), cat(fd, 'vel', 2)
    snap[i] = (mu, T)
    sc = a.tricontourf(z[mx], y[mx], mu[mx], levels=40, cmap='viridis',
                       vmin=mu_lo, vmax=mu_hi)
    a.set_title(f't ≈ {i*DT_OUT:.0f}')
    a.set_xlabel('z'); a.set_aspect('equal')
axs_top[0].set_ylabel('y')
fig.colorbar(sc, ax=axs_top, label=r'$\mu$', shrink=0.9)

# mu vs T collapse (latest snapshot, all GLL points)
mu, T = snap[ilast]
axb = fig.add_subplot(gs[1, :len(picks)])
axb.scatter(T[::7], mu[::7], s=3, c='#1f77b4', alpha=0.25,
            label='GLL points (every 7th)')
Ta = np.linspace(1, 4, 200)
axb.plot(Ta, MU_REF*(Ta/T_REF)**EXPO, 'k-', lw=1.8,
         label=r'$\mu_{ref}(T/T_{ref})^{0.7}$')
dev = np.max(np.abs(mu - MU_REF*(T/T_REF)**EXPO))
axb.set_title(f'property-model collapse (t ≈ {ilast*DT_OUT:.0f}), max dev = {dev:.1e}')
axb.set_xlabel('T'); axb.set_ylabel(r'$\mu$'); axb.grid(alpha=0.3)
axb.legend(fontsize=9)

# plane-averaged <mu>(y) at the picked times
axc = fig.add_subplot(gs[1, len(picks):])
yr = np.round(y, 8)
ylev, inv = np.unique(yr, return_inverse=True)
cnt = np.bincount(inv)
cols = plt.cm.viridis(np.linspace(0, 0.9, len(picks)))
for i, c in zip(picks, cols):
    mui = snap[i][0]
    axc.plot(np.bincount(inv, mui)/cnt, ylev, lw=1.8, color=c,
             label=f't ≈ {i*DT_OUT:.0f}')
axc.set_title(r'plane-averaged  $\langle\mu\rangle(y)$')
axc.set_xlabel(r'$\langle\mu\rangle$'); axc.set_ylabel('y')
axc.grid(alpha=0.3); axc.legend(fontsize=9)

fig.suptitle('Variable viscosity soak: '
             r'$\mu(T) = \mu_{ref}(T/T_{ref})^{0.7}$'
             f'   (mu_ref = {MU_REF:.3e}, hot/cold ratio {(4.0)**EXPO:.2f})',
             fontsize=12)
fig.savefig(RUN + 'viz_mu_evolution.png', dpi=110, bbox_inches='tight')
print('wrote', RUN + 'viz_mu_evolution.png',
      ' mu range [%.3e, %.3e], max power-law dev %.2e' % (mu.min(), mu.max(), dev))
