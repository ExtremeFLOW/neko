# Visual check of the linear-IC (paper-style) convergence run, run_linic_noproj/:
#   1. viz_linic_profiles.png — final u, T, Q_T on every GLL point vs the closed forms
#   2. viz_linic_convergence.png — (a) T(x) snapshots morphing linear -> tanh front,
#      (b) max-error decay vs time, (c) projection on/off comparison (the bug-4 story)
# Run: /home/hochi/projects/matrixenv/bin/python render_linear_ic.py
import re, glob
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysemtools.io.ppymech.neksuite import preadnek
from mpi4py import MPI
comm = MPI.COMM_WORLD
d = 0.2
RUN = 'run_linic_noproj/'

def cat(fd, a, k):
    return np.concatenate([getattr(fd.elem[e], a)[k].ravel() for e in range(fd.nel)])

def errtrace(log, dt):
    """(t, eu, eT, eq) from the [validate] lines; robust to Fortran 3-digit exps."""
    rows = []
    num = re.compile(r'[-+]?[0-9]*\.?[0-9]+(?:[EeDd]?[+-][0-9]+)?')
    for line in open(log, errors='ignore'):
        if '[validate] step' not in line:
            continue
        vals = num.findall(line.replace('E+', 'e+').replace('E-', 'e-')
                               .replace('+1', 'e+1').replace('+0', 'e+0')
                               .replace('ee', 'e'))
        try:
            step = int(float(vals[0])); eu, eT, eq = (float(v) for v in vals[1:4])
        except (ValueError, IndexError):
            continue
        rows.append((step * dt, eu, eT, eq))
    return np.array(rows).T

ua = lambda xx: 0.5*(3 + np.tanh(xx/d))
qa = lambda xx: 0.5/d*(1 - np.tanh(xx/d)**2)

# ---------- figure 1: final profiles ----------
geo = preadnek(RUN + 'field0.f00000', comm)
x = cat(geo, 'pos', 0)
fl = preadnek(sorted(glob.glob(RUN + 'field0.f0*'))[-1], comm)
u, T = cat(fl, 'vel', 0), cat(fl, 'temp', 0)
ql = preadnek(sorted(glob.glob(RUN + 'qt_out0.f0*'))[-1], comm)
qtl_at = qa(x)
QTp, QTt = cat(ql, 'pres', 0), cat(ql, 'temp', 0)
QT = QTt if np.max(np.abs(QTt-qtl_at)) < np.max(np.abs(QTp-qtl_at)) else QTp

eu, eT, eq = (np.max(np.abs(v - a)) for v, a in [(u, ua(x)), (T, ua(x)), (QT, qa(x))])
xa = np.linspace(-1, 1, 600)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
for a, (val, fa, lab, err) in zip(ax, [(u, ua, 'u (velocity)', eu),
                                       (T, ua, 'T (temperature)', eT),
                                       (QT, qa, r'$Q_T$ (thermal divergence)', eq)]):
    a.plot(xa, fa(xa), 'k-', lw=2.2, label='analytic', zorder=1)
    a.scatter(x, val, s=7, c='crimson', alpha=0.4, label='Neko (all GLL pts)', zorder=2)
    a.set_title('%s\nmax err = %.2e' % (lab, err), fontsize=11)
    a.set_xlabel('x'); a.grid(alpha=0.3); a.legend(fontsize=9, loc='best')
fig.suptitle('Linear IC -> converged steady state, t = 4.0 (paper-style test)', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(RUN + 'viz_linic_profiles.png', dpi=110)
print('wrote %sviz_linic_profiles.png  max err: u=%.2e T=%.2e Q_T=%.2e' % (RUN, eu, eT, eq))

# ---------- figure 2: the convergence story ----------
fig2, ax2 = plt.subplots(1, 3, figsize=(17, 4.6))

# (a) T(x) snapshots: linear line morphing into the tanh front (nsamples 20 -> dt_out 0.2)
srt = np.argsort(x)
for i, c in zip([0, 1, 2, 4, 8, 20], plt.cm.viridis(np.linspace(0, 0.9, 6))):
    f = preadnek(RUN + 'field0.f%05d' % i, comm)
    ax2[0].plot(x[srt], cat(f, 'temp', 0)[srt], color=c, lw=1.6, label='t = %.1f' % (0.2*i))
ax2[0].plot(xa, ua(xa), 'k--', lw=1.2, label='analytic')
ax2[0].set_title('T(x): linear IC develops the tanh front')
ax2[0].set_xlabel('x'); ax2[0].grid(alpha=0.3); ax2[0].legend(fontsize=8)

# (b) error decay vs time (projection off)
t, e1, e2, e3 = errtrace(RUN + 'run.log', 1e-3)
for e, lab, c in [(e1, 'u', 'crimson'), (e2, 'T', 'royalblue'), (e3, r'$Q_T$', 'seagreen')]:
    ax2[1].semilogy(t, e, color=c, lw=1.4, label=lab)
ax2[1].axhline(7.2e-3, color='seagreen', ls=':', lw=1); ax2[1].axhline(2.0e-5, color='crimson', ls=':', lw=1)
ax2[1].set_title('max error vs t (projection off)\ndotted = hold-test discretization floor')
ax2[1].set_xlabel('t'); ax2[1].grid(alpha=0.3, which='both'); ax2[1].legend(fontsize=9)

# (c) bug 4: u-error, projection on (blows up, both dt) vs off (converges)
for log, dt, lab, c, ls in [
        ('run_linear_ic/run.log', 1e-3, 'proj 4, dt 1e-3 (NaN)', 'darkorange', '-'),
        ('run_linic_dthalf/run.log', 5e-4, 'proj 4, dt 5e-4 (NaN)', 'firebrick', '-'),
        (RUN + 'run.log', 1e-3, 'proj 0, dt 1e-3 (converges)', 'teal', '-'),
        ('run_linic_reorth/run.log', 1e-3, 'proj 4 + reorthogonalize (converges)', 'purple', '--')]:
    t, e1, _, _ = errtrace(log, dt)
    ax2[2].semilogy(t[t < 1.0], np.clip(e1[t < 1.0], None, 1e3), color=c, lw=1.5, ls=ls, label=lab)
ax2[2].set_title('u error: pressure-projection instability\nunder the density transient (bug 4)')
ax2[2].set_xlabel('t'); ax2[2].set_ylim(1e-6, 2e3)
ax2[2].grid(alpha=0.3, which='both'); ax2[2].legend(fontsize=8, loc='center right')

fig2.tight_layout()
fig2.savefig(RUN + 'viz_linic_convergence.png', dpi=110)
print('wrote %sviz_linic_convergence.png' % RUN)
