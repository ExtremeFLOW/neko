# h1=1.0 verification figure: the same linear-IC + reorth case run with the
# old cref=(min+max)/2 preconditioner build (run_linic_reorth, June-26 lib) vs
# the h1=1.0 Nek-parity build (run_linic_h1verify).
#   left  — u/T error decay overlay: curves must coincide (preconditioner
#           changes the path to the answer, never the answer)
#   right — pressure iterations per step: efficiency comparison
# Run: /home/hochi/projects/matrixenv/bin/python render_h1.py
import re
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def errtrace(log, dt=1e-3):
    rows = []
    for line in open(log, errors='ignore'):
        if '[validate] step' in line:
            v = re.findall(r'[-+]?[0-9]*\.?[0-9]+(?:[Ee][+-][0-9]+)?', line)
            try:
                rows.append((int(float(v[0]))*dt, float(v[1]), float(v[2])))
            except (ValueError, IndexError):
                pass
    return np.array(rows).T

def iters(log, dt=1e-3):
    rows = []
    for line in open(log, errors='ignore'):
        if '| Pressure' in line and 'Projection' not in line:
            f = line.split()
            try:
                rows.append((int(f[0])*dt, int(f[3])))
            except (ValueError, IndexError):
                pass
    return np.array(rows).T

RUNS = [('run_linic_reorth/run.log', 'cref build (old install lib)', '#d62728', '-'),
        ('run_linic_h1verify/run.log', 'h1 = 1.0 build (Nek parity)', '#1f77b4', '--')]

fig, ax = plt.subplots(1, 2, figsize=(13.5, 4.6))
for log, lab, c, ls in RUNS:
    t, eu, eT = errtrace(log)
    ax[0].semilogy(t, eu, color=c, ls=ls, lw=1.8, label=f'u — {lab}')
    ax[0].semilogy(t, eT, color=c, ls=ls, lw=1.0, alpha=0.5, label=f'T — {lab}')
ax[0].set_title('error decay: identical answers\n(solid vs dashed must coincide)')
ax[0].set_xlabel('t'); ax[0].set_ylabel('max error')
ax[0].grid(alpha=0.3, which='both'); ax[0].legend(fontsize=8)

for log, lab, c, ls in RUNS:
    t, it = iters(log)
    ax[1].plot(t, it, color=c, ls=ls, lw=1.0, alpha=0.8, label=lab)
ax[1].set_title('pressure iterations per step (tol 1e-7)')
ax[1].set_xlabel('t'); ax[1].set_ylabel('GMRES iterations')
ax[1].grid(alpha=0.3); ax[1].legend(fontsize=9)

fig.suptitle('h1 = 1.0 preconditioner verification — linear-IC transient, projection 4 + reorthogonalize',
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig('run_linic_h1verify/viz_h1_verify.png', dpi=110)
t1, e1, _ = errtrace(RUNS[0][0]); t2, e2, _ = errtrace(RUNS[1][0])
n = min(len(e1), len(e2))
print('wrote run_linic_h1verify/viz_h1_verify.png')
print('max |u-err difference| between builds over the run: %.2e' % np.max(np.abs(e1[:n]-e2[:n])))
