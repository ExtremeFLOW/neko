# GIF animation of the muT soak: u, v, T, mu in the cross-stream (z,y) plane at
# x = pi, one frame per output (t = 0..200, every 5 tu), fixed color scales.
# Run: /home/hochi/projects/matrixenv/bin/python render_muT_gif.py
import glob, io
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from pysemtools.io.ppymech.neksuite import preadnek
from mpi4py import MPI
comm = MPI.COMM_WORLD
RUN = 'overnight_gravity_muT/'
DT_OUT = 5.0

def cat(fd, a, k):
    return np.concatenate([getattr(fd.elem[e], a)[k].ravel() for e in range(fd.nel)])

geo = preadnek(RUN + 'field0.f00000', comm)
x, y, z = (cat(geo, 'pos', k) for k in range(3))
xlev = np.unique(np.round(x, 8))
mx = np.isclose(x, xlev[np.argmin(np.abs(xlev - np.pi))])
zi, yi = z[mx], y[mx]

n_frames = len(sorted(glob.glob(RUN + 'field0.f0*')))
specs = [('u', 'viridis', np.linspace(0.0, 1.5, 41)),
         ('v', 'RdBu_r', np.linspace(-0.3, 0.3, 41)),
         ('T', 'RdBu_r', np.linspace(1.0, 4.0, 41)),
         (r'$\mu$', 'viridis', np.linspace(3.5e-4, 9.5e-4, 41))]

frames = []
for i in range(n_frames):
    fd = preadnek(RUN + 'field0.f%05d' % i, comm)
    u, v, T = cat(fd, 'vel', 0)[mx], cat(fd, 'vel', 1)[mx], cat(fd, 'temp', 0)[mx]
    pr = preadnek(RUN + 'props0.f%05d' % i, comm)
    mu = cat(pr, 'vel', 0)[mx]
    fig, axs = plt.subplots(2, 2, figsize=(11, 7))
    for a, (lab, cm, lev), val in zip(axs.ravel(), specs, [u, v, T, mu]):
        a.tricontourf(zi, yi, np.clip(val, lev[0], lev[-1]), levels=lev,
                      cmap=cm, extend='both')
        a.set_title(lab, fontsize=11)
        a.set_aspect('equal')
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f'Rayleigh–Bénard–Poiseuille channel, μ(T) soak — cross-plane x = π,  t = {i*DT_OUT:.0f}',
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=85)
    plt.close(fig)
    buf.seek(0)
    frames.append(Image.open(buf).convert('P', palette=Image.ADAPTIVE))
    if i % 10 == 0:
        print(f'frame {i}/{n_frames-1}')

frames[0].save(RUN + 'anim_muT.gif', save_all=True, append_images=frames[1:],
               duration=220, loop=0)
print('wrote', RUN + 'anim_muT.gif', f'({n_frames} frames)')
# also dump one representative still for a quick visual check
frames[len(frames)//3].convert('RGB').save(RUN + 'anim_muT_sample.png')
print('wrote sample still anim_muT_sample.png')
