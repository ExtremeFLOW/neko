# GIF: streamwise (x,y) plane at z = pi/2 for the muT soak, one frame per output.
# Panels: u' = u - <u>(y)  (fluctuation, plane-averaged mean removed per frame),
#         v, T, mu. Axes labeled; mean flow is +x (left -> right).
# Run: /home/hochi/projects/matrixenv/bin/python render_muT_gif_xy.py
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
zlev = np.unique(np.round(z, 8))
mz = np.isclose(z, zlev[np.argmin(np.abs(zlev - np.pi/2))])
xi, yi = x[mz], y[mz]

# y-level machinery for the per-frame plane average <u>(y) over the FULL field
yr = np.round(y, 8)
ylev, yinv = np.unique(yr, return_inverse=True)
ycnt = np.bincount(yinv)

n_frames = len(sorted(glob.glob(RUN + 'field0.f0*')))
specs = [(r"u' = u - $\langle u\rangle(y)$", 'RdBu_r', np.linspace(-0.25, 0.25, 41)),
         ('v', 'RdBu_r', np.linspace(-0.3, 0.3, 41)),
         ('T', 'RdBu_r', np.linspace(1.0, 4.0, 41)),
         (r'$\mu$', 'viridis', np.linspace(3.5e-4, 9.5e-4, 41))]

frames = []
for i in range(n_frames):
    fd = preadnek(RUN + 'field0.f%05d' % i, comm)
    u_full = cat(fd, 'vel', 0)
    u_mean = (np.bincount(yinv, u_full) / ycnt)[yinv]      # <u>(y) mapped back
    up = (u_full - u_mean)[mz]
    v, T = cat(fd, 'vel', 1)[mz], cat(fd, 'temp', 0)[mz]
    pr = preadnek(RUN + 'props0.f%05d' % i, comm)
    mu = cat(pr, 'vel', 0)[mz]

    fig, axs = plt.subplots(4, 1, figsize=(10, 10), sharex=True)
    for a, (lab, cm, lev), val in zip(axs, specs, [up, v, T, mu]):
        a.tricontourf(xi, yi, np.clip(val, lev[0], lev[-1]), levels=lev,
                      cmap=cm, extend='both')
        a.set_title(lab, fontsize=11)
        a.set_ylabel('y')
        a.set_aspect('equal')
    axs[-1].set_xlabel('x  (mean flow →)')
    fig.suptitle('Rayleigh–Bénard–Poiseuille channel, μ(T) soak — streamwise plane z = π/2,  '
                 f't = {i*DT_OUT:.0f}   (hot wall y = 0, cold wall y = 2)',
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=85)
    plt.close(fig)
    buf.seek(0)
    frames.append(Image.open(buf).convert('P', palette=Image.ADAPTIVE))
    if i % 10 == 0:
        print(f'frame {i}/{n_frames-1}')

frames[0].save(RUN + 'anim_muT_streamwise.gif', save_all=True,
               append_images=frames[1:], duration=220, loop=0)
print('wrote', RUN + 'anim_muT_streamwise.gif', f'({n_frames} frames)')
frames[len(frames)//3].convert('RGB').save(RUN + 'anim_muT_streamwise_sample.png')
print('wrote sample still')
