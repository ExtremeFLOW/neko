# GIF of the Re 3000 turbulent heated hump (Old_53K_turbulent frames, order 6):
# u, T, v in the mid-z plane, one frame per output. The triangulation is masked
# below the lower-wall line (recovered from the mesh point cloud), so the hump
# appears as solid geometry instead of interpolated fill.
# Run: /home/hochi/projects/matrixenv/bin/python render_hump_gif.py
import glob, io
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from PIL import Image
from pysemtools.io.ppymech.neksuite import preadnek
from mpi4py import MPI
comm = MPI.COMM_WORLD
RUN = '../heated_hump_res/Old_53K_turbulent/'

def cat(fd, a, k):
    return np.concatenate([getattr(fd.elem[e], a)[k].ravel() for e in range(fd.nel)])

geo = preadnek(RUN + 'field0.f00000', comm)
x, y, z = (cat(geo, 'pos', k) for k in range(3))
zlev = np.unique(np.round(z, 8))
mz = np.isclose(z, zlev[np.argmin(np.abs(zlev - zlev.mean()))])
xi, yi = x[mz], y[mz]

# exact wall line from the mesh generator (hump.geo): cosine bump
#   y = 0.5*Hh*(1 + cos(pi*(x - xc)/wh)) on [x1, x2] = [4, 8], else 0
def ywall_fn(xq):
    yw = np.zeros_like(xq)
    m = (xq > 4.0) & (xq < 8.0)
    yw[m] = 0.5 * 1.0 * (1.0 + np.cos(np.pi * (xq[m] - 6.0) / 2.0))
    return yw

tri = mtri.Triangulation(xi, yi)
xc = xi[tri.triangles].mean(axis=1)
yc = yi[tri.triangles].mean(axis=1)
mask = yc < ywall_fn(xc) - 1e-3
tri.set_mask(mask)
xu = np.linspace(xi.min(), xi.max(), 800)
ywall = ywall_fn(xu)
print(f'slice pts: {xi.size}, triangles masked: {mask.sum()}/{mask.size}')

specs = [('u', 'viridis', np.linspace(-1.0, 2.0, 43)),
         ('T', 'RdBu_r', np.linspace(1.0, 2.0, 43)),
         ('v', 'RdBu_r', np.linspace(-0.8, 0.8, 43))]

files = sorted(glob.glob(RUN + 'field0.f0*'))
frames = []
for i, f in enumerate(files):
    fd = preadnek(f, comm)
    u, v, T = cat(fd, 'vel', 0)[mz], cat(fd, 'vel', 1)[mz], cat(fd, 'temp', 0)[mz]
    fig, axs = plt.subplots(3, 1, figsize=(13, 7.6), sharex=True)
    for a, (lab, cm, lev), val in zip(axs, specs, [u, T, v]):
        a.tricontourf(tri, np.clip(val, lev[0], lev[-1]), levels=lev,
                      cmap=cm, extend='both')
        a.fill_between(xu, 0, ywall, color='0.35', zorder=3)  # solid hump
        a.plot(xu, ywall, color='k', lw=0.6, zorder=4)        # wall outline
        a.set_title(lab, fontsize=11)
        a.set_ylabel('y'); a.set_aspect('equal')
    axs[-1].set_xlabel('x  (mean flow →)')
    fig.suptitle(f'Heated hump, Re 3000, order 6, Dong outflow — mid-z plane,  t = {i:.0f}',
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=75)
    plt.close(fig)
    buf.seek(0)
    frames.append(Image.open(buf).convert('P', palette=Image.ADAPTIVE))
    if i % 10 == 0:
        print(f'frame {i}/{len(files)-1}')

frames[0].save(RUN + 'anim_hump_turb.gif', save_all=True, append_images=frames[1:],
               duration=250, loop=0)
print('wrote', RUN + 'anim_hump_turb.gif', f'({len(frames)} frames)')
frames[-1].convert('RGB').save(RUN + 'anim_hump_turb_sample.png')
print('wrote sample still')
