# Two long GIFs over the FULL 53K story: Old_53K_turbulent (t=0-50) +
# GJP_temperature_stab continuation (t=51-76) = 77 frames, mid-z plane,
# solid-hump masking, fixed scales, "scalar GJP ON" badge from t=51.
#   anim_flow_full.gif : u (top) + v (bottom)
#   anim_temp_full.gif : T
# Run: /home/hochi/projects/matrixenv/bin/python render_story_gifs.py
import glob, io
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from PIL import Image
from pysemtools.io.ppymech.neksuite import preadnek
from mpi4py import MPI
comm = MPI.COMM_WORLD
A = '../heated_hump_res/Old_53K_turbulent/'
B = '../heated_hump_res/GJP_temperature_stab/'
OUT = '../heated_hump_res/'

def cat(fd, a, k):
    return np.concatenate([getattr(fd.elem[e], a)[k].ravel() for e in range(fd.nel)])

geo = preadnek(A + 'field0.f00000', comm)
x, y, z = (cat(geo, 'pos', k) for k in range(3))
zl = np.unique(np.round(z, 8))
mz = np.isclose(z, zl[np.argmin(np.abs(zl - zl.mean()))])
xi, yi = x[mz], y[mz]
def yw(xq):
    w = np.zeros_like(xq); m = (xq > 4.0) & (xq < 8.0)
    w[m] = 0.5*(1.0 + np.cos(np.pi*(xq[m]-6.0)/2.0)); return w
tri = mtri.Triangulation(xi, yi)
xc = xi[tri.triangles].mean(1); yc = yi[tri.triangles].mean(1)
tri.set_mask(yc < yw(xc) - 1e-3)
xs = np.linspace(0, 25, 600)

files = sorted(glob.glob(A + 'field0.f0*')) + sorted(glob.glob(B + 'field0.f0*'))
files = [f for i, f in enumerate(files) if int(f[-5:]) == i]   # contiguous 0..76
print(len(files), 'frames, last =', files[-1])

def solid(a):
    a.fill_between(xs, 0, yw(xs), color='0.35', zorder=3)
    a.plot(xs, yw(xs), color='k', lw=0.5, zorder=4)

flow_frames, temp_frames = [], []
for i, f in enumerate(files):
    fd = preadnek(f, comm)
    u, v, T = cat(fd, 'vel', 0)[mz], cat(fd, 'vel', 1)[mz], cat(fd, 'temp', 0)[mz]
    badge = '   [scalar GJP ON]' if i >= 51 else ''

    fig, ax = plt.subplots(2, 1, figsize=(12.5, 6))
    ax[0].tricontourf(tri, u, levels=np.linspace(-1.0, 2.0, 43), cmap='viridis', extend='both')
    ax[0].set_title('u')
    ax[1].tricontourf(tri, v, levels=np.linspace(-0.8, 0.8, 43), cmap='RdBu_r', extend='both')
    ax[1].set_title('v')
    for a in ax:
        solid(a); a.set_aspect('equal'); a.set_ylabel('y')
    ax[1].set_xlabel('x  (mean flow →)')
    fig.suptitle(f'Heated hump 53K, Re 3000 — flow field,  t = {i}{badge}', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    buf = io.BytesIO(); fig.savefig(buf, format='png', dpi=75); plt.close(fig)
    buf.seek(0); flow_frames.append(Image.open(buf).convert('P', palette=Image.ADAPTIVE))

    fig, a = plt.subplots(figsize=(12.5, 3.4))
    a.tricontourf(tri, T, levels=np.linspace(1.0, 2.0, 43), cmap='RdBu_r', extend='both')
    solid(a); a.set_aspect('equal'); a.set_xlabel('x'); a.set_ylabel('y')
    a.set_title(f'temperature (hot wall x > 11),  t = {i}{badge}', fontsize=11)
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format='png', dpi=75); plt.close(fig)
    buf.seek(0); temp_frames.append(Image.open(buf).convert('P', palette=Image.ADAPTIVE))
    if i % 10 == 0:
        print(f'frame {i}/{len(files)-1}')

flow_frames[0].save(OUT + 'anim_flow_full.gif', save_all=True,
                    append_images=flow_frames[1:], duration=160, loop=0)
temp_frames[0].save(OUT + 'anim_temp_full.gif', save_all=True,
                    append_images=temp_frames[1:], duration=160, loop=0)
print('wrote', OUT + 'anim_flow_full.gif and anim_temp_full.gif')
flow_frames[60].convert('RGB').save(OUT + 'anim_flow_sample.png')
temp_frames[60].convert('RGB').save(OUT + 'anim_temp_sample.png')
print('samples written (frame 60, GJP era)')
