#!/usr/bin/env python3
"""
Slide figure: two-phase turbulent channel flow — case setup overview.

Shows domain geometry, turbulent inlet profile, drop initialisation,
boundary conditions, and key dimensions.

Usage:
    python3 figure_case_setup.py
Output:
    case_setup.png  (200 dpi, white background, slide-ready)
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ── physical parameters ───────────────────────────────────────────────────────
Lx_phys = 4 * np.pi        # streamwise  ≈ 12.57
Lz_phys = 4 * np.pi / 3    # spanwise    ≈ 4.19
R       = 0.3               # drop radius
Re_tau  = 180.0
Re_b    = 2800.0

# ── display coordinate mapping ────────────────────────────────────────────────
# mpl (x, y, z):  x=streamwise, y=spanwise, z=wall-normal [0,2]
cx   = 0.38
Lx   = Lx_phys * cx      # ≈ 4.77  compressed streamwise
Ly   = Lz_phys            # ≈ 4.19  spanwise
Lz   = 2.0                # wall-normal

xc, yc, zc = Lx/2, Ly/2, 1.0      # drop centre

# ── colours ───────────────────────────────────────────────────────────────────
WALL_FC   = '#b0cfe8'
WALL_EC   = '#5a8fb8'
DROP_C    = '#cc1111'
PROF_C    = '#9e1a1a'
ARROW_C   = '#0d3a5c'
DIM_C     = '#1a1a2e'
EDGE_C    = '#4a7fa8'

# ── figure / axes ─────────────────────────────────────────────────────────────
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 12,
                     'figure.facecolor': 'white'})
fig = plt.figure(figsize=(14, 7), facecolor='white')
ax  = fig.add_subplot(111, projection='3d')
ax.set_facecolor('white')
for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
    pane.fill = False
    pane.set_edgecolor((0, 0, 0, 0))
ax.grid(False)
ax.set_axis_off()
ax.view_init(elev=18, azim=210)

# box corners
C = np.array([
    [0,  0,  0 ],  # 0 bottom-front-left
    [Lx, 0,  0 ],  # 1 bottom-front-right
    [Lx, Ly, 0 ],  # 2 bottom-back-right
    [0,  Ly, 0 ],  # 3 bottom-back-left
    [0,  0,  Lz],  # 4 top-front-left
    [Lx, 0,  Lz],  # 5 top-front-right
    [Lx, Ly, Lz],  # 6 top-back-right
    [0,  Ly, Lz],  # 7 top-back-left
])

def face(verts, fc, ec, alpha, lw=0.9):
    ax.add_collection3d(Poly3DCollection(
        [np.array(verts, float)], alpha=alpha,
        facecolors=[fc], edgecolors=[ec], linewidths=lw))

def edge(p1, p2, **kw):
    ax.plot([p1[0],p2[0]], [p1[1],p2[1]], [p1[2],p2[2]], **kw)

# ── walls: light fill, bottom only to anchor the scene ───────────────────────
face(C[[0,1,2,3]], WALL_FC, WALL_EC, alpha=0.22, lw=1.2)   # bottom wall
face(C[[4,5,6,7]], WALL_FC, WALL_EC, alpha=0.07, lw=1.2)   # top wall (almost transparent)

# ── box edges ─────────────────────────────────────────────────────────────────
ekw = dict(color=EDGE_C, lw=1.3)
for i, j in [(0,1),(1,5),(5,4),(4,0),   # front face
             (3,2),(2,6),(6,7),(7,3),   # back face
             (0,3),(1,2),(4,7),(5,6)]:  # horizontal connections
    edge(C[i], C[j], **ekw)

# ── turbulent velocity profile ─────────────────────────────────────────────────
# Drawn on the y=0 front face, protruding OUTWARD in -y direction.
# This makes the profile stand clearly in front of the domain.
k0, C0 = 0.41, 5.17
y_w    = np.linspace(-1, 1, 200)
yp     = (1.0 - np.abs(y_w)) * Re_tau
u_r    = (1/k0 * np.log(1 + k0*yp) +
          (C0 - np.log(k0)/k0) *
          (1 - np.exp(-yp/11) - yp/11 * np.exp(-yp/3))
         ) * Re_tau / Re_b

# Scale factor: converts velocity to display units in the -y direction.
# Tuned so the peak (~1.17) protrudes about 1.2 display units.
scale = 1.03
u_vis  = u_r * scale      # profile offset in -y direction (positive = outward)
z_prof = y_w + 1.0        # wall-normal (0 to 2)

# Profile attached at x=0 (inlet face), front edge (y=0), protrudes into y<0
rects = []
for i in range(len(y_w) - 1):
    rects.append([
        [0,  0,          z_prof[i]  ],
        [0, -u_vis[i],   z_prof[i]  ],
        [0, -u_vis[i+1], z_prof[i+1]],
        [0,  0,          z_prof[i+1]],
    ])
ax.add_collection3d(Poly3DCollection(rects, alpha=0.70,
                                      facecolors=PROF_C, edgecolors='none'))
# Outline curve
ax.plot(np.zeros_like(u_vis), -u_vis, z_prof, color=PROF_C, lw=2.5, zorder=10)
# Vertical baseline at x=0, y=0
ax.plot([0, 0], [0, 0], [0, Lz], color=EDGE_C, lw=1.5)

# ── drop sphere ───────────────────────────────────────────────────────────────
# Draw back hemisphere first, then front hemisphere — forces correct depth order.
u_sp = np.linspace(0, 2*np.pi, 72)
v_sp = np.linspace(0, np.pi,   40)
Xs = xc + R * np.outer(np.cos(u_sp), np.sin(v_sp))
Ys = yc + R * np.outer(np.sin(u_sp), np.sin(v_sp))
Zs = zc + R * np.outer(np.ones(72),  np.cos(v_sp))

# Back half (facing away from viewer — render first so front overwrites)
mask_back = np.cos(u_sp - np.radians(210)) < 0   # rough back-facing rows
Xs_b = Xs[mask_back, :]
Ys_b = Ys[mask_back, :]
Zs_b = Zs[mask_back, :]
ax.plot_surface(Xs_b, Ys_b, Zs_b, color=DROP_C, alpha=0.9,
                rstride=1, cstride=1, linewidth=0, shade=True,
                antialiased=True)

# Front half (facing viewer — render last so it stays on top)
mask_front = ~mask_back
Xs_f = np.vstack([Xs[mask_front, :], Xs[mask_front[-1:], :]])
Ys_f = np.vstack([Ys[mask_front, :], Ys[mask_front[-1:], :]])
Zs_f = np.vstack([Zs[mask_front, :], Zs[mask_front[-1:], :]])
ax.plot_surface(Xs_f, Ys_f, Zs_f, color=DROP_C, alpha=1.0,
                rstride=1, cstride=1, linewidth=0, shade=True,
                antialiased=True)

# Equatorial ring for 3-D cue
t = np.linspace(0, 2*np.pi, 120)
ax.plot(xc + R*np.cos(t), yc + R*np.sin(t), np.full(120, zc),
        color='#ffffff', lw=1.8, alpha=0.80, zorder=9)
# Meridian in the streamwise-wallnormal plane (most visible from this view)
ax.plot(xc + R*np.cos(t), np.full(120, yc), zc + R*np.sin(t),
        color='#ffffff', lw=1.4, alpha=0.65, zorder=9)

# ── drop leader dashed line ───────────────────────────────────────────────────
ax.plot([xc + R*0.4, xc + R + 0.35],
        [yc - R*0.3, yc - R - 0.40],
        [zc + R*0.7, zc + R + 0.55],
        color=DROP_C, lw=1.2, alpha=0.80, linestyle='--', zorder=10)

# ── flow arrow ─────────────────────────────────────────────────────────────────
x0 = Lx * 0.06
dx = Lx * 0.28
yarr = 0.55
zarr = zc + 0.05
ax.plot([x0, x0 + dx], [yarr, yarr], [zarr, zarr],
        color=ARROW_C, lw=3.5, solid_capstyle='round', zorder=11)
ax.quiver(x0 + dx * 0.72, yarr, zarr, dx * 0.28, 0, 0,
          color=ARROW_C, arrow_length_ratio=0.55, linewidth=3.5, zorder=11)

# ── dimension lines ───────────────────────────────────────────────────────────
dz = -0.42
# Lx (streamwise)
ax.plot([0, Lx], [0, 0], [dz]*2,             color=DIM_C, lw=1.1)
ax.plot([0,  0], [0, 0], [dz-0.07, dz+0.07], color=DIM_C, lw=1.1)
ax.plot([Lx,Lx],[0, 0], [dz-0.07, dz+0.07],  color=DIM_C, lw=1.1)
# Ly (spanwise — runs in y direction)
ax.plot([0, 0], [0, Ly], [dz]*2,             color=DIM_C, lw=1.1)
ax.plot([0, 0], [0,  0], [dz-0.07, dz+0.07], color=DIM_C, lw=1.1)
ax.plot([0, 0],[Ly, Ly], [dz-0.07, dz+0.07], color=DIM_C, lw=1.1)
# 2h (wall-normal height)
dx_off = -0.45
ax.plot([dx_off]*2, [0, 0], [0, Lz],                       color=DIM_C, lw=1.1)
ax.plot([dx_off-0.09, dx_off+0.09], [0,0], [0,0],          color=DIM_C, lw=1.1)
ax.plot([dx_off-0.09, dx_off+0.09], [0,0], [Lz,Lz],        color=DIM_C, lw=1.1)

# ── display limits ────────────────────────────────────────────────────────────
prof_max = float(np.max(u_vis))
ax.set_xlim(-0.8, Lx + 0.8)
ax.set_ylim(-prof_max - 0.4, Ly + 0.8)
ax.set_zlim(-0.9, Lz + 0.8)
ax.set_box_aspect([Lx, Ly + prof_max, Lz])

# ── ALL LABELS in 2D figure coordinates ──────────────────────────────────────
# Tuned for view_init(elev=18, azim=210). Adjust if view_init changes.

def lab(x, y, s, **kw):
    defaults = dict(transform=fig.transFigure, ha='center', va='center',
                    fontsize=11, color='#1a1a2e')
    defaults.update(kw)
    fig.text(x, y, s, **defaults)

# No-slip walls
lab(0.44, 0.84, 'no-slip wall', fontsize=12, fontweight='bold', color='#154360')
lab(0.40, 0.41, 'no-slip wall', fontsize=12, fontweight='bold', color='#154360')

# Dimension labels
lab(0.46, 0.30, r'$L_x = 4\pi \approx 12.6\,h$',
    fontsize=11, fontweight='bold', color=DIM_C)
lab(0.17, 0.26, r'$L_z = \dfrac{4\pi}{3} \approx 4.2\,h$',
    fontsize=11, fontweight='bold', color=DIM_C)
lab(0.10, 0.60, r'$2h$',
    fontsize=13, fontweight='bold', color=DIM_C)

# Drop annotation (upper right, aligned with dashed leader)
lab(0.64, 0.84, r'drop:  $R = 0.3\,h$',
    fontsize=12, fontweight='bold', color=DROP_C, ha='left')

# Flow velocity label (above the arrow)
lab(0.35, 0.68, r'$U_b = 1$',
    fontsize=13, fontweight='bold', color=ARROW_C)

# Turbulent IC label — positioned near the protruding profile
# With azim=210 the profile protrudes toward the viewer's right
lab(0.23, 0.27, 'turbulent IC\n(Reichardt profile\n+ perturbations)',
    fontsize=10, color=PROF_C, ha='center', va='top', linespacing=1.6)

# Periodic BCs
lab(0.82, 0.60, 'periodic $(x)$',
    fontsize=10, color='#555', style='italic', ha='left')
lab(0.65, 0.79, 'periodic $(z)$',
    fontsize=10, color='#555', style='italic', ha='left')

# Re / mesh info box
fig.text(0.97, 0.05,
         r'$Re_b = 2800$,  $Re_\tau = 180$'
         '\n'
         r'mesh: $81 \times 18 \times 27$,  $N = 7$',
         transform=fig.transFigure,
         color='#2c2c2c', ha='right', va='bottom', fontsize=11,
         bbox=dict(boxstyle='round,pad=0.45', facecolor='#f0f5fa',
                   edgecolor='#b0c4d8', linewidth=0.9))

plt.savefig('case_setup.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print('Saved case_setup.png')
