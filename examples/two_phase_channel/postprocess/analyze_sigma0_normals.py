"""
analyze_sigma0_normals.py — CDI quality and κ accuracy analysis for σ=0 run.

Investigates three hypotheses for the κ_rms spike (6→64 in 0.4 TU):
  A. Genuine σ=0 drop deformation under turbulent strain (no restoring force)
  B. SEM element-boundary artifact in κ = −∇·n̂
  C. CDI startup transient (interface not fully established at t=20.002)

CDI quality metrics computed per snapshot:
  1. |∇φ| at interface (φ∈[0.3,0.7]): theoretical steady-state = 1/(2ε) = 7.14
  2. Interface thickness: distance between φ=0.1 and φ=0.9 along radial rays = 4ε = 0.28
  3. Volume conservation: ∫φ dV (should be constant across snapshots)
  4. Sphericity: std of r = |x − x_c| on φ≈0.5 iso-surface
  5. n̂ alignment: angle between computed n̂ and radial direction n̂_radial = (x−x_c)/r
  6. Numerical κ at interface: compare to analytical 2/R = 6.67 for sphere

Usage:
    python3 analyze_sigma0_normals.py
    python3 analyze_sigma0_normals.py --run channel_test_sigma0_diag --fps 4
    python3 analyze_sigma0_normals.py --run channel_test_sigma0 --no-animation
"""

import argparse
import glob
import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek
from pysemtools.datatypes.msh import Mesh as msh_c
from pysemtools.datatypes.field import Field as field_c

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--run', default='channel_test_sigma0',
                    help='Run directory name under /lscratch/sieburgh/simulations/')
parser.add_argument('--fps', type=int, default=3,
                    help='Frames per second in output GIF')
parser.add_argument('--dpi', type=int, default=120)
parser.add_argument('--no-animation', action='store_true',
                    help='Skip animation output (just produce metric plots)')
args = parser.parse_args()

RUN_DIR  = f'/lscratch/sieburgh/simulations/{args.run}'
OUT_DIR  = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
comm = MPI.COMM_WORLD

# ---------------------------------------------------------------------------
# Physical / CDI parameters (from the case file)
# ---------------------------------------------------------------------------
LLX = 4.0 * math.pi          # streamwise length ≈ 12.566
LLY = 2.0                    # wall-normal half-channel  [-1, 1]
LLZ = 4.0 / 3.0 * math.pi   # spanwise length ≈ 4.189
X_C = LLX / 2.0              # drop centre x ≈ 6.283
Y_C = 0.0                    # drop centre y = 0 (centreline)
Z_C = LLZ / 2.0              # drop centre z ≈ 2.094

R     = 0.3                  # drop radius
EPS   = 0.07                 # interface thickness parameter
GAMMA = 0.05                 # CDI compression coefficient
GRAD_PHI_THEORY = 1.0 / (4.0 * EPS)   # = 3.571  (CDI steady-state: |∇φ|=φ(1-φ)/ε, at φ=0.5: 0.25/0.07)
THICKNESS_THEORY = 4.0 * EPS          # = 0.280  (φ=0.1→0.9 thickness)
KAPPA_SPHERE = 2.0 / R                # = 6.667  (analytical curvature of sphere)

# ---------------------------------------------------------------------------
# Discover field files
# ---------------------------------------------------------------------------
field_files = sorted(glob.glob(os.path.join(RUN_DIR, 'field0.f[0-9]*')))
if not field_files:
    raise FileNotFoundError(f'No field0.f* files in {RUN_DIR}')
print(f'Found {len(field_files)} field files in {RUN_DIR}')

# ---------------------------------------------------------------------------
# Read mesh from first file
# ---------------------------------------------------------------------------
print('Reading mesh ...')
xyz_data = preadnek(field_files[0], comm)
msh = msh_c(comm, data=xyz_data)

n_elem = msh.x.shape[0]
lz_    = msh.x.shape[1]   # GLL pts along k (z) axis
ly_    = msh.x.shape[2]   # GLL pts along j (y) axis
lx_    = msh.x.shape[3]   # GLL pts along i (x) axis
print(f'  n_elem={n_elem}, lx=ly=lz={lx_}  (N={lx_-1})')

# Flatten all nodal positions once (shape: [n_pts, 3])
x_all = msh.x.reshape(-1)
y_all = msh.y.reshape(-1)
z_all = msh.z.reshape(-1)

# Element centre z, y for slice masks
iz_mid = lz_ // 2
iy_mid = ly_ // 2
ix_mid = lx_ // 2
z_elem_centers = msh.z[:, iz_mid, iy_mid, ix_mid]
y_elem_centers = msh.y[:, iz_mid, iy_mid, ix_mid]
x_elem_centers = msh.x[:, iz_mid, iy_mid, ix_mid]

# Approximate element widths
nz_elems = 27
nx_elems = 81
dz_approx = LLZ / nz_elems
dx_approx = LLX / nx_elems
dy_approx = LLY / 18.0 * 2.0   # 18 elements across full channel

# Slice masks (same logic as animate_two_phase_channel.py)
mask_xy = np.abs(z_elem_centers - Z_C) < 0.6 * dz_approx   # x-y at z=z_c
mask_xz = np.abs(y_elem_centers)        < 0.6 * dy_approx   # x-z at y=0

print(f'  x-y slice (z≈{Z_C:.2f}): {mask_xy.sum()} elements')
print(f'  x-z slice (y≈0):        {mask_xz.sum()} elements')


# ---------------------------------------------------------------------------
# Utility: build 2D triangulation for a slice
# ---------------------------------------------------------------------------
def build_triang(xs, ys):
    """Build matplotlib Triangulation from (n_el, nj, ni) GLL point arrays."""
    n_el, nj, ni = xs.shape
    pts_x, pts_y, tris = [], [], []
    offset = 0
    for e in range(n_el):
        pts_x.append(xs[e].flatten())
        pts_y.append(ys[e].flatten())
        for j in range(nj - 1):
            for i in range(ni - 1):
                p0 = offset + j * ni + i
                p1 = p0 + 1
                p2 = offset + (j + 1) * ni + i
                p3 = p2 + 1
                tris.append([p0, p1, p3])
                tris.append([p0, p3, p2])
        offset += nj * ni
    return tri.Triangulation(
        np.concatenate(pts_x),
        np.concatenate(pts_y),
        np.array(tris),
    )


triang_xy = build_triang(
    msh.x[mask_xy, iz_mid, :, :],
    msh.y[mask_xy, iz_mid, :, :],
)
triang_xz = build_triang(
    msh.x[mask_xz, :, iy_mid, :],
    msh.z[mask_xz, :, iy_mid, :],
)


# ---------------------------------------------------------------------------
# Utility: element-local gradient (∇φ) using np.gradient
#
# Assumes tensor-product element (valid for genmeshbox hex mesh): x varies
# only with i, y with j, z with k. Returns arrays of shape (n_elem, lz, ly, lx).
# ---------------------------------------------------------------------------
def compute_gradient(phi_arr, msh):
    """Compute element-local ∇φ on GLL grid. Returns (gphi_x, gphi_y, gphi_z)."""
    n_el, lz, ly, lx = phi_arr.shape
    gphi_x = np.zeros_like(phi_arr)
    gphi_y = np.zeros_like(phi_arr)
    gphi_z = np.zeros_like(phi_arr)

    for e in range(n_el):
        # 1-D coordinate arrays along each axis (tensor product assumption)
        x1d = msh.x[e, 0, 0, :]   # x varies with i (innermost index)
        y1d = msh.y[e, 0, :, 0]   # y varies with j
        z1d = msh.z[e, :, 0, 0]   # z varies with k (outermost index)

        phi_e = phi_arr[e]  # shape (lz, ly, lx) = (k, j, i)
        # np.gradient(f, z, y, x) → (df/dz, df/dy, df/dx) for f[k,j,i]
        gz, gy, gx = np.gradient(phi_e, z1d, y1d, x1d)
        gphi_x[e] = gx
        gphi_y[e] = gy
        gphi_z[e] = gz

    return gphi_x, gphi_y, gphi_z


# ---------------------------------------------------------------------------
# Utility: numerical κ = −∇·(∇φ/|∇φ|)
#
# Computes element-local div(n̂) using np.gradient on n̂ components.
# ---------------------------------------------------------------------------
def compute_kappa(phi_arr, msh, threshold=1e-3):
    """Compute κ = −∇·n̂ where n̂ = ∇φ/|∇φ|. Returns kappa array same shape."""
    gphi_x, gphi_y, gphi_z = compute_gradient(phi_arr, msh)
    mag = np.sqrt(gphi_x**2 + gphi_y**2 + gphi_z**2)
    mag = np.where(mag < threshold, threshold, mag)   # avoid /0

    nx_ = gphi_x / mag
    ny_ = gphi_y / mag
    nz_ = gphi_z / mag

    n_el, lz, ly, lx = phi_arr.shape
    div_n = np.zeros_like(phi_arr)

    for e in range(n_el):
        x1d = msh.x[e, 0, 0, :]
        y1d = msh.y[e, 0, :, 0]
        z1d = msh.z[e, :, 0, 0]

        _, _, dnx_dx = np.gradient(nx_[e], z1d, y1d, x1d)
        _, dny_dy, _ = np.gradient(ny_[e], z1d, y1d, x1d)
        dnz_dz, _, _ = np.gradient(nz_[e], z1d, y1d, x1d)

        div_n[e] = dnx_dx + dny_dy + dnz_dz

    return -div_n  # κ = −∇·n̂


# ---------------------------------------------------------------------------
# Utility: compute CDI quality metrics for one snapshot
# ---------------------------------------------------------------------------
def compute_metrics(phi_arr, msh, label=''):
    """
    Returns dict of CDI quality metrics.
    phi_arr: shape (n_elem, lz, ly, lx)
    """
    # Flatten
    phi_flat = phi_arr.reshape(-1)
    x_f = msh.x.reshape(-1)
    y_f = msh.y.reshape(-1)
    z_f = msh.z.reshape(-1)

    # Compute current drop centroid (mass-weighted centroid of φ field).
    # This is needed because the drop is advected by turbulence away from
    # the injection point (X_C, Y_C, Z_C).
    phi_sum = phi_flat.sum()
    if phi_sum > 0:
        xc_now = np.dot(phi_flat, x_f) / phi_sum
        yc_now = np.dot(phi_flat, y_f) / phi_sum
        zc_now = np.dot(phi_flat, z_f) / phi_sum
    else:
        xc_now, yc_now, zc_now = X_C, Y_C, Z_C
    centroid_displacement = math.sqrt(
        (xc_now - X_C)**2 + (yc_now - Y_C)**2 + (zc_now - Z_C)**2
    )

    # Radial distance from CURRENT drop centre (not injection point)
    r_f = np.sqrt((x_f - xc_now)**2 + (y_f - yc_now)**2 + (z_f - zc_now)**2)

    # Interface mask: φ ∈ [0.3, 0.7]
    interface_mask = (phi_flat >= 0.3) & (phi_flat <= 0.7)
    # Iso-surface mask: φ ≈ 0.5
    iso_mask = np.abs(phi_flat - 0.5) < 0.05

    # --- Metric 1: |∇φ| at interface ---
    gphi_x, gphi_y, gphi_z = compute_gradient(phi_arr, msh)
    mag_flat = np.sqrt(
        gphi_x.reshape(-1)**2 + gphi_y.reshape(-1)**2 + gphi_z.reshape(-1)**2
    )
    grad_mean = mag_flat[interface_mask].mean() if interface_mask.any() else np.nan
    grad_std  = mag_flat[interface_mask].std()  if interface_mask.any() else np.nan

    # --- Metric 2: Sphericity on φ=0.5 iso-surface ---
    r_iso = r_f[iso_mask]
    sph_mean = r_iso.mean() if iso_mask.any() else np.nan
    sph_std  = r_iso.std()  if iso_mask.any() else np.nan

    # --- Metric 3: Volume integral ∫φ dV (approximate: sum over all nodes) ---
    # Use element Jacobian: for the box mesh, vol_weight ≈ element_volume / (lx*ly*lz)
    # Approximate via trapezoidal: just use the nodal sum weighted by dx*dy*dz per element
    # (rough but sufficient for conservation check across snapshots)
    vol_phi = _volume_integral(phi_arr, msh)

    # --- Metric 4: n̂ alignment (angle between n̂_SEM and n̂_radial) ---
    mag_safe = np.where(mag_flat < 1e-3, 1e-3, mag_flat)
    nx_f = gphi_x.reshape(-1) / mag_safe
    ny_f = gphi_y.reshape(-1) / mag_safe
    nz_f = gphi_z.reshape(-1) / mag_safe

    # Use CURRENT centroid (not injection point) for radial normal direction
    r_safe = np.where(r_f < 1e-6, 1e-6, r_f)
    nr_x = (x_f - xc_now) / r_safe
    nr_y = (y_f - yc_now) / r_safe
    nr_z = (z_f - zc_now) / r_safe

    dot = nx_f * nr_x + ny_f * nr_y + nz_f * nr_z
    dot_clipped = np.clip(dot, -1, 1)
    angle_deg = np.degrees(np.arccos(np.abs(dot_clipped)))

    align_mean = angle_deg[interface_mask].mean() if interface_mask.any() else np.nan
    align_std  = angle_deg[interface_mask].std()  if interface_mask.any() else np.nan

    # --- Metric 5: Numerical κ on interface ---
    kappa_arr = compute_kappa(phi_arr, msh)
    kappa_flat = kappa_arr.reshape(-1)
    kappa_mean = kappa_flat[iso_mask].mean() if iso_mask.any() else np.nan
    kappa_rms  = np.sqrt((kappa_flat[iso_mask]**2).mean()) if iso_mask.any() else np.nan
    kappa_err  = kappa_mean - KAPPA_SPHERE  # bias vs spherical

    # --- Metric 6: tanh profile residual along radial rays ---
    # Sample φ along N_ray uniformly-spaced rays from current centroid.
    # Compare to φ_th(r) = 0.5*(1 + tanh((R-r)/(2*EPS))).
    # Reports mean L2 residual and interface width (r at φ=0.1 to φ=0.9).
    tanh_l2, width_mean = _tanh_profile_residual(
        phi_flat, x_f, y_f, z_f, xc_now, yc_now, zc_now
    )

    metrics = {
        'grad_mean': grad_mean, 'grad_std': grad_std,
        'sph_mean':  sph_mean,  'sph_std':  sph_std,
        'vol_phi':   vol_phi,
        'align_mean': align_mean, 'align_std': align_std,
        'kappa_mean': kappa_mean, 'kappa_rms': kappa_rms, 'kappa_err': kappa_err,
        'tanh_l2':    tanh_l2,    'width_mean': width_mean,
        'n_interface': interface_mask.sum(),
        'centroid': (xc_now, yc_now, zc_now),
        'centroid_disp': centroid_displacement,
        # Pass through for plotting
        'phi_flat': phi_flat, 'kappa_flat': kappa_flat,
        'mag_flat': mag_flat,
        'phi_arr': phi_arr, 'kappa_arr': kappa_arr, 'mag_arr': gphi_x,
        '_gphi': (gphi_x, gphi_y, gphi_z),
    }
    return metrics


def _tanh_profile_residual(phi_flat, x_f, y_f, z_f, xc, yc, zc,
                            n_ray=36, n_sample=60):
    """
    Sample φ along n_ray radial rays from drop centroid (xc, yc, zc).
    Fit against the CDI steady-state tanh profile:
        φ_th(r) = 0.5 * (1 + tanh((R - r) / (2*EPS)))
    Returns:
        tanh_l2   — mean L2 residual ||φ - φ_th||_2 averaged over rays
        width_mean — mean radial width between φ=0.1 and φ=0.9 (theory = 4*EPS)
    """
    r_max = 2.0 * R + 6.0 * EPS   # only need nodes within 2R+6ε of centroid
    r_all = np.sqrt((x_f - xc)**2 + (y_f - yc)**2 + (z_f - zc)**2)

    # Pre-filter to a small sphere — reduces search from ~20M to ~200k nodes
    near_mask = r_all <= r_max
    x_near  = x_f[near_mask]
    y_near  = y_f[near_mask]
    z_near  = z_f[near_mask]
    phi_near = phi_flat[near_mask]

    if x_near.size == 0:
        return np.nan, np.nan

    r_sample = np.linspace(0.0, r_max, n_sample)
    phi_theory = 0.5 * (1.0 + np.tanh((R - r_sample) / (2.0 * EPS)))

    # Ray directions: distribute uniformly on sphere using golden spiral
    golden = (1.0 + math.sqrt(5.0)) / 2.0
    l2_vals, width_vals = [], []
    for k in range(n_ray):
        theta = 2.0 * math.pi * k / golden
        cos_phi_s = 1.0 - 2.0 * (k + 0.5) / n_ray
        sin_phi_s = math.sqrt(max(0.0, 1.0 - cos_phi_s**2))
        dx = sin_phi_s * math.cos(theta)
        dy = sin_phi_s * math.sin(theta)
        dz = cos_phi_s

        # For each sample radius, find nearest node (search only in pre-filtered set)
        phi_ray = np.zeros(n_sample)
        for i, r in enumerate(r_sample):
            px = xc + r * dx;  py = yc + r * dy;  pz = zc + r * dz
            dist2 = (x_near - px)**2 + (y_near - py)**2 + (z_near - pz)**2
            phi_ray[i] = phi_near[dist2.argmin()]

        l2_vals.append(math.sqrt(((phi_ray - phi_theory)**2).mean()))

        r09 = _find_crossing(r_sample, phi_ray, 0.9)
        r01 = _find_crossing(r_sample, phi_ray, 0.1)
        if r09 is not None and r01 is not None and r01 > r09:
            width_vals.append(r01 - r09)

    tanh_l2   = float(np.mean(l2_vals)) if l2_vals else np.nan
    width_mean = float(np.mean(width_vals)) if width_vals else np.nan
    return tanh_l2, width_mean


def _find_crossing(r_arr, phi_arr, phi_target):
    """Linear interpolation: find r where phi_arr crosses phi_target (decreasing)."""
    for i in range(len(phi_arr) - 1):
        if (phi_arr[i] >= phi_target >= phi_arr[i + 1] or
                phi_arr[i] <= phi_target <= phi_arr[i + 1]):
            if abs(phi_arr[i + 1] - phi_arr[i]) < 1e-10:
                return r_arr[i]
            t = (phi_target - phi_arr[i]) / (phi_arr[i + 1] - phi_arr[i])
            return r_arr[i] + t * (r_arr[i + 1] - r_arr[i])
    return None


def _volume_integral(phi_arr, msh):
    """
    Approximate ∫φ dV using element-local trapezoidal quadrature.
    For a tensor-product hex mesh this is element_vol_sum ≈ volume_fraction × total_vol.
    """
    vol = 0.0
    n_el, lz, ly, lx = phi_arr.shape
    for e in range(n_el):
        x1d = msh.x[e, 0, 0, :]
        y1d = msh.y[e, 0, :, 0]
        z1d = msh.z[e, :, 0, 0]
        # Weight = product of np.gradient spacings (trapezoidal in each dir)
        dx_weights = np.gradient(x1d)
        dy_weights = np.gradient(y1d)
        dz_weights = np.gradient(z1d)
        # Outer product to get 3D volume weights
        ww = np.einsum('k,j,i->kji', dz_weights, dy_weights, dx_weights)
        vol += (phi_arr[e] * ww).sum()
    return vol


# ---------------------------------------------------------------------------
# Slice helper: extract field values on a 2D slice
# ---------------------------------------------------------------------------
def slice_values(arr, mask, gll_idx, axis):
    """
    Extract values from arr (n_elem, lz, ly, lx) on a slice.
    axis='z': fix k=gll_idx, result shape (n_elem_in_slice, ly, lx)
    axis='y': fix j=gll_idx, result shape (n_elem_in_slice, lz, lx)
    """
    sub = arr[mask]
    if axis == 'z':
        return sub[:, gll_idx, :, :]
    else:  # 'y'
        return sub[:, :, gll_idx, :]


# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------
records = []   # list of dicts with t + metrics

for fname in field_files:
    print(f'\nLoading {os.path.basename(fname)} ...')
    data = preadnek(fname, comm)
    fld  = field_c(comm, data=data)
    t    = fld.t

    phi_arr = fld.fields['scal'][0]   # (n_elem, lz, ly, lx)
    phi_flat = phi_arr.reshape(-1)
    print(f'  t={t:.4f}  phi_max={phi_flat.max():.4f}  phi_min={phi_flat.min():.4f}')

    # CDI quality metrics (slow: ~20s per snapshot due to element-loop gradient)
    print('  Computing gradient and metrics ...')
    m = compute_metrics(phi_arr, msh)
    m['t'] = t

    cx, cy, cz = m['centroid']
    print(f'  Drop centroid         = ({cx:.3f}, {cy:.3f}, {cz:.3f})  inj=({X_C:.3f},{Y_C:.3f},{Z_C:.3f})')
    print(f'  Centroid displacement = {m["centroid_disp"]:.4f}  (from injection point)')
    print(f'  |∇φ|_mean(interface)  = {m["grad_mean"]:.3f}  (theory={GRAD_PHI_THEORY:.3f})')
    print(f'  |∇φ|_std              = {m["grad_std"]:.3f}')
    print(f'  Sphericity r_mean     = {m["sph_mean"]:.4f}  (R={R})')
    print(f'  Sphericity r_std      = {m["sph_std"]:.4f}  (0=perfect sphere)')
    print(f'  ∫φ dV                 = {m["vol_phi"]:.5f}')
    print(f'  n̂ align angle_mean   = {m["align_mean"]:.2f}°  (0°=perfect radial)')
    print(f'  n̂ align angle_std    = {m["align_std"]:.2f}°')
    print(f'  κ_mean(iso-surface)  = {m["kappa_mean"]:.3f}  (sphere={KAPPA_SPHERE:.3f})')
    print(f'  κ_rms(iso-surface)   = {m["kappa_rms"]:.3f}')
    print(f'  κ error (bias)       = {m["kappa_err"]:+.3f}')
    print(f'  tanh L2 residual     = {m["tanh_l2"]:.4f}  (0=perfect CDI profile)')
    print(f'  interface width      = {m["width_mean"]:.4f}  (theory=4ε={THICKNESS_THEORY:.3f})')

    # ---------- Slice plots (one figure per snapshot) ----------
    phi_xy = slice_values(phi_arr,              mask_xy, iz_mid, 'z').reshape(-1)
    kap_xy = slice_values(m['kappa_arr'],       mask_xy, iz_mid, 'z').reshape(-1)
    mag_xy = slice_values(
        np.sqrt(m['_gphi'][0]**2 + m['_gphi'][1]**2 + m['_gphi'][2]**2),
        mask_xy, iz_mid, 'z'
    ).reshape(-1)

    phi_xz = slice_values(phi_arr,              mask_xz, iy_mid, 'y').reshape(-1)
    kap_xz = slice_values(m['kappa_arr'],       mask_xz, iy_mid, 'y').reshape(-1)
    mag_xz = slice_values(
        np.sqrt(m['_gphi'][0]**2 + m['_gphi'][1]**2 + m['_gphi'][2]**2),
        mask_xz, iy_mid, 'y'
    ).reshape(-1)

    fig, axes = plt.subplots(2, 3, figsize=(16, 7), dpi=args.dpi)
    fig.subplots_adjust(wspace=0.30, hspace=0.35)

    titles_xy = ['φ  (x–y, z=z_c)', '|∇φ|  (x–y)', 'κ_numerical  (x–y)']
    titles_xz = ['φ  (x–z, y=0)',   '|∇φ|  (x–z)', 'κ_numerical  (x–z)']

    # (values, colormap, vmin, vmax, colorbar_label)
    data_xy = [
        (phi_xy, 'RdBu_r',  0.0,   1.0,  'φ'),
        (mag_xy, 'plasma',  0.0,  15.0,  '|∇φ|'),
        (kap_xy, 'seismic', -35.0, 35.0, 'κ'),
    ]
    data_xz = [
        (phi_xz, 'RdBu_r',  0.0,   1.0,  'φ'),
        (mag_xz, 'plasma',  0.0,  15.0,  '|∇φ|'),
        (kap_xz, 'seismic', -35.0, 35.0, 'κ'),
    ]

    tcf_xy = []
    for col, (vals, cmap, vmin, vmax, _) in enumerate(data_xy):
        ax = axes[0, col]
        ax.set_aspect('equal')
        tcf = ax.tricontourf(triang_xy, vals, levels=60, cmap=cmap, vmin=vmin, vmax=vmax)
        tcf_xy.append(tcf)
        try:
            ax.tricontour(triang_xy, phi_xy, levels=[0.5], colors='k', linewidths=1.0)
        except Exception:
            pass
        ax.set(title=titles_xy[col], xlabel='x', ylabel='y')

    for col, (vals, cmap, vmin, vmax, _) in enumerate(data_xz):
        ax = axes[1, col]
        ax.set_aspect('equal')
        ax.tricontourf(triang_xz, vals, levels=60, cmap=cmap, vmin=vmin, vmax=vmax)
        try:
            ax.tricontour(triang_xz, phi_xz, levels=[0.5], colors='k', linewidths=1.0)
        except Exception:
            pass
        ax.set(title=titles_xz[col], xlabel='x', ylabel='z')

    # One colorbar per column, shared between xy and xz rows
    for col, (_, _, vmin, vmax, cb_label) in enumerate(data_xy):
        cbar = fig.colorbar(tcf_xy[col], ax=[axes[0, col], axes[1, col]],
                            shrink=0.85, pad=0.03, aspect=25)
        cbar.set_label(cb_label, fontsize=9)

    # Add element boundary lines on φ column (x-y slice)
    # Vertical lines at element x-boundaries
    n_elem_xy = mask_xy.sum()
    x_corners = msh.x[mask_xy, iz_mid, 0, 0]   # left face x for each element in slice
    y_corners_bot = msh.y[mask_xy, iz_mid, 0, :]  # bottom face y
    y_corners_top = msh.y[mask_xy, iz_mid, -1, :]
    x_unique = np.unique(np.round(msh.x[mask_xy, iz_mid, :, 0].reshape(-1), 4))
    y_unique = np.unique(np.round(msh.y[mask_xy, iz_mid, 0, :].reshape(-1), 4))
    for xv in x_unique[::1]:
        axes[0, 0].axvline(xv, color='gray', lw=0.3, alpha=0.5)
    for yv in y_unique:
        axes[0, 0].axhline(yv, color='gray', lw=0.3, alpha=0.5)

    run_name = os.path.basename(RUN_DIR)
    fn_idx = os.path.basename(fname).replace('field0.', '')
    fig.suptitle(
        f'{run_name}  t={t:.4f} TU   '
        f'κ_rms={m["kappa_rms"]:.2f}  r_std={m["sph_std"]:.4f}  '
        f'|∇φ|_mean={m["grad_mean"]:.2f} (th={GRAD_PHI_THEORY:.2f})',
        fontsize=10,
    )
    out_png = os.path.join(OUT_DIR, f'normals_{run_name}_{fn_idx}.png')
    plt.savefig(out_png, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_png}')

    # Store for summary plots
    records.append({
        't': t,
        'grad_mean': m['grad_mean'], 'grad_std': m['grad_std'],
        'sph_mean': m['sph_mean'],   'sph_std':  m['sph_std'],
        'vol_phi':  m['vol_phi'],
        'align_mean': m['align_mean'],
        'kappa_rms': m['kappa_rms'], 'kappa_mean': m['kappa_mean'],
        'centroid_disp': m['centroid_disp'],
        'tanh_l2': m['tanh_l2'],     'width_mean': m['width_mean'],
        'phi_xy': phi_xy, 'kap_xy': kap_xy, 'mag_xy': mag_xy,
        'phi_xz': phi_xz, 'kap_xz': kap_xz, 'mag_xz': mag_xz,
    })

# ---------------------------------------------------------------------------
# Summary time-series plot (if multiple snapshots)
# ---------------------------------------------------------------------------
if len(records) > 1:
    t_arr   = np.array([r['t']         for r in records])
    gm_arr  = np.array([r['grad_mean'] for r in records])
    sph_arr = np.array([r['sph_std']   for r in records])
    vol_arr = np.array([r['vol_phi']   for r in records])
    ali_arr = np.array([r['align_mean'] for r in records])
    kr_arr  = np.array([r['kappa_rms'] for r in records])
    l2_arr  = np.array([r['tanh_l2']   for r in records])
    wid_arr = np.array([r['width_mean'] for r in records])

    fig2, axs = plt.subplots(4, 2, figsize=(10, 11), dpi=args.dpi)
    fig2.subplots_adjust(hspace=0.40, wspace=0.30)

    axs[0, 0].plot(t_arr, gm_arr, 'o-')
    axs[0, 0].axhline(GRAD_PHI_THEORY, ls='--', color='red', label=f'theory={GRAD_PHI_THEORY:.2f}')
    axs[0, 0].set(ylabel='|∇φ|_mean (interface)', xlabel='t', title='CDI sharpness')
    axs[0, 0].legend(fontsize=8)

    axs[0, 1].plot(t_arr, sph_arr, 'o-')
    axs[0, 1].set(ylabel='r_std on φ=0.5 surface', xlabel='t',
                  title='Sphericity (0=perfect sphere)')

    axs[1, 0].plot(t_arr, vol_arr, 'o-')
    axs[1, 0].set(ylabel='∫φ dV', xlabel='t', title='Volume conservation')

    axs[1, 1].plot(t_arr, ali_arr, 'o-')
    axs[1, 1].set(ylabel='mean angle (°)', xlabel='t',
                  title='n̂ vs n̂_radial misalignment')

    axs[2, 0].plot(t_arr, kr_arr, 'o-')
    axs[2, 0].axhline(KAPPA_SPHERE, ls='--', color='red', label=f'sphere={KAPPA_SPHERE:.2f}')
    axs[2, 0].set(ylabel='κ_rms (iso-surface)', xlabel='t',
                  title='Curvature accuracy')
    axs[2, 0].legend(fontsize=8)

    # κ histogram (last snapshot)
    axs[2, 1].hist(records[-1]['kap_xy'], bins=80, range=(-40, 40))
    axs[2, 1].axvline(KAPPA_SPHERE, color='red', ls='--', label=f'κ_sphere={KAPPA_SPHERE:.2f}')
    axs[2, 1].set(xlabel='κ (x-y slice)', ylabel='count',
                  title=f'κ histogram  t={records[-1]["t"]:.3f}')
    axs[2, 1].legend(fontsize=8)

    # CDI profile quality
    axs[3, 0].plot(t_arr, l2_arr, 'o-')
    axs[3, 0].set(ylabel='tanh L2 residual', xlabel='t',
                  title='CDI profile quality (0=perfect tanh)')

    axs[3, 1].plot(t_arr, wid_arr, 'o-')
    axs[3, 1].axhline(THICKNESS_THEORY, ls='--', color='red',
                      label=f'4ε={THICKNESS_THEORY:.3f}')
    axs[3, 1].set(ylabel='interface width (φ:0.1→0.9)', xlabel='t',
                  title='CDI interface thickness')
    axs[3, 1].legend(fontsize=8)

    run_name = os.path.basename(RUN_DIR)
    fig2.suptitle(f'{run_name}: CDI quality metrics time series')
    out_ts = os.path.join(OUT_DIR, f'normals_timeseries_{run_name}.png')
    plt.savefig(out_ts, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig2)
    print(f'\nSaved time-series: {out_ts}')

# ---------------------------------------------------------------------------
# Animation: φ and κ on x-y slice over time
# ---------------------------------------------------------------------------
if not args.no_animation and len(records) > 1:
    run_name = os.path.basename(RUN_DIR)
    out_gif  = os.path.join(OUT_DIR, f'animate_sigma0_kappa_{run_name}.gif')

    kap_global_max = max(np.abs(r['kap_xy']).max() for r in records)
    kap_scale = min(kap_global_max, 30.0)   # clip for visual clarity

    fig_a, axes_a = plt.subplots(1, 2, figsize=(12, 4), dpi=args.dpi)
    fig_a.subplots_adjust(right=0.90, wspace=0.10)

    ax_phi_a, ax_kap_a = axes_a
    for ax in axes_a:
        ax.set_aspect('equal')
    ax_phi_a.set(title='φ  (x–y, z=z_c)', xlabel='x', ylabel='y')
    ax_kap_a.set(title='κ  (x–y, z=z_c)', xlabel='x', ylabel='y')

    title_a = fig_a.suptitle('', fontsize=10)

    def _update_anim(frame):
        r = records[frame]
        for ax in axes_a:
            for coll in ax.collections:
                coll.remove()
        axes_a[0].tricontourf(triang_xy, r['phi_xy'], levels=60,
                               cmap='RdBu_r', vmin=0, vmax=1)
        axes_a[1].tricontourf(triang_xy, r['kap_xy'], levels=60,
                               cmap='seismic', vmin=-kap_scale, vmax=kap_scale)
        for ax, phi_vals in zip(axes_a, [r['phi_xy'], r['phi_xy']]):
            try:
                ax.tricontour(triang_xy, phi_vals, levels=[0.5],
                              colors='k', linewidths=1.0)
            except Exception:
                pass
        title_a.set_text(
            f't = {r["t"]:.4f}   κ_rms(x-y)={np.sqrt((r["kap_xy"]**2).mean()):.2f}   '
            f'r_std={r["sph_std"]:.4f}'
        )

    anim = FuncAnimation(fig_a, _update_anim, frames=len(records),
                         interval=int(1000 / args.fps), repeat=True)
    anim.save(out_gif, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    plt.close(fig_a)
    print(f'Saved animation: {out_gif}')

# ---------------------------------------------------------------------------
# Decision gate summary
# ---------------------------------------------------------------------------
print('\n' + '='*60)
print('HYPOTHESIS ASSESSMENT (from field analysis)')
print('='*60)
for r in records:
    vol_ref = records[0]['vol_phi'] if records else 0
    print(f"\n  t = {r['t']:.4f}")
    print(f"    |∇φ|_mean = {r['grad_mean']:.3f}  (theory {GRAD_PHI_THEORY:.3f}, "
          f"dev={r['grad_mean']-GRAD_PHI_THEORY:+.3f})")
    print(f"    r_std     = {r['sph_std']:.4f}  (0=sphere; large → deformation [Hyp A])")
    print(f"    ∫φ dV     = {r['vol_phi']:.5f}  (drift from t0: {r['vol_phi']-vol_ref:+.5f})")
    print(f"    n̂ angle  = {r['align_mean']:.2f}°  (0°=radial alignment)")
    print(f"    κ_rms     = {r['kappa_rms']:.3f}  (sphere κ = {KAPPA_SPHERE:.3f})")
    print(f"    tanh_L2   = {r['tanh_l2']:.4f}  (0=perfect CDI profile)")
    print(f"    width     = {r['width_mean']:.4f}  (theory 4ε = {THICKNESS_THEORY:.3f})")
