"""Error norms and exact solutions from Saini et al. (2026), JCP 561, 114961.

**This file is a snapshot**, trimmed to what the three `visualize.ipynb`
notebooks in this repo actually use. The living version is
`../neko-multiphase/examples/saini_benchmarks/norms.py`; if the two drift, that
repo is the source of truth. Same convention as
`references/saini_2026_test_cases.md`.

See `references/saini_2026_test_cases.md` section 2 for what each norm measures
and why `E_r` and `E_v` decouple.

Naming, fixed project-wide and deliberately *not* Saini's: **`phi` is the phase
field** (0-1, tanh profile) and **`psi` is the signed-distance field**. Saini use
those two symbols the other way round, so their equation numbers read transposed
against this code. See `CDI_METHOD.md` section 1.

Every case here is quasi-2D: one element in z, with a z-invariant field.
`coef.B` is a 3D mass matrix, so each integral is divided by the slab thickness
to recover the area integral the paper's definitions assume.
"""

import numpy as np
import matplotlib.tri as mtri
from matplotlib.figure import Figure


# --------------------------------------------------------------------------
# integration
# --------------------------------------------------------------------------

def slab_thickness(msh):
    """z-extent of the one-element-thick slab."""
    return float(msh.z.max() - msh.z.min())


def area_integral(f, coef, msh):
    """int f dA over the 2D domain, from the 3D GLL mass matrix."""
    return float((f * coef.B).sum()) / slab_thickness(msh)


# --------------------------------------------------------------------------
# Saini Eqs. (79)-(81)
# --------------------------------------------------------------------------

def E_r(phi, phi_e, coef, msh):
    """Relative L1 norm, Eq. (79). Dominated by over/undershoots at the
    interface -- this is the norm that sees dispersion error."""
    return (area_integral(np.abs(phi - phi_e), coef, msh)
            / area_integral(phi_e, coef, msh))


def E_v(phi, phi_0, phi_e, coef, msh):
    """Volume (mass) conservation error, Eq. (80). Signed."""
    return ((area_integral(phi, coef, msh) - area_integral(phi_0, coef, msh))
            / area_integral(phi_e, coef, msh))


def E_s(phi, phi_e, coef, msh, perimeter):
    """Shape error, Eq. (81): how much enclosed area crossed to the wrong side
    of the interface, normalised by 2 * perimeter * exact area.

    The Heaviside is sharp, as written in the paper, so this is a crude
    quadrature of a discontinuous integrand at the GLL points. Kept verbatim so
    the values stay comparable with theirs; it is also why the paper reports
    non-monotone behaviour for this norm.
    """
    num = area_integral(np.abs(_heav(phi - 0.5) - _heav(phi_e - 0.5)), coef, msh)
    den = 2.0 * perimeter * area_integral(_heav(phi_e - 0.5), coef, msh)
    return num / den


def _heav(g):
    return np.where(g >= 0.0, 1.0, 0.0)


# --------------------------------------------------------------------------
# interface traces
# --------------------------------------------------------------------------

def l_avg(phi, coef, msh, perimeter, lo=0.05, hi=0.95):
    """Mean interface thickness: area of the band {lo <= phi <= hi} over the
    length of the 0.5 isocontour.

    Saini's Eq. (87) puts their diffuse delta function in the denominator. We do
    not have that field, so the isocontour length is used instead -- the same
    quantity (interface perimeter) by a route we can compute. The band covers
    90% of the diffused interface, as in their definition.
    """
    band = np.where((phi >= lo) & (phi <= hi), 1.0, 0.0)
    return area_integral(band, coef, msh) / perimeter


def contour_length(msh, f, level=0.5, iz=None):
    """Length of the `level` isocontour of `f`, from the GLL points directly.

    Uses a Delaunay triangulation of the point cloud (the domain is a convex
    box, so this is well posed) rather than resampling onto a uniform grid.
    """
    if iz is None:
        iz = msh.z.shape[1] // 2
    x = msh.x[:, iz, :, :].ravel()
    y = msh.y[:, iz, :, :].ravel()
    v = f[:, iz, :, :].ravel()

    # Figure() rather than pyplot: no backend, no open window, nothing to close.
    ax = Figure().add_subplot(111)
    cset = ax.tricontour(mtri.Triangulation(x, y), v, levels=[level])

    total = 0.0
    for path in cset.get_paths():
        seg = np.diff(path.vertices, axis=0)
        total += float(np.hypot(seg[:, 0], seg[:, 1]).sum())
    return total


def n_unbounded(phi):
    """Boundedness violations: nodes with phi outside [0, 1]."""
    return int(np.count_nonzero((phi < 0.0) | (phi > 1.0)))


# --------------------------------------------------------------------------
# advecting_slab_1d
# --------------------------------------------------------------------------

SLAB_HALFWIDTH = 0.2


def slab_distance(x, t=0.0, u=1.0, center=0.5, halfwidth=SLAB_HALFWIDTH):
    """Signed distance to the slab, positive inside, advected to time `t`.

    The domain is periodic on [0, 1] and u = 1, so the slab returns to its
    initial position at every integer t -- including the t = 20 endpoint the
    error norms are quoted at.
    """
    dx = np.abs((x - (center + u * t)) % 1.0)
    return halfwidth - np.minimum(dx, 1.0 - dx)


def slab_phase(x, eps, t=0.0, u=1.0, center=0.5, halfwidth=SLAB_HALFWIDTH):
    return 0.5 * (1.0 + np.tanh(slab_distance(x, t, u, center, halfwidth)
                                / (2.0 * eps)))


# --------------------------------------------------------------------------
# zalesak_disk -- slotted disk, Saini Eq. (77)
# --------------------------------------------------------------------------

SLOT_CX, SLOT_CY, SLOT_R = 0.5, 0.75, 0.15
SLOT_HALFWIDTH, SLOT_TOP = 0.025, 0.85


def slotted_disk_distance(x, y):
    """Signed distance to the slotted disk, positive inside (their Eq. 77).

    Same construction as the .f90 initial condition; duplicated here so the
    notebook can build the exact solution without reading fields.
    """
    d0 = np.sqrt((x - SLOT_CX) ** 2 + (y - SLOT_CY) ** 2) - SLOT_R
    d1 = SLOT_HALFWIDTH - np.abs(x - SLOT_CX)
    d2 = SLOT_TOP - y
    return -np.maximum(d0, np.minimum(d1, d2))


def slotted_disk_distance_periodic(x, y):
    """Eq. (77) made periodic on [0,1]^2: the nearest of the nine images.

    `slotted_disk_distance` is the infinite-plane distance and the disk is
    centred at y = 0.75, so that field jumps by 0.5 across the y = 0/1 seam of a
    periodic mesh. That is invisible in `phi` (both sides saturate to 0) and
    fatal to a transported `psi`, whose gather-scatter would average across the
    jump. Use this one wherever `psi` itself is carried as a field.
    """
    best = None
    for dx in (-1.0, 0.0, 1.0):
        for dy in (-1.0, 0.0, 1.0):
            cand = slotted_disk_distance(x + dx, y + dy)
            best = cand if best is None else np.where(
                np.abs(cand) < np.abs(best), cand, best)
    return best


def slotted_disk_phase(x, y, eps):
    return 0.5 * (1.0 + np.tanh(slotted_disk_distance(x, y) / (2.0 * eps)))


def slotted_disk_phase_at(x, y, t, eps, omega=np.pi):
    """The exact solution of the Zalesak case at time `t`.

    The flow u = pi(0.5 - y), v = pi(x - 0.5) is a rigid counter-clockwise
    rotation at omega = pi about the box centre, so the exact phase field is
    just the initial shape carried round: evaluate it at the pre-image of each
    point. Comparing against the *unrotated* disk instead is the easy mistake
    -- it is only correct at whole rotations (t = 2, 4, ... 20).
    """
    xr, yr = rotate_about(x, y, omega * t)
    return slotted_disk_phase(xr, yr, eps)


def rotate_about(x, y, theta, cx=0.5, cy=0.5):
    """Rigid rotation by `theta` about (cx, cy) -- the exact solution map for
    Zalesak's solid-body rotation u = pi(0.5 - y), v = pi(x - 0.5)."""
    c, s = np.cos(theta), np.sin(theta)
    xr, yr = x - cx, y - cy
    return cx + c * xr + s * yr, cy - s * xr + c * yr


# --------------------------------------------------------------------------
# rider_kothe -- plain disk, same position and radius, no slot
# --------------------------------------------------------------------------

def disk_distance(x, y):
    """Signed distance to the plain disk, positive inside."""
    return SLOT_R - np.sqrt((x - SLOT_CX) ** 2 + (y - SLOT_CY) ** 2)


def disk_distance_periodic(x, y):
    """Nearest of the nine images -- the disk sits 0.10 from the y = 1 seam, so
    `slotted_disk_distance_periodic`'s reasoning applies here unchanged."""
    best = None
    for dx in (-1.0, 0.0, 1.0):
        for dy in (-1.0, 0.0, 1.0):
            cand = disk_distance(x + dx, y + dy)
            best = cand if best is None else np.where(
                np.abs(cand) < np.abs(best), cand, best)
    return best


def disk_phase(x, y, eps):
    return 0.5 * (1.0 + np.tanh(disk_distance(x, y) / (2.0 * eps)))
