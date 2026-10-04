"""Checks the fluid statistics averaged in one or two directions against
exact averages of polynomial fields.

The user file sets u = x + 2y + 3z, v = xyz, w = x^2 y + z^2 x + 1 and
p = x^2 + y^2 + z^2 at every step on a frozen fluid, so that the time
averages are the fields themselves and every statistic is a polynomial
whose averages are computed here by Gauss-Legendre quadrature. With
polynomial order 5 the spectral element quadrature of every product (degree
up to 8 in one direction) is exact as well. The statistics are accumulated
directly in the averaged space, on the device when one is used, which is
what this test checks: the basic set averaged in x and y and the full set
averaged in z (2D output), and the full set averaged in xz and the basic
set in yz (1D output). Two more components keep the statistics as 3D
fields (`keep_3d_fields`) and average them when writing, which must give
the same results.
"""
import glob
import subprocess
from os import remove
from os.path import isfile, join

import numpy as np

import conftest
from testlib import configure_nprocs, get_genmeshbox, get_makeneko, run_neko

TEST_DIR = join("tests", "test_fluid_stats_avg")
EXTENT = {1: (0.0, 2.0), 2: (0.0, 1.0), 3: (0.0, 4.0)}
# The 2D output is written in single precision, the 1D one in the
# precision of the run. In a single precision run the gradients of the
# fields carry errors of up to 1e-4 relative to the largest values.
TOL_2D = {"dp": 1e-5, "sp": 5e-4}
TOL_1D = {"dp": 1e-9, "sp": 5e-4}


def fields(x, y, z):
    """The fields, their gradients and the 44 statistics in output order."""
    u = x + 2.0 * y + 3.0 * z
    v = x * y * z
    w = x * x * y + z * z * x + 1.0
    p = x * x + y * y + z * z
    du = (1.0 + 0.0 * x, 2.0 + 0.0 * x, 3.0 + 0.0 * x)
    dv = (y * z, x * z, x * y)
    dw = (2.0 * x * y + z * z, x * x, 2.0 * z * x)
    return [p, u, v, w, p * p, u * u, v * v, w * w, u * v, u * w, v * w,
            u**3, v**3, w**3, u * u * v, u * u * w, u * v * v, u * v * w,
            v * v * w, u * w * w, v * w * w, u**4, v**4, w**4, p**3, p**4,
            p * u, p * v, p * w,
            p * du[0], p * du[1], p * du[2], p * dv[0], p * dv[1], p * dv[2],
            p * dw[0], p * dw[1], p * dw[2],
            sum(d * d for d in du), sum(d * d for d in dv),
            sum(d * d for d in dw), sum(a * b for a, b in zip(du, dv)),
            sum(a * b for a, b in zip(du, dw)),
            sum(a * b for a, b in zip(dv, dw))]


def gauss(d, n=8):
    """Gauss-Legendre nodes and weights of the average over direction d."""
    xi, wi = np.polynomial.legendre.leggauss(n)
    lo, hi = EXTENT[d]
    return 0.5 * (hi - lo) * xi + 0.5 * (hi + lo), 0.5 * wi


def expected_2d(d, c1, c2, n_stats):
    """Exact averages in direction d at the written 2D coordinates, which
    are (z, y) for d=1, (x, z) for d=2 and (x, y) for d=3."""
    nodes, weights = gauss(d)
    out = np.zeros((n_stats, c1.size))
    for q, wq in zip(nodes, weights):
        if d == 1:
            vals = fields(q, c2, c1)
        elif d == 2:
            vals = fields(c1, q, c2)
        else:
            vals = fields(c1, c2, q)
        for i in range(n_stats):
            out[i] += wq * vals[i]
    return out


def expected_1d(dirs, coord, n_stats):
    """Exact averages over the two directions `dirs` at the level
    coordinate in the remaining direction."""
    level = ({1, 2, 3} - set(dirs)).pop()
    n1, w1 = gauss(dirs[0])
    n2, w2 = gauss(dirs[1])
    out = np.zeros((n_stats, coord.size))
    for q1, wq1 in zip(n1, w1):
        for q2, wq2 in zip(n2, w2):
            xyz = {level: coord, dirs[0]: q1, dirs[1]: q2}
            vals = fields(xyz[1], xyz[2], xyz[3])
            for i in range(n_stats):
                out[i] += wq1 * wq2 * vals[i]
    return out


def read_fld_2d(path, mesh_path):
    """Reads a 2D statistics fld: returns the coordinates and the
    statistics in output order. The mesh is only in the first file of the
    series, `mesh_path`."""
    def read(fn):
        with open(fn, "rb") as f:
            header = f.read(132).decode("ascii", errors="replace").split()
            wdsize = int(header[1])
            lx, ly, lz, nelv = (int(header[i]) for i in range(2, 6))
            rdcode = header[11]
            f.read(4)
            f.read(4 * nelv)
            dtype = "<f4" if wdsize == 4 else "<f8"
            data = np.frombuffer(f.read(), dtype=dtype).astype(np.float64)
        return rdcode, nelv, lx * ly * lz, data

    rdcode, nelv, lxyz, data = read(path)
    pos = 0
    if rdcode.startswith("X"):
        xy = data[pos:pos + 2 * nelv * lxyz].reshape(nelv, 2, lxyz)
        pos += 2 * nelv * lxyz
    else:
        mrdcode, mnelv, mlxyz, mdata = read(mesh_path)
        assert mrdcode.startswith("X") and mnelv == nelv and mlxyz == lxyz
        xy = mdata[:2 * nelv * lxyz].reshape(nelv, 2, lxyz)
    n_scalars = int(rdcode[rdcode.index("S") + 1:])
    uv = data[pos:pos + 2 * nelv * lxyz].reshape(nelv, 2, lxyz)
    pos += 2 * nelv * lxyz
    p = data[pos:pos + nelv * lxyz]
    pos += nelv * lxyz
    t = data[pos:pos + nelv * lxyz]
    pos += nelv * lxyz
    s = []
    for _ in range(n_scalars):
        s.append(data[pos:pos + nelv * lxyz])
        pos += nelv * lxyz
    assert pos == data.size, (pos, data.size)
    # The output puts <u>, <v> in the velocity, <p> in the pressure, <pp>
    # in the temperature, the remaining statistics in the scalars and <w>
    # in the last scalar.
    stats = [p, uv[:, 0, :].ravel(), uv[:, 1, :].ravel(), s[-1], t] + s[:-1]
    return xy[:, 0, :].ravel(), xy[:, 1, :].ravel(), np.array(stats)


def check_2d(name, d, n_stats):
    files = sorted(glob.glob(join(TEST_DIR, f"{name}0*.f0*")))
    assert len(files) >= 1, files
    c1, c2, stats = read_fld_2d(files[-1], files[0])
    assert stats.shape[0] == n_stats, stats.shape
    ref = expected_2d(d, c1, c2, n_stats)
    for i in range(n_stats):
        error = np.abs(stats[i] - ref[i]).max() / (np.abs(ref[i]).max() + 1.0)
        assert error < TOL_2D[conftest.RP], \
            f"{name}: statistic {i + 1} relative error {error}"


def check_1d(name, dirs, n_stats):
    files = sorted(glob.glob(join(TEST_DIR, f"{name}[0-9].csv")))
    assert len(files) == 1, files
    data = np.genfromtxt(files[0], delimiter=",")
    assert data.shape[1] == 2 + n_stats, data.shape
    coord = data[:, 1]
    ref = expected_1d(dirs, coord, n_stats)
    for i in range(n_stats):
        error = np.abs(data[:, 2 + i] - ref[i]).max() / (np.abs(ref[i]).max() + 1.0)
        assert error < TOL_1D[conftest.RP], \
            f"{name}: statistic {i + 1} relative error {error}"


def test_fluid_stats_avg(launcher_script, request, log_file, tmp_path):
    del request, tmp_path

    for path in glob.glob(join(TEST_DIR, "avg_*")):
        if isfile(path):
            remove(path)

    result = subprocess.run(
        [get_makeneko(), join(TEST_DIR, "test_fluid_stats_avg.f90")],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert result.returncode == 0, \
        f"makeneko process failed with exit code {result.returncode}"

    # A non-periodic box with non-cubic elements: 5x3x7 elements on
    # [0, 2] x [0, 1] x [0, 4].
    result = subprocess.run(
        [get_genmeshbox(), "0", "2", "0", "1", "0", "4", "5", "3", "7",
         ".false.", ".false.", ".false."],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert result.returncode == 0, \
        f"genmeshbox process failed with exit code {result.returncode}"

    nprocs = configure_nprocs(3)
    case_file = join(TEST_DIR, "test_fluid_stats_avg.case")
    result = run_neko(launcher_script, nprocs, case_file, "./neko", log_file)
    assert result.returncode == 0, \
        f"neko process failed with exit code {result.returncode}"

    check_2d("avg_x", 1, 11)
    check_2d("avg_y", 2, 11)
    check_2d("avg_z", 3, 44)
    check_1d("avg_xz", (1, 3), 44)
    check_1d("avg_yz", (2, 3), 11)
    check_2d("avg_z3d", 3, 44)
    check_1d("avg_xz3d", (1, 3), 11)
