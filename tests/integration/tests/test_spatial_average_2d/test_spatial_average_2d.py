"""Checks the averaging in one direction (map_2d) of the spatial_average
simcomp against exact averages of polynomial fields.

The user file sets u = x + 2y + 3z, v = xyz, w = x^2 y + z^2 x + 1 and
p = x^2 + y^2 + z^2, whose averages in x, y or z are known in closed form.
The fields are averaged on a non-periodic 5x3x7 box with non-cubic elements
and on the same box with its interior elements rotated, so that elements in
a column have different local orientations, and the 2D output is compared
with the exact averages at the written 2D coordinates.
"""
import glob
import json
import os
import struct
import subprocess
from os.path import join

import numpy as np

from testlib import configure_nprocs, get_genmeshbox, get_makeneko, run_neko

TEST_DIR = join("tests", "test_spatial_average_2d")
EXTENT = {1: (0.0, 2.0), 2: (0.0, 1.0), 3: (0.0, 4.0)}
# The 2D output is written in single precision.
TOLERANCE = 1e-5


def exact_average(field, d, c1, c2):
    """Exact average of `field` in direction d (1=x, 2=y, 3=z) at the 2D
    coordinates written by map_2d: (z, y) for d=1, (x, z) for d=2 and
    (x, y) for d=3."""
    lo, hi = EXTENT[d]
    m1 = 0.5 * (lo + hi)
    m2 = (hi**3 - lo**3) / (3.0 * (hi - lo))
    if d == 1:
        x, x2, y, y2, z, z2 = m1, m2, c2, c2**2, c1, c1**2
    elif d == 2:
        x, x2, y, y2, z, z2 = c1, c1**2, m1, m2, c2, c2**2
    else:
        x, x2, y, y2, z, z2 = c1, c1**2, c2, c2**2, m1, m2
    if field == "u":
        return x + 2.0 * y + 3.0 * z
    if field == "v":
        return x * y * z
    if field == "w":
        return x2 * y + z2 * x + 1.0
    return x2 + y2 + z2


def read_fld_2d(path, mesh_path):
    """Reads a 2D fld written from the fields [u, v, w, p]. The mesh is only
    stored in the first file of the series, `mesh_path`."""
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
        assert lz == 1, header
        return rdcode, nelv, lx * ly * lz, data

    rdcode, nelv, lxyz, data = read(path)
    assert rdcode in ("XUPS01", "UPS01"), rdcode
    pos = 0
    if rdcode.startswith("X"):
        xy = data[pos:pos + 2 * nelv * lxyz].reshape(nelv, 2, lxyz)
        pos += 2 * nelv * lxyz
    else:
        mrdcode, mnelv, mlxyz, mdata = read(mesh_path)
        assert mrdcode.startswith("X") and mnelv == nelv and mlxyz == lxyz
        xy = mdata[:2 * nelv * lxyz].reshape(nelv, 2, lxyz)
    # Velocity (u, v), pressure (p) and, as the last scalar, w.
    uv = data[pos:pos + 2 * nelv * lxyz].reshape(nelv, 2, lxyz)
    pos += 2 * nelv * lxyz
    p = data[pos:pos + nelv * lxyz]
    pos += nelv * lxyz
    w = data[pos:pos + nelv * lxyz]
    pos += nelv * lxyz
    assert pos == data.size, (pos, data.size)
    return (xy[:, 0, :].ravel(), xy[:, 1, :].ravel(),
            {"u": uv[:, 0, :].ravel(), "v": uv[:, 1, :].ravel(),
             "w": w, "p": p})


# Vertex slots of an nmsh hex element in the cyclic ordering: 1-4 counter-
# clockwise around the bottom face, 5-8 the corresponding top vertices.
_HEX = struct.Struct("<i" + "i3d" * 8)
_POS = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                 (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)], dtype=float)
_RX = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
_RY = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
_RZ = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
_ROTATIONS = [_RZ, _RZ @ _RZ, _RZ @ _RZ @ _RZ, _RX @ _RX, _RY @ _RY,
              _RX, _RY, _RX @ _RZ]


def _rotation_perm(R):
    """perm[k] = old vertex slot that becomes slot k when the element is
    relabelled by the proper rotation R about its centre."""
    c = np.array([0.5, 0.5, 0.5])
    perm = []
    for k in range(8):
        found = [m for m in range(8)
                 if np.allclose(R @ (_POS[m] - c) + c, _POS[k])]
        assert len(found) == 1
        perm.append(found[0])
    return perm


def rotate_interior_elements(src, dst):
    """Writes a copy of the nmsh mesh `src` in which the elements that do
    not touch the bounding box are relabelled with different rotations, so
    that their local orientation differs from their neighbours'."""
    with open(src, "rb") as f:
        head = f.read(8)
        nelv = struct.unpack("<ii", head)[0]
        records = [bytearray(f.read(_HEX.size)) for _ in range(nelv)]
        tail = f.read()
    verts = []
    for rec in records:
        r = _HEX.unpack(rec)
        verts.append(np.array([r[2 + 4 * k:5 + 4 * k] for k in range(8)]))
    allv = np.concatenate(verts)
    lo, hi = allv.min(axis=0), allv.max(axis=0)
    rotated = 0
    with open(dst, "wb") as f:
        f.write(head)
        for rec, v in zip(records, verts):
            inner = (v.min(axis=0) > lo + 1e-12).all() and \
                (v.max(axis=0) < hi - 1e-12).all()
            if inner:
                r = _HEX.unpack(rec)
                perm = _rotation_perm(_ROTATIONS[rotated % len(_ROTATIONS)])
                vals = [r[0]]
                for k in range(8):
                    vals += list(r[1 + 4 * perm[k]:5 + 4 * perm[k]])
                rec = _HEX.pack(*vals)
                rotated += 1
            f.write(rec)
        f.write(tail)
    return rotated


def check_outputs(prefix):
    for d, name in ((1, "x"), (2, "y"), (3, "z")):
        files = sorted(glob.glob(join(TEST_DIR, f"{prefix}sa_{name}0*.f0*")))
        assert len(files) >= 2, files
        c1, c2, values = read_fld_2d(files[-1], files[0])
        for field in ("u", "v", "w", "p"):
            ref = exact_average(field, d, c1, c2)
            error = np.abs(values[field] - ref).max() / (np.abs(ref).max() + 1.0)
            assert error < TOLERANCE, \
                f"{prefix}avg_{name} of {field}: relative error {error}"


def test_spatial_average_2d(launcher_script, request, log_file, tmp_path):
    del request

    # Remove stale outputs (sa_x0.f00000, ..., rot_sa_x0.nek5000).
    for pattern in ("sa_?0*", "rot_sa_?0*"):
        for path in glob.glob(join(TEST_DIR, pattern)):
            os.remove(path)

    result = subprocess.run(
        [get_makeneko(), join(TEST_DIR, "test_spatial_average_2d.f90")],
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
    case_file = join(TEST_DIR, "test_spatial_average_2d.case")
    result = run_neko(launcher_script, nprocs, case_file, "./neko", log_file)
    assert result.returncode == 0, \
        f"neko process failed with exit code {result.returncode}"
    check_outputs("")

    # The same box with its 15 interior elements rotated.
    rotated = rotate_interior_elements("box.nmsh", "box_rot.nmsh")
    assert rotated == 15, rotated
    with open(case_file) as f:
        case = json.load(f)
    case["case"]["mesh_file"] = "box_rot.nmsh"
    for comp in case["case"]["simulation_components"]:
        comp["output_filename"] = "rot_" + comp["output_filename"]
    rot_case_file = str(tmp_path / "test_spatial_average_2d_rot.case")
    with open(rot_case_file, "w") as f:
        json.dump(case, f, indent=2)
    result = run_neko(launcher_script, nprocs, rot_case_file, "./neko",
                      log_file)
    assert result.returncode == 0, \
        f"neko process failed with exit code {result.returncode}"
    check_outputs("rot_")
