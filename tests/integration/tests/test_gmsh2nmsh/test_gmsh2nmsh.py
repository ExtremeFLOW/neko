"""Integration tests for the gmsh2nmsh mesh converter.

The fixtures are small Gmsh meshes generated from the .geo files in this
directory, with the commands given at the top of each .geo file. Together they
cover MSH 4.1 and 2.2 in ASCII and binary form, first and second order
hexahedra, and a 2D mesh of second order quadrilaterals.
"""
import math
import re
import struct
import subprocess
from pathlib import Path

import pytest

from testlib import get_genmeshbox, get_gmsh2nmsh, get_neko, which_command

HERE = Path(__file__).resolve().parent

# The nmsh writer stores the vertices of an element in cyclic order. This is
# the symmetric vertex number of each stored vertex.
CYC_TO_SYM = [1, 2, 4, 3, 5, 6, 8, 7]

# End vertices, in symmetric numbering, of the 12 curve edges of a hexahedron
HEX_EDGE = [(1, 2), (2, 4), (4, 3), (3, 1), (5, 6), (6, 8), (8, 7), (7, 5),
            (1, 5), (2, 6), (4, 8), (3, 7)]

# Labeled zones (type 7) only define the element, facet and label fields
LABELED_ZONE = 7
PERIODIC_ZONE = 5


def read_nmsh(path):
    """Decode an .nmsh file into elements, zones and curves."""
    data = Path(path).read_bytes()
    nelv, gdim = struct.unpack_from("<ii", data, 0)
    off = 8
    nverts = 8 if gdim == 3 else 4
    elements = []
    for _ in range(nelv):
        el = struct.unpack_from("<i", data, off)[0]
        off += 4
        verts = []
        for _ in range(nverts):
            vid, x, y, z = struct.unpack_from("<iddd", data, off)
            off += 28
            verts.append((vid, (x, y, z)))
        elements.append((el, verts))

    nzones = struct.unpack_from("<i", data, off)[0]
    off += 4
    zones = []
    for _ in range(nzones):
        zone = struct.unpack_from("<9i", data, off)
        off += 36
        if zone[8] == LABELED_ZONE:
            zone = (zone[0], zone[1], 0, zone[3], 0, 0, 0, 0, LABELED_ZONE)
        zones.append(zone)

    ncurves = struct.unpack_from("<i", data, off)[0]
    off += 4
    curves = []
    for _ in range(ncurves):
        el = struct.unpack_from("<i", data, off)[0]
        curve_data = struct.unpack_from("<60d", data, off + 4)
        curve_type = struct.unpack_from("<12i", data, off + 484)
        off += 532
        curves.append((el, curve_data, curve_type))

    return {"nelv": nelv, "gdim": gdim, "elements": elements,
            "zones": zones, "curves": curves}


def resolve_mesh_checker():
    """Locate mesh_checker on PATH or next to Neko."""
    checker = which_command("mesh_checker")
    if checker is not None:
        return Path(checker).resolve()
    candidate = Path(get_neko()).resolve().with_name("mesh_checker")
    assert candidate.is_file(), "The mesh_checker executable could not be found"
    return candidate


def convert(workdir, msh, *options, nmsh="out.nmsh"):
    """Run gmsh2nmsh on a fixture."""
    return subprocess.run(
        [str(Path(get_gmsh2nmsh()).resolve()), str(HERE / msh), nmsh,
         *options],
        cwd=workdir,
        capture_output=True,
        text=True,
        errors="replace",
    )


def check_mesh(workdir, nmsh):
    """Run mesh_checker and return the report from its size section on."""
    result = subprocess.run(
        [str(resolve_mesh_checker()), nmsh],
        cwd=workdir,
        capture_output=True,
        text=True,
        errors="replace",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout[result.stdout.index("--------------Size"):]


def zone_sizes(report):
    """Number of facets of each labeled zone in a mesh_checker report."""
    return {int(zone): int(n)
            for zone, n in re.findall(r"Zone\s+(\d+):\s+(\d+) faces", report)}


def symmetric_vertices(element):
    """Vertex coordinates of an element keyed by symmetric vertex number."""
    return {CYC_TO_SYM[j]: xyz for j, (_, xyz) in enumerate(element[1])}


def test_labeled_zones(tmp_path):
    """Each physical surface of the box becomes one labeled zone."""
    result = convert(tmp_path, "box_v41.msh")
    assert result.returncode == 0, result.stdout + result.stderr

    report = check_mesh(tmp_path, "out.nmsh")
    assert re.search(r"Number of elements:\s+24\n", report)
    assert zone_sizes(report) == {1: 6, 2: 6, 3: 8, 4: 8, 5: 12, 6: 12}


def test_formats_agree(tmp_path):
    """MSH 4.1 ASCII and MSH 2.2 binary give the same mesh."""
    for msh, nmsh in (("box_v41.msh", "a.nmsh"), ("box_v22_bin.msh", "b.nmsh")):
        result = convert(tmp_path, msh, nmsh=nmsh)
        assert result.returncode == 0, result.stdout + result.stderr

    a = read_nmsh(tmp_path / "a.nmsh")
    b = read_nmsh(tmp_path / "b.nmsh")
    assert a["nelv"] == b["nelv"] == 24
    assert sorted(a["zones"]) == sorted(b["zones"])
    for (ea, va), (eb, vb) in zip(a["elements"], b["elements"]):
        assert ea == eb
        assert [v[0] for v in va] == [v[0] for v in vb]
        for (_, xa), (_, xb) in zip(va, vb):
            assert max(abs(p - q) for p, q in zip(xa, xb)) < 1e-12


def test_periodic_matches_genmeshbox(tmp_path):
    """A box made periodic in x and y matches the same box from genmeshbox."""
    result = convert(tmp_path, "box_v41.msh",
                     "--periodic=xmin:xmax,ymin:ymax", nmsh="gmsh.nmsh")
    assert result.returncode == 0, result.stdout + result.stderr
    report = check_mesh(tmp_path, "gmsh.nmsh")
    assert "Number of periodic faces: 28" in report
    assert zone_sizes(report) == {5: 12, 6: 12}

    ref_dir = tmp_path / "ref"
    ref_dir.mkdir()
    result = subprocess.run(
        [str(Path(get_genmeshbox()).resolve()), "0", "2", "0", "1.5", "0",
         "1", "4", "3", "2", ".true.", ".true.", ".false."],
        cwd=ref_dir,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    assert check_mesh(ref_dir, "box.nmsh") == report


def test_curved_annulus(tmp_path):
    """Curved walls of a second order mesh become midpoint curves."""
    result = convert(tmp_path, "annulus_v41_bin.msh", "--periodic=bottom:top")
    assert result.returncode == 0, result.stdout + result.stderr

    mesh = read_nmsh(tmp_path / "out.nmsh")
    elements = {el[0]: el for el in mesh["elements"]}
    assert len(mesh["curves"]) == 16

    ncurved = 0
    for el, curve_data, curve_type in mesh["curves"]:
        verts = symmetric_vertices(elements[el])
        for k, ctype in enumerate(curve_type):
            if ctype == 0:
                continue
            assert ctype == 4
            ncurved += 1
            mid = curve_data[5 * k:5 * k + 3]
            a, b = verts[HEX_EDGE[k][0]], verts[HEX_EDGE[k][1]]
            # The midpoint and both end points lie on one of the walls, and
            # the midpoint belongs to this edge
            radius = math.hypot(mid[0], mid[1])
            assert min(abs(radius - 0.5), abs(radius - 1.0)) < 1e-12
            assert abs(math.hypot(a[0], a[1]) - radius) < 1e-12
            assert abs(math.hypot(b[0], b[1]) - radius) < 1e-12
            chord_mid = [0.5 * (p + q) for p, q in zip(a, b)]
            assert math.dist(mid, chord_mid) < 0.25 * math.dist(a, b)
    assert ncurved == 32

    report = check_mesh(tmp_path, "out.nmsh")
    assert "Number of periodic faces: 32" in report
    assert zone_sizes(report) == {1: 8, 2: 8}

    result = convert(tmp_path, "annulus_v41_bin.msh", "--linear",
                     nmsh="linear.nmsh")
    assert result.returncode == 0, result.stdout + result.stderr
    assert read_nmsh(tmp_path / "linear.nmsh")["curves"] == []


def test_2d_mesh(tmp_path):
    """A 2D mesh is reoriented, and tags above 20 are renumbered."""
    result = convert(tmp_path, "disk_v22.msh")
    assert result.returncode == 0, result.stdout + result.stderr
    assert re.search(r"Reoriented \.+ 4 left-handed elements", result.stdout)
    assert re.search(r"101\s+8\s+1\s+inner", result.stdout)
    assert re.search(r"102\s+8\s+2\s+outer", result.stdout)

    mesh = read_nmsh(tmp_path / "out.nmsh")
    assert mesh["gdim"] == 2 and mesh["nelv"] == 16
    assert len(mesh["curves"]) == 16
    for element in mesh["elements"]:
        # Stored in cyclic order, so the shoelace area is the element area
        xy = [xyz[:2] for _, xyz in element[1]]
        area = sum(xy[i][0] * xy[i - 3][1] - xy[i - 3][0] * xy[i][1]
                   for i in range(4))
        assert area > 0

    report = check_mesh(tmp_path, "out.nmsh")
    assert zone_sizes(report) == {1: 8, 2: 8}


@pytest.mark.parametrize("periodic, message", [
    ("xmin:ymin", "must match one to one"),
    ("xmin:nothere", "is not a boundary physical group"),
    ("xmin", "expects pairs"),
])
def test_invalid_periodic(tmp_path, periodic, message):
    """Invalid periodic pairs are rejected with a clear message."""
    result = convert(tmp_path, "box_v41.msh", f"--periodic={periodic}")
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
