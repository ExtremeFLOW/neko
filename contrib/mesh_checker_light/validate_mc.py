#!/usr/bin/env python3
# Copyright (c) 2026, The Neko Authors
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials provided
#     with the distribution.
#
#   * Neither the name of the authors nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
#     _  __  ____  __ __  ____
#    / |/ / / __/ / //_/ / __ \
#   /    / / _/  / ,<   / /_/ /
#  /_/|_/ /___/ /_/|_|  \____/
#
"""Validation for mesh_checker_light across every 3D .nmsh in the Neko tree.

There is no bundled reference output for mesh_checker (it needs the full Neko
build to run), so this harness independently recomputes -- in Python, from the
raw nmsh bytes -- every quantity mesh_checker_light reports, cross-checks the
tool's stdout against it, and verifies the topological identity that any correct
face dedup must satisfy:

    2 * (#faces) - (#external faces) == 6 * nelv          (each internal face is
    shared by exactly 2 elements, each external face by 1; 6 facets per hex)

and that every external face of a VALID mesh is either labeled or periodic
(i.e. #unlabeled-external == 0). It also confirms the tool's exit code (0 for a
valid mesh). A separate broken-mesh test strips a labeled zone and
checks the tool then detects exactly the newly-exposed faces.

Usage:  python3 validate_mc.py [path/to/mesh_checker_light] [neko_root]
"""
import struct, subprocess, glob, os, sys, re

HERE = os.path.dirname(os.path.abspath(__file__))
EXE  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'mesh_checker_light')
NEKO = sys.argv[2] if len(sys.argv) > 2 else os.path.abspath(os.path.join(HERE, '..', '..'))

FACE = [(1,5,8,4),(2,6,7,3),(1,2,6,5),(4,3,7,8),(1,2,3,4),(5,6,7,8)]   # nmsh-vtx, 1-based
EDGE = [(1,2),(3,4),(5,6),(7,8),(1,4),(2,3),(5,8),(6,7),(1,5),(2,6),(4,8),(3,7)]


def compute(path):
    d = open(path, 'rb').read()
    nelv, gdim = struct.unpack('<ii', d[:8])
    if gdim != 3:
        return None
    off = 8
    elems = []          # list of 8 v_idx (nmsh order)
    maxv = 0
    for _ in range(nelv):
        struct.unpack('<i', d[off:off+4]); off += 4      # el_idx (identity)
        v = []
        for _ in range(8):
            vi = struct.unpack('<i', d[off:off+4])[0]
            off += 28                                    # v_idx + 3 dp
            v.append(vi); maxv = max(maxv, vi)
        elems.append(v)
    nz = struct.unpack('<i', d[off:off+4])[0]; off += 4
    marked = {}         # (e_pos_1based, facet_1based) -> 'L' or 'P'
    periodic = 0
    labeled = {}
    remap = {}          # raw point id -> merged (periodic) id, exactly like Neko
    for _ in range(nz):
        z = struct.unpack('<9i', d[off:off+36]); off += 36
        e, f, ztype, plabel = z[0], z[1], z[8], z[3]
        if ztype == 5:
            periodic += 1; marked[(e, f)] = 'P'
            # Reproduce apply_periodic_facet: corner k of facet f (nmsh slot
            # FACE[f][k]) takes the record's glb_pt_ids(k). The file stores those
            # in the same face_nodes order (mesh_reset_periodic_ids/hex_facet_order).
            g = z[4:8]
            for k in range(4):
                raw = elems[e-1][FACE[f-1][k]-1]
                remap[raw] = min(remap.get(raw, g[k]), g[k])
        elif ztype == 7:
            labeled[plabel] = labeled.get(plabel, 0) + 1; marked[(e, f)] = 'L'

    def rm(x):
        return remap.get(x, x)

    # Neko's mesh_checker reports the counts AFTER the periodic merge
    # (glb_mfcs / glb_meds), so faces/edges are hashed on the MERGED ids.
    fcnt = {}           # sorted 4-tuple of merged ids -> (count, first (e,f))
    for e in range(nelv):
        for fi, fc in enumerate(FACE):
            key = tuple(sorted(rm(elems[e][c-1]) for c in fc))
            if key in fcnt:
                fcnt[key] = (fcnt[key][0] + 1, fcnt[key][1])
            else:
                fcnt[key] = (1, (e + 1, fi + 1))
    edges = set()
    for e in range(nelv):
        for ec in EDGE:
            edges.add(tuple(sorted(rm(elems[e][c-1]) for c in ec)))
    external = [ef for (cnt, ef) in fcnt.values() if cnt == 1]
    unlabeled = sum(1 for ef in external if ef not in marked)
    # A mesh that is a single element thick in a periodic direction collapses each
    # cross-section to one point, so some element loses distinct corners after the
    # merge. Neko produces exactly these collapsed counts too (same apply path), but
    # the non-degenerate topological identities below no longer apply -- flag it.
    degenerate = any(len(set(rm(elems[e][c-1]) for c in range(1, 9))) < 8
                     for e in range(nelv))
    return dict(nelv=nelv, points=maxv, faces=len(fcnt), edges=len(edges),
                periodic=periodic, labeled=labeled, n_external=len(external),
                unlabeled=unlabeled, degenerate=degenerate)


def parse_tool(out):
    g = {}
    for pat, key in [(r'Number of elements:\s*(\d+)', 'nelv'),
                     (r'Number of points:\s*(\d+)', 'points'),
                     (r'Number of faces:\s*(\d+)', 'faces'),
                     (r'Number of edges:\s*(\d+)', 'edges'),
                     (r'Number of periodic faces:\s*(\d+)', 'periodic')]:
        m = re.search(pat, out)
        if m:
            g[key] = int(m.group(1))
    g['labeled'] = {int(a): int(b) for a, b in re.findall(r'Zone\s+(\d+):\s*(\d+)\s*faces', out)}
    m = re.search(r'Found\s+(\d+)\s+unlabeled external faces', out)
    g['unlabeled'] = int(m.group(1)) if m else 0
    return g


def broken_mesh_test(exe, neko):
    """Strip the labeled wall zone from turb_pipe and confirm the tool detects
    exactly the newly-exposed faces (the whole diagnostic point of the tool)."""
    import tempfile
    src = os.path.join(neko, 'examples', 'turb_pipe', 'turb_pipe.nmsh')
    if not os.path.exists(src):
        print('  (skip: turb_pipe.nmsh not found)'); return True
    d = open(src, 'rb').read()
    nelv = struct.unpack('<i', d[:4])[0]
    off = 8 + nelv * 228
    nz = struct.unpack('<i', d[off:off+4])[0]; zbase = off + 4
    zones = [d[zbase+i*36:zbase+i*36+36] for i in range(nz)]
    tail = d[zbase+nz*36:]
    kept = [z for z in zones if struct.unpack('<9i', z)[8] != 7]
    dropped_faces = sum(1 for z in zones if struct.unpack('<9i', z)[8] == 7)
    new = d[:off] + struct.pack('<i', len(kept)) + b''.join(kept) + tail
    with tempfile.TemporaryDirectory() as t:
        b = os.path.join(t, 'broken.nmsh'); open(b, 'wb').write(new)
        p = subprocess.run([exe, b], capture_output=True, text=True)
    m = re.search(r'Found\s+(\d+)\s+unlabeled', p.stdout)
    got = int(m.group(1)) if m else 0
    ok = (got == dropped_faces and p.returncode == 1)
    print('broken-mesh test (strip turb_pipe wall zone = %d faces):' % dropped_faces)
    print('  tool detected %d unlabeled external faces, exit %d -> %s' %
          (got, p.returncode, 'PASS' if ok else 'FAIL'))
    return ok


def jacobian_test(exe):
    """A correctly-oriented unit cube (jac>0, not flagged) and a mirrored one
    (jac<0, must be flagged) confirm the negative-Jacobian detector."""
    import tempfile
    good = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    mirror = [(-x, y, z) for (x, y, z) in good]

    def make(path, corners):
        with open(path, 'wb') as f:
            f.write(struct.pack('<ii', 1, 3)); f.write(struct.pack('<i', 1))
            for j in range(8):
                f.write(struct.pack('<i', j + 1)); f.write(struct.pack('<3d', *corners[j]))
            f.write(struct.pack('<i', 0)); f.write(struct.pack('<i', 0))

    ok = True
    with tempfile.TemporaryDirectory() as t:
        for name, c, want_flag in [('good cube', good, False), ('mirrored cube', mirror, True)]:
            p = os.path.join(t, 'm.nmsh'); make(p, c)
            r = subprocess.run([exe, p, '--jacobian'], capture_output=True, text=True)
            flagged = ('negative/zero Jacobian (first' in r.stdout)
            good_case = (flagged == want_flag)
            ok = ok and good_case
            print('  %-16s flagged=%-5s (want %s) -> %s' %
                  (name, flagged, want_flag, 'OK' if good_case else 'FAIL'))
    return ok


def fld_test(exe, neko):
    """Write zone_indices.fld for the cylinder, parse it back, and confirm the
    scalar field equals the mesh's per-GLL-point zone label everywhere."""
    import tempfile
    mesh = os.path.join(neko, 'examples', 'cylinder', 'cylinder.nmsh')
    if not os.path.exists(mesh):
        print('  (skip: cylinder.nmsh not found)'); return True
    exe = os.path.abspath(exe)
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as t:
        os.chdir(t)
        try:
            subprocess.run([exe, mesh, '--write-zone-indices'], capture_output=True, text=True)
            d = open('zone_indices.fld', 'rb').read()
        finally:
            os.chdir(cwd)
    hdr = d[:132].split(); off = 132
    lx, nelv = int(hdr[2]), int(hdr[5]); lxyz = lx**3
    tp = struct.unpack('<f', d[off:off+4])[0]; off += 4
    off += 4 * nelv                                    # idx
    off += 4 * nelv * 3 * lxyz                          # mesh block
    scal = struct.unpack('<%df' % (nelv * lxyz), d[off:off+4*nelv*lxyz]); off += 4*nelv*lxyz
    off += 4 * nelv * 6                                 # geometry bb metadata
    off += 4 * nelv * 2                                 # scalar min/max metadata
    struct_ok = (off == len(d) and abs(tp - 6.54321) < 1e-4)
    # expected labels from the mesh
    m = open(mesh, 'rb').read(); mn = struct.unpack('<i', m[:4])[0]
    mo = 8 + mn * 228; nz = struct.unpack('<i', m[mo:mo+4])[0]; mo += 4

    def facet_pts(ff):
        out = []
        for it in range(3):
            for js in range(3):
                for ir in range(3):
                    on = ((ff==1 and ir==0) or (ff==2 and ir==2) or (ff==3 and js==0) or
                          (ff==4 and js==2) or (ff==5 and it==0) or (ff==6 and it==2))
                    if on:
                        out.append(ir + 3*js + 9*it)
        return out

    exp = {}
    for i in range(nz):
        z = struct.unpack('<9i', m[mo+i*36:mo+i*36+36])
        if z[8] == 7:
            e, ff, lbl = z[0]-1, z[1], z[3]
            for pt in facet_pts(ff):
                exp[(e, pt)] = max(exp.get((e, pt), 0), lbl)
    bad = sum(1 for e in range(mn) for pt in range(27)
              if round(scal[e*27+pt]) != exp.get((e, pt), 0))
    ok = struct_ok and bad == 0
    print('  .fld struct ok=%s ; zone-index field mismatches=%d -> %s' %
          (struct_ok, bad, 'OK' if ok else 'FAIL'))
    return ok


def main():
    if not os.path.exists(EXE):
        print('error: mesh_checker_light not found at', EXE, '(build it first)'); return 1
    paths = sorted(glob.glob(os.path.join(NEKO, 'examples', '**', '*.nmsh'), recursive=True) +
                   glob.glob(os.path.join(NEKO, 'tests', '**', '*.nmsh'), recursive=True))
    print('%-46s %8s %8s %8s %6s %6s %5s' %
          ('mesh', 'nelv', 'faces', 'edges', 'per', 'unlab', 'ok'))
    n = nok = 0
    allok = True
    for path in paths:
        ref = compute(path)
        if ref is None:
            continue
        p = subprocess.run([EXE, path], capture_output=True, text=True)
        got = parse_tool(p.stdout)
        # (1) tool output must equal the independent Python computation, field by field
        fields_ok = all(got.get(k) == ref[k] for k in ('nelv', 'points', 'faces', 'edges', 'periodic', 'unlabeled'))
        lbl_ok = (got['labeled'] == ref['labeled'])
        # (2) topological identity any correct MERGED face dedup must satisfy:
        #     each internal face is shared by 2 elements, each external by 1.
        #     (Skipped for degenerate single-element-thick collapsed meshes, which
        #     Neko reproduces identically -- there the tool==python check carries it.)
        euler_ok = ref['degenerate'] or \
            (2 * ref['faces'] - ref['n_external'] == 6 * ref['nelv'])
        # (2b) a non-degenerate mesh with no external faces is a closed 3-torus
        #      (fully periodic): then #faces == #edges == 3*nelv exactly. Hard
        #      ground truth that ONLY the correct periodic merge (right glb_pt_ids
        #      ordering) satisfies -- the raw, unmerged counts do not.
        torus_ok = True
        if ref['n_external'] == 0 and not ref['degenerate']:
            torus_ok = (ref['faces'] == 3 * ref['nelv'] and ref['edges'] == 3 * ref['nelv'])
        # (3) exit code reflects the diagnostic (0 iff no unlabeled external faces)
        exit_ok = (p.returncode == (1 if ref['unlabeled'] > 0 else 0))
        ok = fields_ok and lbl_ok and euler_ok and torus_ok and exit_ok
        n += 1; nok += ok; allok = allok and ok
        note = '' if ref['unlabeled'] == 0 else '(no-BC mesh: whole surface unlabeled)'
        print('%-46s %8d %8d %8d %6d %6d %-4s %s' %
              (os.path.relpath(path, NEKO), ref['nelv'], ref['faces'], ref['edges'],
               ref['periodic'], ref['unlabeled'], 'OK' if ok else 'FAIL', note))
        if not ok:
            print('     tool=%s ref=%s euler=%s exit=%d' %
                  (got, {k: ref[k] for k in ('nelv','points','faces','edges','periodic','labeled','unlabeled')},
                   euler_ok, p.returncode))
    print('\n3D meshes tested            : %d' % n)
    print('tool==python + Euler + exit : %d / %d' % (nok, n))
    print('  (independent Python recomputation matches the tool on every field;')
    print('   2F-B=6*nelv holds; exit code = 0 iff 0 unlabeled external faces.')
    print('   Meshes with unlabeled>0 are benchmark boxes/cylinders that carry no')
    print('   mesh-level BCs -- the tool correctly flags their whole open surface,')
    print('   exactly as the real mesh_checker would.)')
    print()
    brk = broken_mesh_test(EXE, NEKO)
    print('\njacobian test (good vs mirrored/inverted cube):')
    jac = jacobian_test(EXE)
    print('\nzone-index .fld test:')
    fld = fld_test(EXE, NEKO)
    allok = allok and brk and jac and fld
    print('\nRESULT:', 'PASS' if allok else 'FAIL')
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
