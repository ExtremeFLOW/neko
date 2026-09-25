# icem2re2 — In-depth guide

This document walks through what happens when you convert an ICEM/Fluent
`.msh` mesh into something Neko can simulate on. It's meant for people who
want to understand the *why* of each step (debugging, contributing, or just
trusting the output), not just the *how*. For a quick-reference summary, see
[README.md](README.md).

## The conversion chain

```
   ICEM/Fluent .msh  ──►  Nek5000 .rea  ──►  .re2  ──►  Neko .nmsh
                       │                 │           │
                       │                 │           └── rea2nbin (Neko binary)
                       │                 └────────────── rea2re2.py + pymech
                       └──────────────────────────────── mshconvert.py
   |─────────────── icem2re2.py wraps these two ─────────|
```

`icem2re2.py` is a CLI wrapper that runs **stage 1 (mshconvert)** and
**stage 2 (rea2re2)** in sequence, hiding the intermediate `.rea` in a temp
directory by default. **Stage 3 (rea2nbin)** is a Neko-shipped tool you run
separately.

You don't strictly need this wrapper — you could call mshconvert and rea2re2
yourself — but the wrapper handles a few mechanical chores: locating files,
parsing JSON config, cleaning up intermediates, validating inputs.

## What each format is and why we go through them

### `.msh` (ANSYS Fluent ASCII mesh)

This is what ICEM CFD exports when you "Output > Fluent V6". It's a text
format consisting of S-expression-style sections like:

```
(0 "comment")
(2 3)                                              ; dimension = 3
(10 (1 1 7c5c 1 3) ( ... node coords ... ))        ; node block
(13 (1 1 4b8a3 2 4) ( ... face -> nodes ... ))     ; face block
(39 (13 wall PERIODIC1)())                         ; zone metadata
(18 ( ... periodic face pairings ... ))            ; optional periodic block
```

Key things to know:

- **All integers are hex.** Zone IDs, node/face counts, indices — everything.
  `(39 (13 wall PERIODIC1))` declares zone 0x13 (decimal 19), not decimal 13.
- The mesh stores **faces**, not elements. Element connectivity gets
  reconstructed from face-cell adjacency by mshconvert.
- Boundary types are encoded as small integer codes on each face record
  (3 = wall, 4 = pressure outlet, 8 = periodic, etc.).

### `.rea` (Nek5000 ASCII mesh + parameters)

The historical Nek5000 input format. ASCII, line-oriented, sections marked
by ALL-CAPS headers like `**MESH DATA**`, `***** CURVED SIDE DATA *****`,
`***** FLUID BOUNDARY CONDITIONS *****`.

For each hex element it writes 6 lines: x/y/z of the 4 bottom corners, then
x/y/z of the 4 top corners. For each face of each element it writes one BC
line. This format is human-inspectable, which is enormously useful for
debugging.

### `.re2` (Nek5000 binary mesh)

Same information as `.rea`, just packed binary. Much smaller and faster to
read. pymech writes this from a `HexaData` Python object.

### `.nmsh` (Neko native binary mesh)

Neko's own binary format, written by `rea2nbin`. This is what Neko actually
loads at simulation start. The conversion does two important things:

1. Remaps face numbering from Nek5000-classic convention to Neko's internal
   (Morton-order) convention via a fixed permutation
   `facet_map = (3, 2, 4, 1, 5, 6)`.
2. Verifies periodic BCs: for every face marked `'P'` it checks that the
   declared partner face's corners actually match this face's corners
   under the periodic translation (tolerance `1e-7`).

This verification is **strict** — if even one corner can't be matched,
`rea2nbin` aborts with `"Cannot find matching periodic point"`. Most issues
in the conversion chain surface here.

## Stage 1: `mshconvert.py` (`.msh` → `.rea`)

Upstream `mshconvert.py` is a Python script by Mikael Mortensen
(https://github.com/mikaem/tools), included here with a Python-3
compatibility port and a security fix (replaced `eval()` with `float()` for
node-coordinate parsing). Apart from that we don't modify it.

What it does:

1. Reads the `.msh` line-by-line, regex-matching section headers.
2. Builds tables: `nodes`, `cell_face_map`, `face_list`, `cell_map`, etc.
3. For each cell, runs `create_cell_face_map_3D` which figures out the
   Nek-standard 1-8 corner ordering by matching faces against their
   neighbours. This is the most subtle bit — it ensures every element is
   laid out consistently regardless of how Fluent ordered its faces.
4. If you pass `periodic_dx={(zone_a, zone_b): (dx, dy, dz), ...}`, runs
   `create_periodic_face_map` which:
   - Walks every node on `zone_a` and finds its partner on `zone_b` at the
     translated position (within `1e-7 * |dx|`).
   - For each face on `zone_a`, finds the face on `zone_b` whose 4 nodes
     all map to its 4 nodes under the translation.
   - Writes the result into `periodic_face_map`, then
     `create_periodic_cell_face_map` converts it to (cell, local_face) →
     (shadow_cell, shadow_face).
5. Writes the `.rea`. Periodic faces get BC letter `'P'` with the shadow
   element + face in the connection slots; other boundary faces get
   whatever letter `bcs[zone_id]` says.

**Gotcha:** mshconvert only supports meshes up to ~1M elements (the format
string for the element+face header doesn't cover larger sizes).

**Note on what's in `bcs.json`:** non-periodic zones only. Periodic zones
listed in `periodic.json` get handled by mshconvert directly and are
auto-labelled `'P'`. If you also list them in `bcs.json` they will be
ignored for the periodic side but you'll see misleading output.

## Stage 2: `rea2re2.py` (`.rea` → `.re2`)

Reads the ASCII `.rea` line by line into a pymech `HexaData` structure,
then asks pymech to dump it as binary `.re2`.

The tricky part is parsing BC lines, because **mshconvert writes the
element+face header with a size-dependent format**:

| `nel`           | format                 | example                         | split() yields |
|-----------------|------------------------|---------------------------------|----------------|
| `< 1000`        | `"{elem:3d}{face:3d}"` | `" 14  3 1.0e+02 3.0e+00 ..."`  | 2 separate tokens |
| `1000–99999`    | `"{elem:5d}{face:1d}"` | `"   143 1.0e+02 3.0e+00 ..."`  | 1 concatenated  |
| `100000–999999` | `"{elem:6d}{face:1d}"` | `"    143 1.0e+02 3.0e+00 ..."` | 1 concatenated  |

So when splitting on whitespace, you get 2 leading tokens for small meshes
and 1 leading token for medium/large. This is a real foot-gun: a fix that
works for a tiny test mesh silently corrupts BC parameters on a big
production mesh. The code detects which regime via `nel < 1000` and skips
the right number of leading tokens before reading the 5 numeric BC params.

(We don't actually need the element/face from the line — they're determined
by the iteration counter `ibc`. We only need to skip them to reach the
numeric params.)

## Stage 3: `rea2nbin` (`.re2` → `.nmsh`)

This is a Neko-shipped Fortran tool. In recent Neko versions (≥ 1.99.3) it
takes a `.re2` as input — older versions accepted `.rea` but that's been
removed. Usage:

```bash
LD_LIBRARY_PATH=/path/to/json-fortran-install/lib \
    /path/to/neko-install/bin/rea2nbin in.re2 out.nmsh
```

The `LD_LIBRARY_PATH` is needed for `libjsonfortran.so.9.0` which is
typically not on the system's default loader path.

At this stage Neko validates the mesh end-to-end. Common errors here:

- **"Cannot find matching periodic point"** — a periodic face's corner
  doesn't match its declared partner's corner under the implied
  translation (within `1e-7`). Usually means the upstream tools wrote
  incorrect partner info (e.g. wrong `(elem, face)` because of a parsing
  bug like the one above), or the mesh genuinely isn't node-exactly
  periodic.
- **"rea is no longer supported"** — your Neko is recent enough to reject
  `.rea`; convert to `.re2` first using `rea2re2.py`.

## End-to-end example: a periodic hill

For a mesh `hill.msh` with:

- zone 12 = inlet, zone 15 = outlet, paired x-periodic (displacement +9 in x)
- zone 13 = top, zone 14 = bottom, paired z-periodic (displacement −4.5 in z)
- zone 16 = top wall, zone 17 = hill (lower wall) — both no-slip

`bcs.json` (only non-periodic zones; both walls share the same letter so
Neko will see them as a single labelled zone — split them if you need
per-wall force tracking):

```json
{
    "16": "W",
    "17": "W"
}
```

`periodic.json`:

```json
[
    {"zones": [13, 14], "displacement": [0.0, 0.0, -4.5]},
    {"zones": [12, 15], "displacement": [9.0, 0.0, 0.0]}
]
```

Run:

```bash
./icem2re2.py hill.msh hill.re2 --bcs bcs.json --periodic periodic.json --keep-rea
LD_LIBRARY_PATH=...lib rea2nbin hill.re2 hill.nmsh
```

Sanity checks before going to the cluster:

```bash
# Number of periodic faces in the .rea (should be 2 * sum(faces per pair))
grep -c "^ P " hill.rea

# None of them should have (0, 0) for the partner — that means unpaired
grep "^ P " hill.rea | awk '$4 == "0.0000000e+00" && $5 == "0.0000000e+00"' | wc -l
# expect 0
```

## How to figure out the periodic displacement

If you don't already know the vector between two periodic zones, the
quickest way is to peek at the node bounding boxes. A small helper using
mshconvert's parsing:

```python
import mshconvert as mc

with open('hill.msh') as f:
    mc.scan_fluent_mesh(f)

for z in (13, 14):
    nids = sorted(set(mc.boundary_nodes[z]))
    pts = mc.nodes[:, [n - 1 for n in nids]]
    print(f"zone {z}: bbox {pts.min(axis=1)} - {pts.max(axis=1)}")
    print(f"         centroid {pts.mean(axis=1)}")

dx = (mc.nodes[:, [n - 1 for n in sorted(set(mc.boundary_nodes[14]))]].mean(axis=1) -
      mc.nodes[:, [n - 1 for n in sorted(set(mc.boundary_nodes[13]))]].mean(axis=1))
print(f"displacement zone 13 -> zone 14: {dx}")
```

For a planar periodic surface, the centroid difference *is* the
displacement to use in `periodic.json`.

## Limitations and known gotchas

- **Only translational periodicity.** Rotational (e.g. axisymmetric)
  periodicity is not supported by the code in mshconvert that this tool
  exposes.
- **Mesh size cap ~1M elements** — set by mshconvert's format strings.
- **Single BC label per wall group.** Two zones with the same letter (e.g.
  HILL=`W`, TOP=`W`) merge into a single labelled zone in the `.nmsh`. If
  you want to track per-wall forces or apply different conditions later,
  give them distinct letters in `bcs.json`.
- **Float comparison is strict.** Both mshconvert (when pairing nodes) and
  Neko (when verifying periodic faces) use `1e-7` absolute tolerance.
  Meshes generated by ICEM are usually fine; meshes from other tools may
  need pre-cleaning if node positions drift on supposedly-periodic surfaces.
- **Boundary type from Fluent.** mshconvert uses the face block's bc-type
  field to decide if a face is wall (3), interior (2), periodic (8/12),
  etc. If your `.msh` mis-declares a periodic surface as a wall (which
  ICEM exports often do), you have to declare the periodicity yourself via
  `periodic.json` — mshconvert won't infer it.
