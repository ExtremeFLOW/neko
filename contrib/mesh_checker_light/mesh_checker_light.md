```text
     _  __  ____  __ __  ____
    / |/ / / __/ / //_/ / __ \
   /    / / _/  / ,<   / /_/ /
  /_/|_/ /___/ /_/|_|  \____/
```

# mesh_checker_light

A standalone, low-memory diagnostics tool for Neko `.nmsh` meshes — a lightweight,
CPU-only alternative to the **default** checks of `contrib/mesh_checker`.

Needs only a Fortran compiler (no MPI, no Neko library):

```sh
make                                   # gfortran -O2 -o mesh_checker_light mesh_checker_light.f90
./mesh_checker_light mesh.nmsh
./mesh_checker_light mesh.nmsh --jacobian              # + negative/zero Jacobian check
./mesh_checker_light mesh.nmsh --write-zone-indices    # + write zone_indices.fld (implies --jacobian)
make check                             # cross-check vs independent Python on every 3D mesh
```

> `make check` reads the meshes shipped with Neko, so it expects this folder to
> live inside a Neko source checkout (e.g. dropped into `contrib/`, so
> `../../examples` and `../../tests` resolve). To run it from elsewhere, point the
> validator at your Neko tree directly:
> `python3 validate_mc.py ./mesh_checker_light /path/to/neko`. Checking your own
> mesh (`./mesh_checker_light mesh.nmsh`) needs no Neko tree at all.

It reports — exactly like `mesh_checker` — the element / point / face / edge
counts, the bounding box, the number of periodic faces, the labeled-zone face
counts, and the key diagnostic: the number of **unlabeled external faces** (an
external face that is neither a labeled zone nor periodic — must be `0` for a
mesh whose boundaries are fully specified). It exits non-zero if any are found.
With `--jacobian` it also flags **negative/zero Jacobians** (inverted/tangled
elements); with `--write-zone-indices` it also writes the `zone_indices.fld`
visualization field. Both are **streamed element-by-element** — O(1) extra memory,
no `dofmap`/`coef`.

## Why it exists

Neko's `contrib/mesh_checker` reads the `.nmsh` into a full `mesh_t` with
connectivity generation on, because it needs the face/edge neighbour information
for the unlabeled-face check and a `dofmap` for the bounding box. Building that
goes through Neko's general-purpose, production hash tables — the same robust,
well-tested structures reused all over the code base. They're sized generously
and store polymorphic keys, which is exactly what you want inside the solver, but
it means a lot of memory per element when all you want is a quick sanity check of
a large mesh.

`mesh_checker_light` is a lightweight reimplementation of just those checks. It
works directly from the deduplicated point ids already stored in the `.nmsh`,
builds a couple of small typed integer hashes instead of the general tables, and
does **no Neko initialisation** — giving the same answers in a small fraction of
the memory.

One practical note: `mesh_checker` is MPI-parallel, so its footprint divides
across ranks — with enough ranks per-rank memory is fine. The awkward case is
running a quick check on **one or a few ranks** for a very large mesh, where
you'd otherwise have to spin up a big MPI job just to confirm the boundaries are
labeled. That's the gap this fills: a **serial** check of a large mesh on a
single node.

A couple of nice side effects of not linking Neko:

- **No GPU rebuild needed.** The `contrib` tools are CPU-only. If your Neko is
  configured for GPUs, you'd normally have to rebuild it for CPU just to run
  them. `mesh_checker_light` needs only a Fortran compiler, so that step goes
  away.
- **Runs anywhere** — a login node, a laptop, a small VM — since there's no MPI
  and no Neko library to carry along.

> **Tested, but please check your own case.** This tool is cross-checked field
> for field against an independent Python recomputation on every 3D mesh shipped
> with Neko, plus topological-identity and broken-mesh tests (see
> [Validation](#validation)). It is still your responsibility to make sure the
> result is correct for *your* mesh — `make check`, or cross-check against a real
> `mesh_checker` run, on a case you trust before relying on it.

## What `mesh_checker_light` does instead

The mesh's `.nmsh` already stores the deduplicated point ids (`v_idx`), so no
point table is needed. The tool reads the element connectivity (8 ids per
element), applies the periodic merge, and builds just **two small typed integer
hashes** (structure-of-arrays, no polymorphism, no per-entry heap):

- a **face hash** keyed by the 4 sorted `v_idx` of each facet, tracking an
  occurrence count (external = seen once, internal = seen twice) — `~17 B/slot`;
- an **edge hash** keyed by the 2 sorted `v_idx` of each edge (membership only).

Faces/edges are extracted with the **same `FACE_RE2` corner table validated
byte-exact for periodic ids in `rea2nbin_light`** (Neko's `face_nodes`/`edge_nodes`
composed with the nmsh→mesh vertex swap). Zones are read straight from the `.nmsh`
to mark labeled/periodic facets; an external face whose facet is unmarked is an
unlabeled external face.

**Periodic merge.** Neko's `mesh_checker` reports the face/edge counts *after*
merging periodic point ids: on read, `apply_periodic_facet` overwrites each
periodic facet's corner ids with the record's `glb_pt_ids` before the
connectivity tables are built, so a periodic facet-pair collapses to a single
shared face and `glb_mfcs`/`glb_meds` count it once. This tool does exactly the
same — it applies the stored `glb_pt_ids` as a point-id remap before hashing — so
its `Number of faces`/`Number of edges` match `glb_mfcs`/`glb_meds` on periodic
meshes (channels, pipes, TGV, …), not just the raw per-facet totals. (The stored
`glb_pt_ids` are in `face_nodes` order, matching `FACE_RE2`; see
`mesh_reset_periodic_ids`/`hex_facet_order`.)

| | `mesh_checker` (75e6 hex, 1 rank) | `mesh_checker_light` (75e6 hex) |
|---|--:|--:|
| Dominant structures | general polymorphic point/connectivity tables + curve/zone scratch + dofmap | one face hash + one edge hash |
| Peak memory (default) | large (needs a big-memory node when run serially) | modest — a small multiple of the mesh's point data |

The `--jacobian` and `--write-zone-indices` options stay cheap because both are
streamed element-by-element (O(1) extra memory) rather than materialising a
`coef`/`dofmap`.

## Validation

There is no bundled reference output for `mesh_checker` (it needs the full Neko
build to run), so `validate_mc.py` validates several independent ways, across every
3D `.nmsh` in the tree (43 meshes):

1. **Independent recomputation.** It recomputes *every* reported quantity — points,
   faces, edges, periodic count, per-label counts, unlabeled-external count — in
   Python straight from the raw nmsh bytes, **applying the same periodic
   `glb_pt_ids` merge Neko applies on read**, and checks the tool matches **field
   for field**. Result: **43 / 43, exact.** (For the periodic meshes the merged
   `Number of faces`/`Number of edges` equal Neko's `glb_mfcs`/`glb_meds` — e.g.
   `cylinder` 5712 / 5952, `turb_channel` 17820 / 18144 — not the larger raw
   per-facet totals.)
2. **Topological identity (ground truth, independent of the tool).** Any correct
   *merged* face dedup must satisfy `2·(#faces) − (#external) = 6·nelv` (each
   internal face shared by 2 elements, each external by 1); and a fully-periodic
   mesh is a closed 3-torus, so `#faces = #edges = 3·nelv` exactly. Both hold on
   every non-degenerate mesh (e.g. `tgv/512`: 1536 faces = 1536 edges = 3·512;
   `turb_pipe`: `2·111720 − 4560 = 218880 = 6·36480`). Only the correct periodic
   merge — right `glb_pt_ids` ordering — satisfies these; the raw unmerged counts
   do not. (A single-element-thick periodic slab such as `euler_1d_sod` collapses
   each cross-section to a point; Neko produces those degenerate counts too, so
   the identity is skipped there while the field-for-field check still holds.)
3. **The diagnostic actually fires.** Stripping turb_pipe's labeled wall zone
   (4560 faces) makes the tool report **exactly 4560** unlabeled external faces and
   exit 1 — it detects precisely the newly-exposed boundary.
4. **Jacobian.** A correctly-oriented unit cube (`jac = +0.125`) is not flagged; a
   mirrored one (`jac = −0.125`) is flagged. No shipped mesh has a negative Jacobian.
5. **Zone-index field.** `zone_indices.fld` is written, parsed back, and its scalar
   field is checked to equal the mesh's per-GLL-point zone label everywhere (0
   mismatches; header/`test_pattern`/section byte-offsets all consistent).

The exit code is `0` iff there are no unlabeled external faces. Several shipped
benchmark meshes (`poisson/`, `nekbone/` boxes, some cylinders) carry **no
mesh-level BCs at all**, so their entire open surface is unlabeled external — the
tool flags them (e.g. a 128-element box → its full 160-face surface), exactly as
the real `mesh_checker` would.

> If you have a real `mesh_checker` run for one of your meshes, cross-checking its
> `Number of faces/edges` against this tool is a nice extra confirmation of the
> `glb_mfcs`/`glb_meds` semantics — the same way the real pipe `.nmsh` nailed down
> `rea2nbin_light`.

## Scope / limitations

- **3D hex meshes** (the common case). 2D meshes exit with a clear error (the full
  `mesh_checker` handles 2D via a separate path).
- **Straight-sided geometry** for `--jacobian` and `--write-zone-indices`. Both use
  the trilinear map at the GLL nodes, which is **exact for non-curved elements** (a
  trilinear field is inside the lx=3 space, so it equals `mesh_checker`'s dofmap
  geometry and Jacobian there). For **curved** (`'C'`/`'m'`) elements the curve
  deformation is *not* applied, so: the Jacobian is a valid **corner-geometry**
  inversion check but won't match `mesh_checker`'s exact GLL Jacobian; and the
  `zone_indices.fld` renders element edges as straight (the **zone labels — the
  point of the file — are exact**, only the drawn geometry is faceted). Full
  curve-aware geometry would be a natural extension.
- **Jacobian severity.** `--jacobian` flags **negative *or* zero** Jacobians and
  exits non-zero (a hard check). Neko's `mesh_checker` is slightly
  more lenient — it *warns* on strictly-negative Jacobians only (`jac < 0`) and
  does not fail. So a perfectly degenerate (exactly-zero-Jacobian) element is
  flagged here but not by `mesh_checker`; the verdict for genuinely inverted
  (`jac < 0`) elements is identical.
- The size counts match `mesh_checker`'s **serial** output; because the merged
  unique faces/edges/points are partition-independent, they also equal its
  parallel **global** counts (`glb_mfcs`/`glb_meds`/`glb_mpts`).

## Files

- `mesh_checker_light.f90` — the tool (hash + geometry modules + program).
- `validate_mc.py` — the independent cross-check, broken-mesh, Jacobian and
  zone-index-`.fld` harness (`make check`).
- `Makefile` — `make`, `make check`, `make clean`.
