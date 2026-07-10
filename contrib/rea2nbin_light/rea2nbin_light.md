```text
     _  __  ____  __ __  ____
    / |/ / / __/ / //_/ / __ \
   /    / / _/  / ,<   / /_/ /
  /_/|_/ /___/ /_/|_|  \____/
```

# rea2nbin_light

A standalone, low-memory `re2 → nmsh` mesh converter — a lightweight, CPU-only
alternative to Neko's `contrib/rea2nbin` for the common case of a **serial**
conversion.

Needs only a Fortran compiler (no MPI, no Neko library):

```sh
make                                   # gfortran -O2 -o rea2nbin_light rea2nbin_light.f90
./rea2nbin_light mesh.re2 [mesh.nmsh]
make check                             # byte-exactness vs examples/hemi (golden pair)
make check-roundtrip                   # geometry+topology check vs every 3D mesh in Neko
make check-periodic                    # periodic-BC identification vs every periodic mesh
make check-curved                      # curved-element ('C'/'m') records vs every curved mesh
make check-all                         # the four in-tree suites above
make check-golden RE2=p.re2 NMSH=p.nmsh # TRUE oracle: byte-diff vs a REAL Neko .nmsh
```

> The `make check*` targets read the meshes shipped with Neko, so they expect
> this folder to live inside a Neko source checkout (e.g. dropped into
> `contrib/`, so `../../examples` resolves). To run them from elsewhere, point
> the validator at your Neko tree directly, e.g.
> `python3 validate_roundtrip.py ./rea2nbin_light /path/to/neko`. Converting your
> own mesh (`./rea2nbin_light mesh.re2`) needs no Neko tree at all.

## Why it exists

Neko's `contrib/rea2nbin` reads the `.re2` into a full `mesh_t` and writes it
out. Building that `mesh_t` goes through Neko's general-purpose, production hash
tables — the same robust, well-tested structures reused all over the code base.
They're sized generously and store polymorphic keys, which is exactly what you
want for a solver, but it means a fair amount of memory per element when all you
actually want is to convert a mesh once.

`rea2nbin_light` is a lightweight reimplementation of just that conversion. It
uses plain typed arrays instead of the general hash tables and does **no Neko
initialisation**, so it produces the same `.nmsh` in a small fraction of the
memory — comfortably converting meshes that would otherwise need a very large
node.

A couple of nice side effects of not linking Neko:

- **No GPU rebuild needed.** The `contrib` tools are CPU-only. If your Neko is
  configured for GPUs, you'd normally have to rebuild it for CPU just to run
  them. `rea2nbin_light` needs only a Fortran compiler, so that whole step goes
  away.
- **Runs anywhere** — a login node, a laptop, a small VM — since there's no MPI
  and no Neko library to carry along.

> **Tested, but please check your own case.** This tool has been validated
> against the meshes shipped with Neko and against a real, wall-resolved pipe
> (see [Validation](#validation)), and it is byte-exact with the real
> `rea2nbin` on those. It is nonetheless your responsibility to make sure the
> result is correct for *your* mesh — run one of the checks below on a case you
> trust before relying on it in production.

## What `rea2nbin_light` does instead

- A single **typed** open-addressing coordinate hash (structure-of-arrays:
  `real(dp) kx/ky/kz`, `integer id`), **no polymorphism, no per-entry heap**,
  sized to the number of *unique* points (~`nelv` for a hex mesh, not `8·nelv`)
  and **grown on demand**. ~28 bytes per slot.
- Elements are **streamed** straight from the re2 record to the nmsh record —
  the mesh is never materialised.
- **Periodic (`P`) BCs** are reproduced by a **line-for-line replica of Neko's
  `mesh_create_periodic_ids`**: the same **fixed 3 sweeps** over the `P` facets
  (Neko's `do j = 1, 3` in `re2_file_read_bcs`), each setting `id = min(id_i,
  id_j)` in place on the coordinate-matched corner pairs, with Neko's
  `match /= 1` abort guard. Using Neko's exact 3-sweep (rather than a full
  union-find) makes the merged `glb_pt_ids` byte-identical for *any* mesh —
  including pathological >3-hop merge chains — not just the ≤3-direction case.
  This handles 1, 2 and 3 simultaneous periodic directions (channel, TGV/HIT, …);
  see validation §3 and the real-pipe oracle §5.
- **Curved elements** (`'C'` circle / `'m'` midside re2 curve records) are
  **passed through** to nmsh curve records exactly as Neko does — essential for
  curved geometry such as a pipe/cylinder wall, which the solver bends onto the
  arcs in `dofmap_generate_xyz`. Records are stored compactly (not the dense
  `5·12·nelv` array Neko allocates) and re-emitted one per curved element in
  ascending element order via an O(n) counting sort; a single unsupported type
  (`'s'`/`'e'`/other) makes Neko treat the mesh as non-curved, and the tool
  matches that (writes `ncurves = 0`). See validation §4.

Measured on `examples/hemi`: **2718 unique points** vs `8·2042 = 16336` vertex
slots (true unique count ~1.33·nelv), confirming the `8·nelv` sizing is ~6×
oversized before the general table rounds up.

| | `rea2nbin` (75e6 hex) | `rea2nbin_light` (75e6 hex) |
|---|--:|--:|
| Dominant structures | general polymorphic point/connectivity tables + per-element stacks + dense `5·12·nelv` curve buffer | one `~2·uniq` typed coordinate hash + compact per-record curve store |
| Peak memory | very large (needs a big-memory node) | **~10–15 GB** straight-sided; **+~5 GB** if fully curved |

(The curve store is compact — one entry per re2 curve record (`~52 B`), not
Neko's dense `5·12·nelv` (`528 B·nelv`, ~40 GB at 75e6) that exists whether or
not an element is curved. So even a fully-curved 75e6 pipe stays well under
~20 GB.)

## Validation

Reproduce everything with `make check-all`.

### 1. Byte-exactness vs the golden pair (`examples/hemi/hemi.re2` → `hemi.nmsh`)

`hemi` is the only `.re2` shipped with Neko, so it is the one true end-to-end
reference. Converting it and comparing to the bundled `hemi.nmsh`:

- **Element section (header + hex records, 465 584 B): byte-for-byte identical.**
  The coordinate dedup, the first-appearance `v_idx` numbering, the (net-identity)
  vertex ordering, and the packed binary layout all match Neko exactly.
- **`nzones` (1232), `ncurves` (0), total file size, and every meaningful zone
  field (`e`, `f`, `p_f`, `type`) match** position-by-position for all zones.
- The **only** differing bytes (8587 / 509 944 = 1.68 %) are the `p_e` and
  `glb_pt_ids` fields of labeled zones, which Neko's own writer leaves
  uninitialised (`nmsh_file.f90` sets only `%e/%f/%p_f/%type` for `type = 7`
  zones, so it writes stale heap into the file). Those bytes are ignored on read
  and are not reproducible by construction; `rea2nbin_light` writes deterministic
  zeros there, producing a *cleaner* but functionally identical file.

### 2. Geometry + topology across every 3D mesh in the Neko tree (43 meshes)

Since there is only one `.re2`, `validate_roundtrip.py` round-trips each shipped
3D `.nmsh` instead: `nmsh → synthesize a .re2 → rea2nbin_light → nmsh'`, and
checks the two invariants that define a correct mesh import. A bad import can only
show up as **wrong coordinates** (geometry corrupted) or **wrong sharing** (two
coincident vertices not merged, or distinct ones merged = topology corrupted).

Corpus: 43 meshes, **25 → 262 144 elements** — cylinders, boxes, turbulent pipe,
hemisphere, TGV, Poisson/nekbone stacks, lid, Rayleigh–Bénard, forward-facing
step, and more.

| Invariant | Result |
|---|---|
| **Coordinates bit-exact** (geometry) | **43 / 43** |
| **Point dedup partition identical** (topology / sharing) | **43 / 43** |
| `v_idx` *values* exact | 34 / 43 |
| remaining 9 (`genmeshbox` boxes): id map is a clean **bijection** | verified |

The 9 meshes whose `v_idx` *values* differ are all `genmeshbox`-generated boxes.
Their `.nmsh` numbers points in lexicographic order because they were never made
from a `.re2`; `rea2nbin_light` (like Neko's re2 reader) numbers points in
first-appearance order. The two id sets are in exact one-to-one correspondence
(a benign **relabeling**, not a topology change) — so the mesh is identical. Run
the original `rea2nbin` on a box's `.re2` and it would produce the same
first-appearance numbering, which is exactly why the tool is byte-exact against
the genuine re2-derived reference (`hemi`).

### 3. Periodic-BC identification across every periodic mesh in the tree (19 meshes)

Periodicity is the one piece of information a `.nmsh` carries that isn't in the
element records: Neko merges the point-ids of matched periodic facets
(`mesh_create_periodic_ids`) and writes the *merged* ids into each periodic
(type-5) zone's `glb_pt_ids`. `rea2nbin_light` must rebuild exactly that merge
from the raw re2 `P` records. `validate_periodic.py` synthesizes the equivalent
`.re2` (`P` records reconstructed from each mesh's periodic zones), runs the
tool, and checks three things against the shipped `.nmsh`:

| Check | What it proves | Result |
|---|---|--:|
| **Zone records** (`e`, `f`, `p_e`, `p_f`) match position-by-position | facets paired correctly | **19 / 19** |
| **Merge partition** of physical coordinates identical | the *right* nodes are fused (numbering-independent) | **19 / 19** |
| **One merged id per coordinate** in both files | no split/contradictory merges | **19 / 19** |

Corpus: 19 periodic meshes covering all three multiplicities —
**1 direction** (cylinders, `turb_pipe`, `rayleigh_benard`), **2 directions**
(`TS_channel`, `turb_channel`), and **3 directions** (`tgv/512`, `tgv/32768`,
`tgv/262144` — 24 576 periodic facets). The 3-sweep in-place min-id merge
reproduces Neko's `mesh_create_periodic_ids` exactly on every one, so a channel
(two periodic pairs) or a TGV/HIT box (three periodic pairs) converts identically
to what Neko would build.

### 4. Curved-element records across every curved mesh in the tree (11 meshes)

Curve records are the geometry a straight-edged hex can't carry: Neko reads each
re2 curve record (`elem`, `edge`, 5 values, type char) into
`curve_data(5,12,elem)`/`curve_type(12,elem)`, writes one nmsh curve record per
curved element, and later **bends the element faces onto them** in
`dofmap_generate_xyz` — circle arcs for `'C'`, explicit midside nodes for `'m'`.
Dropping them would flatten a pipe/cylinder wall into facets. `rea2nbin_light`
passes them through. There is no shipped `.re2` with *supported* curves (`hemi`'s
are `'s'`, which Neko itself drops), so `validate_curved.py` reconstructs the
equivalent re2 curve records from each curved `.nmsh` and compares the tool's
nmsh curve section to the reference, record by record:

| Check | Result |
|---|--:|
| **curve record count** matches | **11 / 11** |
| **element ids** match (ascending order) | **11 / 11** |
| **`curve_data(5,12)` bit-exact** (all 60 doubles/record) | **11 / 11** |
| **`type(12)` exact** | **11 / 11** |

Corpus: 11 curved meshes — the cylinders, `cyl_boundary_layer`,
`rayleigh_benard_cylinder` (100 % curved), and **`turb_pipe`** (36 480 elements,
17 860 curved — the pipe case directly analogous to a wall-resolved DNS pipe).
A separate synthetic test confirms the two paths the shipped meshes don't
exercise: mixed `'C'`+`'m'` records on one element round-trip per-edge, and an
`'s'` record correctly forces `ncurves = 0` (Neko-identical drop-all).

### 5. True oracle — head-to-head vs a real Neko `.nmsh` (a genuine pipe)

The §1–§4 harnesses synthesize their `.re2` from a `.nmsh`, so they could in
principle share a blind spot with the tool. `validate_golden.py` removes that: it
byte-diffs the tool's output against a **`.nmsh` produced by the real Neko
`rea2nbin`**, section by section. Run on a genuine **Re_τ=550 pipe**
(`P550…Nc2198_Nz5`, 10 990 elements = 2198 × 5; streamwise-periodic wall-resolved,
all features at once):

| Section | Differing bytes |
|---|--:|
| header | **0** |
| elements (2.5 MB: dedup ids, vertex order, coords) | **0** |
| **periodic (type-5) zones — incl. merged `glb_pt_ids`** | **0** |
| **curves** (4.6 MB, 8640 curved elements) | **0** |
| labeled (type-7) zone `p_e`/`glb_pt_ids` | 2621 (Neko-uninitialised heap) |

**Total file size: identical. The only differing bytes are the type-7 labeled-zone
`p_e`/`glb_pt_ids` fields that Neko's writer never initialises** (it leaves stale
heap — values like `134239120` — while the tool writes deterministic `0`); the
nmsh reader never reads those fields, so the meshes are functionally identical and
the tool's output is the cleaner of the two. Point `make check-golden` at any real
pair you have (`RE2=…  NMSH=…`).

The converter has additionally been cross-checked by section-by-section source
comparison against Neko's re2 reader and nmsh writer, plus deliberate
break-attempt tests on malformed inputs. Verdict: **byte-exact on well-formed 3D
`#v002` meshes, modulo the type-7 heap fields.**

### Robustness

Where the real `rea2nbin` calls `neko_error` and writes no file, the tool does
the same rather than crash or silently emit a wrong mesh:

- **Bounds guards** — curve/BC element `∉ [1,nelv]`, curve edge `∉ [1,12]`, BC
  face `∉ [1,6]`, boundary label `∉ [1,20]`, or a periodic partner out of range
  give a clean error and exit 1.
- **Curve type** is decoded on the **first character only** (matching Neko's
  `character(len=1)` field), so a NUL/junk-padded type field from a non-Nek5000
  writer does not silently drop all curves.
- **`NEKO_PERIODIC_TOL`** environment variable honoured (default `1e-7`), exactly
  as Neko does.

## Scope / limitations

- **3D hex meshes** (the common case). **2D meshes are not handled** — the tool
  exits with a clear error (the real `rea2nbin` *does* convert 2D quads; if you
  need 2D, use it). Other element types are likewise out of scope.
- **Curved elements** *are* supported: `'C'` (circle) and `'m'` (midside) re2
  curve records are converted to nmsh curve records bit-exact (validation §4).
  Only the two curve types Neko itself supports are handled; `'s'`/`'e'`/unknown
  types make Neko drop all curves, and the tool matches that. (`hemi`'s 700 curve
  records are `'s'`, so both it and the tool write `ncurves = 0` — which is why
  the golden-pair byte-exactness in §1 still holds.)
- **Periodic (`P`) BCs** *are* supported (exact 3-sweep min-id replica, 1/2/3
  directions; validation §3, §5). Labeled/named BCs (`MSH`/`EXO`, `W`, `v/V`,
  `O/o`, `SYM`, `ON`, `s*`) are fully ported, including Neko's user-label offset
  and first-appearance ordering.
- **Header versions**: `#v002`/`#v003` (double precision) are validated
  end-to-end by the shipped suites. `#v001` (single precision, same layout) and
  `#v004` (EXO, i16 element counts) parsing is implemented but not exercised by
  the in-tree meshes — if you rely on either, validate a real pair once with
  `make check-golden RE2=… NMSH=…`.

For meshes within scope this is a validated, much lighter serial converter,
byte-exact with the real `rea2nbin` (modulo the type-7 fields Neko leaves
uninitialised). As always, confirm it on a case you trust before relying on it.

## Files

- `rea2nbin_light.f90` — the converter (module + program).
- `validate_roundtrip.py` — the geometry+topology harness (`make check-roundtrip`).
- `validate_periodic.py` — the periodic-BC harness (`make check-periodic`).
- `validate_curved.py` — the curved-element harness (`make check-curved`).
- `validate_golden.py` — the real-`.nmsh` byte-diff oracle (`make check-golden`).
- `Makefile` — `make`, `make check`, `make check-roundtrip`, `make check-periodic`,
  `make check-curved`, `make check-all`, `make check-golden`, `make clean`.
