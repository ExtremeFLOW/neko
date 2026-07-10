```text
     _  __  ____  __ __  ____
    / |/ / / __/ / //_/ / __ \
   /    / / _/  / ,<   / /_/ /
  /_/|_/ /___/ /_/|_|  \____/
```

# prepart_light

A lightweight, CPU-only **mesh partitioner** for Neko `.nmsh` files, in pure
scientific Python — it reads an `.nmsh`, partitions it into `nparts` parts, and
writes a reordered `.nmsh` whose elements are grouped into contiguous
per-partition blocks, exactly like `contrib/prepart` (a later *linear* read
reproduces the partition). Curves and zones are first-class. No MPI, no ParMETIS
build, no Fortran compiler — just Python.

```sh
pip install numpy scipy                        # required
pip install pymetis                            # optional: --backend metis
pip install pyamg                              # optional: faster spectral on large meshes

./prepart_light.py mesh.nmsh 8                 # -> mesh_8.nmsh (spectral, default)
./prepart_light.py mesh.nmsh 64 --metis        # multilevel k-way (best cut)
./prepart_light.py mesh.nmsh 1024 --geometric --no-stats   # near-instant fast path
make check                                     # validation across a Neko tree
```

Usage mirrors `prepart`:

```
prepart_light.py mesh.nmsh nparts [out.nmsh] [--backend spectral|metis|geometric] [--no-stats]
```

`nparts` is any positive integer; the default output is `<base>_<nparts>.nmsh`.

> `make check` reads meshes from a Neko checkout, so it expects this folder
> inside a Neko source tree (`../../examples`, `../../tests`); otherwise:
> `python3 validate_prepart.py /path/to/neko [nparts]`. Partitioning your own
> mesh needs no Neko tree.

## Why it exists — and why Python

Neko partitions meshes with **ParMETIS** (via `prepart` and the runtime load
balancer) — excellent, but an external MPI library you must build and link.
`prepart_light` is the dependency-light offline alternative: on most systems
`python + numpy/scipy` is *more* readily available than a Fortran compiler or a
ParMETIS build, and this is preprocessing — run once, cache the result.

The performance case for Python here is real, not hand-waving: the element dual
graph is built as a single **compiled sparse matrix product** `A = E·Eᵀ`
(`E` = element×vertex incidence; `A[i,j]` = number of shared vertices — exactly
the weighted dual graph Nek5000's `genmap` uses), and the partitioners are
compiled library code (ARPACK/LOBPCG, METIS). Measured on a 1M-element box into
16 parts, against the Fortran spectral tool from this same tool family:

| backend | time (1M box, 16 parts) | weighted edge cut |
|---|--:|--:|
| Fortran spectral (genmap_light) | 76 s | 3.95 % |
| `--backend metis` | **16.5 s** | **2.96 %** |
| `--backend geometric` | **4.5 s** (2.6 s with `--no-stats`) | 2.86 % |
| `--backend spectral` (default) | comparable to Fortran | ≈ Fortran |

Extrapolated to a 75M-element mesh: **minutes** with `--geometric`
(sort-dominated, no graph), **~20 min** with `--metis` (the ~2×10⁹-edge dual
graph costs ~25 GB — a single fat node, feasible). The dual-graph build itself
is ~1.3 s per million elements.

> **Tested, but please check your own case.** All three backends are validated
> across the meshes shipped with Neko (see [Validation](#validation)), but it
> remains your responsibility to confirm the result is correct for *your* mesh —
> `make check`, or run the validator on a case you trust, first.

## The three backends

All partition the **same weighted element dual graph** (edge weight = number of
shared vertices; **periodic facets are merged** through the stored `glb_pt_ids`
so periodic neighbours are graph-adjacent), and are scored by the same weighted
edge cut:

- **`spectral` (default; numpy/scipy only).** Recursive spectral bisection —
  Nek5000 `genmap`'s algorithm with a robust eigensolver. Per bisection it takes
  the few lowest non-trivial eigenvectors of the sub-graph Laplacian (dense for
  tiny subsets; shift-invert Lanczos for medium; pyamg-preconditioned LOBPCG for
  large subsets when `pyamg` is installed), forms the balanced split for each
  mode, and keeps the lowest-cut split that leaves both halves connected (BFS
  region-growing fallback otherwise). Parts differ by ≤ 1 element. Above ~50 k
  elements per subset it uses pyamg-preconditioned LOBPCG (with `pyamg`,
  `tgv/262144` partitions in ~22 s); without `pyamg` it falls back to
  shift-invert with a printed warning — fine to ~10⁵ elements, slow beyond
  (use `--metis` or `--geometric` there).
- **`metis` (needs `pip install pymetis`).** Multilevel k-way METIS on the same
  weighted graph, so it minimises the same objective. Usually the best cut and
  much faster than spectral. Balance follows METIS's semantics: the **largest**
  part stays within METIS's ~3 % tolerance of the ideal size (plus ±1 element of
  integer granularity on tiny parts) — e.g. 25 elements into 8 parts gives
  `[3,3,3,3,3,3,3,4]`. Spectral/geometric instead guarantee max−min ≤ 1. When
  `nparts` approaches `nelv` (tiny mean part size) k-way METIS can collapse
  parts; the tool detects this, retries with METIS recursive bisection, and
  errors out loudly rather than ever writing fewer parts than requested.
- **`geometric` (numpy only, no graph).** Recursive coordinate bisection on
  element centroids: split along the longest bounding-box axis at the balanced
  point, recurse. O(n log n), near-instant at any size; with `--no-stats` the
  dual graph is never built at all. **Caveat:** centroids know nothing of
  periodic wrap-around, so on strongly periodic meshes the cut can be worse than
  the graph-based backends (in practice it is often surprisingly good — 5.10 %
  cut on the periodic turbulent pipe, on par with the graph-based backends
  there).

Determinism: all three are deterministic — fixed eigensolver starting vectors
plus an eigenvector sign convention (spectral), deterministic METIS defaults,
and stable sorts (geometric). Same input + nparts ⇒ byte-identical output.

## What the reorder does (the prepart contract)

Elements are renumbered into contiguous per-partition blocks (`el_idx` becomes
`1..nelv` in file order); point ids and coordinates are **never changed**. Every
**curve record** is carried through with its element reference renumbered and
the geometric payload (`curve_data(5,12)`, `curve_type(12)`) copied verbatim.
**Periodic (type-5) zones** keep `f`/`p_f`/`glb_pt_ids` verbatim (those
reference point ids, which don't change) with `e`/`p_e` renumbered; **labeled
(type-7) zones** keep their label (`p_f`), with the `p_e`/`glb_pt_ids` fields
Neko leaves uninitialised written as deterministic zeros. Zone records of
**other (legacy) types** — e.g. the type-1/2 zones in older meshes such as the
shipped poisson/nekbone boxes, which Neko's current reader ignores — are carried
through verbatim with only the element id renumbered (with a note printed).
A linear read on `nparts` ranks then lands each part on its rank.

The reader validates section framing (truncated zone/curve sections, implausible
zone types, out-of-range curve element ids all abort with a clear error before
anything is written) and tolerates the trailing bytes some real meshes carry
past the curve section (an MPI-IO no-truncate artifact; Neko ignores them too).

I/O is numpy-vectorized (packed structured dtypes matching the 228/36/532-byte
nmsh records; one `fromfile`/`tofile` per section), so reading + writing a
1M-element mesh takes ~2 s.

## Validation

`validate_prepart.py` runs every (mesh × backend) combination across a Neko tree
and checks, independently of the tool's internals:

1. **Pure permutation** — the output's element payloads are exactly a
   permutation of the input's; `el_idx` contiguous `1..nelv`.
2. **Curves carried and bound** — every curve payload survives byte-verbatim
   *and* is attached to the geometrically identical element (matched by the
   element's vertex ids + coordinates), not merely preserved as a multiset.
3. **Zones carried and bound** — periodic zones keep `f`/`p_f`/`glb_pt_ids` on
   the geometrically identical element; labeled zones keep their label/binding.
4. **Contiguous blocks** — the file order equals the stable
   partition-grouped order (the prepart contract).
5. **Balance** — parts within ≤ 1 element (spectral/geometric) or METIS's
   tolerance; plus a weighted-cut report vs a naive linear split.

Additionally, every reordered output re-reads with `mesh_checker_light` giving
byte-identical element/point/face/edge/periodic/labeled-zone counts to the
input — confirming periodic fusion and topology survive the reorder.

## Scope / limitations

- **3D hex meshes** (as `prepart`). 2D (`gdim=2`) exits with a clear error.
- **Serial / offline** (like `prepart`): the whole mesh is loaded. Memory is
  dominated by the element records (~230 B/element) plus, for the graph-based
  backends, the weighted dual graph (~24 GB at 75M; `--geometric --no-stats`
  skips it entirely).
- **Requires numpy + scipy**; `pymetis`/`pyamg` optional. If you cannot have
  Python at all, the Fortran sibling (`genmap_light`, gfortran-only) covers the
  spectral tier.
- For production-scale *parallel* partitioning within a running solve, Neko's
  ParMETIS load balancer / `prepart` remains the tool of choice; this is the
  offline, dependency-light complement.

## Files

- `prepart_light.py` — the partitioner (single file; spectral / metis /
  geometric backends).
- `validate_prepart.py` — the permutation/binding/blocks/balance harness
  (`make check`).
- `Makefile` — `make check`, `make clean`.
