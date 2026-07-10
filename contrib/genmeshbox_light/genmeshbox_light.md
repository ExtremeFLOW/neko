```text
     _  __  ____  __ __  ____
    / |/ / / __/ / //_/ / __ \
   /    / / _/  / ,<   / /_/ /
  /_/|_/ /___/ /_/|_|  \____/
```

# genmeshbox_light

A standalone, **O(1)-memory** box-mesh generator for Neko — a lightweight,
CPU-only alternative to `contrib/genmeshbox` that writes the `.nmsh` by
**streaming**, without ever building a `mesh_t`.

Needs only a Fortran compiler (no MPI, no Neko library):

```sh
make                                        # gfortran -O2 -o genmeshbox_light genmeshbox_light.f90
./genmeshbox_light 0 1 0 1 0 1 8 8 8 .true. .true. .false.
make check                                  # byte-exactness vs every genmeshbox box in Neko
```

> `make check` reads the boxes shipped with Neko, so it expects this folder to
> live inside a Neko source checkout (e.g. dropped into `contrib/`, so
> `../../examples` and `../../tests` resolve). To run it from elsewhere, point the
> validator at your Neko tree directly:
> `python3 validate_gmb.py ./genmeshbox_light /path/to/neko`. Generating a box
> (`./genmeshbox_light ...`) needs no Neko tree at all.

The command line follows `genmeshbox`, with two additions — an optional
trailing output filename and the `--direct` flag:

```
genmeshbox_light x0 x1 y0 y1 z0 z1 nelx nely nelz \
    [px py pz] [dist_x dist_y dist_z] [out.nmsh] [--direct]
```

`px/py/pz` (`.true./.false.`) make that direction periodic; `dist_*` is `uniform`
or a CSV file of the `nel+1` grid coordinates, read **exactly as Neko reads it**:
either a single line of comma/space-separated values, or a header line followed
by the values (when the file has more than one line, Neko treats the first line
as a header and skips it). A distribution file that works with `genmeshbox` works
here and yields the same grid. Default output is `box.nmsh`. `--direct` is
explained under *Coordinates* below.

## Why it exists

Neko's `contrib/genmeshbox` assembles the box as a full `mesh_t` and writes it
out. That goes through Neko's general-purpose, production hash tables — the same
robust, well-tested structures used throughout the solver. They're sized
generously and store polymorphic keys, which is the right trade-off inside Neko,
but it means a lot of memory for a large box when every point id and vertex
coordinate is really just a closed-form function of the grid indices.

`genmeshbox_light` is a lightweight reimplementation of just the box generator.
It computes each point id and coordinate directly, streams the records to disk,
and does **no Neko initialisation** — so it produces the same `.nmsh` in a tiny
fraction of the memory (a few kilobytes of working set regardless of box size).

A couple of nice side effects of not linking Neko:

- **No GPU rebuild needed.** The `contrib` tools are CPU-only. If your Neko is
  configured for GPUs, you'd normally have to rebuild it for CPU just to run
  them. `genmeshbox_light` needs only a Fortran compiler, so that step goes away.
- **Runs anywhere** — a login node, a laptop, a small VM — since there's no MPI
  and no Neko library to carry along.

> **Tested, but please check your own case.** This tool is validated
> byte-for-byte against every genmeshbox-produced box shipped with Neko and
> head-to-head against a freshly-built current `genmeshbox` (see
> [Validation](#validation)). It is still your responsibility to make sure the
> result is correct for *your* box — `make check` on a case you trust before
> relying on it.

## What `genmeshbox_light` does instead

For a box, *everything* is a closed form of the grid indices, so no mesh, no hash
tables, and no dedup are needed:

- **Point ids** are lexicographic: `id = 1 + ix + iy·(nx+1) + iz·(nx+1)(ny+1)`.
- **Element records** are streamed in genmeshbox's loop order (`ez` outer … `ex`
  inner), each with the 8 corner ids/coordinates in the standard vertex order.
- **Periodic zones** carry the *merged* point ids, computed by wrapping each
  periodic direction's high boundary to its low one (`px == nx → 0`) — the exact
  min-id result of `create_periodic_ids`, in closed form.
- **Labeled zones** get the label in `p_f`; their `p_e`/`glb_pt_ids` are written
  zero (Neko leaves those uninitialised — stale heap).

Memory is just the three 1-D coordinate arrays (`nx+ny+nz+3` doubles) plus a
streaming write buffer — a few kilobytes for a 75M box.

| | `genmeshbox` (75e6 box) | `genmeshbox_light` (75e6 box) |
|---|--:|--:|
| Dominant structures | general polymorphic mesh/connectivity tables + points | three 1-D coord arrays |
| Peak memory | large (needs a big-memory node) | **~a few KB** |

## Coordinates — `--direct` vs the default

genmeshbox builds the 1-D grid line coordinates by **cumulative summation**
(`cumm(1)=x0; cumm(k+1)=cumm(k)+el_len`, with `el_len=(x1-x0)/n`), together with
the distribution-file support added in **Sep 2024**. This *drifts*: for
`x0=0, x1=4, n=18`, grid line 6 is `1.3333333333333335` and the last line is
`4.000000000000001` (not exactly `4.0`). `genmeshbox_light`'s **default reproduces
that cumulative-sum arithmetic exactly** (identical operations, order and types →
identical bits) — **verified byte-for-byte against the real, freshly-built current
genmeshbox** (see Validation).

Every box **shipped in the Neko examples predates that change** and instead stores
`x0 + ((x1-x0)/n)·k` (element length formed once, multiplied by `k`, landing on
exact values like `2.0`). Pass `--direct` to reproduce those; with a distribution
file, `--direct` uses the file's values verbatim as the grid lines.

## Validation

`validate_gmb.py` regenerates every genmeshbox-produced box in the Neko tree
(detected by the lexicographic point-id scheme — benchmark boxes from other
generators, which number points sequentially per element, are skipped) and
byte-compares. It feeds the exact grid lines as `--direct` distribution files, so
the check is independent of the coordinate-formula detail and exercises the part
that actually matters — element ordering, point ids, and the zone / periodic-merge
logic.

Result: **9 / 9 byte-exact** — header, every element record, and every periodic
(type-5) zone identical; labeled (type-7) zones differ **only** in the
`p_e`/`glb_pt_ids` fields Neko leaves uninitialised (the tool writes zeros).

Corpus covers **1, 2 and 3 periodic directions**: `lid`/`lid2d`/`recycling`
(z), `rayleigh_benard` (xy), `turb_channel`/`TS_channel` (xz),
`euler_1d_sod` (yz), `advecting_cone`/`euler_2d_smooth` (xyz) — up to
`euler_2d_smooth` at 10 000 elements.

**Head-to-head against the real tool (definitive).** The current `genmeshbox` was
built from source and run on three cases; `genmeshbox_light` (default arithmetic)
matched its output **byte-for-byte** — header, all element records, all periodic
zones, and total file size — with the only differences being the uninitialised
labeled-zone fields:

| case | args | result |
|---|---|---|
| 18³, 1-periodic | `0 4 -1 1 0 1.5 18 18 18 .false. .true. .false.` | **byte-exact** (the accumulate-drift case) |
| 8³, 3-periodic | `0 1 0 1 0 1 8 8 8 .true. .true. .true.` | **byte-exact** |
| 5×5×1, 2-periodic | `0 4.5 0 4.5 0 1 5 5 1 .true. .true. .false.` | **byte-exact** |

The 18³ case is the one where accumulate (`1.333…35`) and the closed form
(`1.333…33`) diverge — the default matches the accumulating original exactly.

## Scope / limitations

- **3D hex boxes** (genmeshbox's purpose). Non-uniform grids via a distribution
  file are supported (uniform or per-direction).
- Byte-identical to the **current** genmeshbox by default (cumulative-sum
  coordinates); `--direct` reproduces pre-2024 / shipped-example boxes.
- The only bytes that can differ from genmeshbox are the labeled-zone
  `p_e`/`glb_pt_ids` that Neko writes from uninitialised memory.

## Files

- `genmeshbox_light.f90` — the generator (single program).
- `validate_gmb.py` — the byte-exactness harness (`make check`).
- `Makefile` — `make`, `make check`, `make clean`.
