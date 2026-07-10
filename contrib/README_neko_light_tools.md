```text
     _  __  ____  __ __  ____
    / |/ / / __/ / //_/ / __ \
   /    / / _/  / ,<   / /_/ /
  /_/|_/ /___/ /_/|_|  \____/
```

# Neko light tools

Five standalone, CPU-only companions to Neko's `contrib/` tools. Each folder is
self-contained (tool + validation harness + Makefile + doc) and is designed to
sit at `contrib/<tool>/` inside a Neko source checkout, so that `make check` can
find the meshes shipped with Neko (`../../examples`, `../../tests`). Using the
tools on your own meshes needs no Neko tree at all.

The common idea: Neko's own `contrib` tools build a full `mesh_t` through Neko's
general-purpose, production data structures — robust and reused everywhere, but
memory-hungry on large meshes, and CPU-only (so a GPU-configured Neko must be
rebuilt for CPU just to run them). These reimplementations use plain typed
arrays and no Neko initialisation: the same results in a small fraction of the
memory, needing only a Fortran compiler or a Python interpreter.

| Tool | What it does | Build / run needs | `make check` needs |
|---|---|---|---|
| `rea2nbin_light` | Serial `.re2` → `.nmsh` converter, byte-exact with `contrib/rea2nbin` (periodic BCs, curved elements) | gfortran | python3 |
| `genmeshbox_light` | O(1)-memory box-mesh generator, byte-exact with `contrib/genmeshbox` | gfortran | python3 |
| `mesh_checker_light` | Low-memory `.nmsh` diagnostics (sizes, BC zones, Jacobian, zone-index field) | gfortran | python3 |
| `genmap_light` | Spectral-bisection mesh partitioner (Nek5000 `genmap` algorithm), reorders like `prepart` | gfortran | python3 + numpy + scipy |
| `prepart_light` | Mesh partitioner in scientific Python: spectral / METIS / geometric backends | python3 + numpy + scipy (optional: pymetis, pyamg) | python3 + numpy + scipy |

Two partitioners are included on purpose: `genmap_light` when a Fortran compiler
is all you have (deterministic, dependency-free), `prepart_light` when Python
with NumPy/SciPy is available (faster; best cut with `pip install pymetis`; a
near-instant geometric mode for very large meshes). Both write a reordered
`.nmsh` whose linear read reproduces the partition — the same contract as
Neko's `contrib/prepart`.

Every tool prints a `--help`/usage message, validates its inputs with clean
one-line errors, and is deterministic (same input, same output bytes). See each
folder's `.md` for the algorithm, the validation evidence, and the scope and
limitations. All tools are validated against the meshes shipped with Neko, but
it remains your responsibility to confirm the result is correct for your own
mesh — each folder's `make check` (or validator script) makes that easy to do
on a case you trust.
