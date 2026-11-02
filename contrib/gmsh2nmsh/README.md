## `gmsh2nmsh`

`gmsh2nmsh` converts a Gmsh `.msh` mesh directly to Neko's `.nmsh` format. It
replaces the route through Nek5000's `.re2` format with Nek5000's `gmsh2nek`
and `rea2nbin`, and needs no interactive input.

## Usage

```text
gmsh2nmsh mesh.msh [mesh.nmsh] [--periodic=PAIRS] [--tol=VALUE]
          [--curve-tol=VALUE] [--linear]
```

- `--periodic=PAIRS` makes pairs of physical groups periodic. Groups are given
  by tag or name, for example `--periodic=1:2` or
  `--periodic=inlet:outlet,front:back`. The option can be repeated.
- `--tol=VALUE` sets the absolute periodic matching tolerance. The default is
  `1e-6` times the shortest facet edge of the pair, but at least `1e-12` times
  the largest coordinate magnitude of the pair.
- `--curve-tol=VALUE` sets how far an edge midpoint, relative to the edge
  length, must be from the straight edge to be stored as curved. The default is
  `1e-4`.
- `--linear` ignores high-order nodes and writes straight-sided elements.

## Supported input

- MSH format 4.1 or 2.2, ASCII or binary. Older and partitioned files are
  rejected.
- 3D meshes of 8, 20 or 27 node hexahedra, with quadrilateral boundary facets.
- 2D meshes of 4, 8 or 9 node quadrilaterals in a plane of constant z, with
  line boundary facets. Neko extrudes these into one layer of hexahedra.

## Behaviour

- Boundary facets in physical groups (surfaces in 3D, curves in 2D) become
  labeled zones. If all their tags are in `[1, 20]`, the zone index is the
  physical tag, as with `gmsh2nek` + `rea2nbin`. Otherwise the groups are
  numbered `1, 2, ...` in order of their tags. The mapping is printed.
- Periodic groups must be related by a translation, which is taken from the
  centroids of the two groups. Every facet must have exactly one partner.
- Curved edges of second order elements are stored as midpoint curves, the
  same representation `gmsh2nek` produces. Face and volume centre nodes are
  not used.
- Left-handed elements are mirrored, so all elements are right-handed.

## Building

The tool is built and installed with the other contrib tools. It has its own
configure script, which Neko's configure runs with the same options, and it is
compiled with the plain Fortran compiler `FC` rather than `MPIFC`. It therefore
links neither MPI nor any of Neko's dependencies. It can also be compiled on
its own from the single source file, for example on the machine where the mesh
is generated:

```text
gfortran -O2 gmsh2nmsh.F90 -o gmsh2nmsh
```

Built this way, the banner printed at start shows the version as `unknown`.

## Recommended workflow

1. Mesh with hexahedra (3D) or quadrilaterals (2D), and put the volume
   (or surface) and every boundary in physical groups. When physical groups
   exist, Gmsh only saves elements that belong to one.
2. Use `-order 2` if the boundary is curved, for example
   `gmsh pipe.geo -3 -order 2`.
3. Run `gmsh2nmsh pipe.msh --periodic=inlet:outlet`.
4. Run `mesh_checker pipe.nmsh` on the result.
