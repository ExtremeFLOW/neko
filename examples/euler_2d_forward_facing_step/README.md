# Wind tunnel with a forward-facing step

## Mesh generation

Read `step.geo` in `gmsh`

```bash
gmsh step.geo
```

Generate a mesh and export it to `step.msh` in a format compatible with `gmsh2nek`, see [this](https://nek5000.github.io/NekDoc/tools/gmsh2nek.html).

Run `gmsh2nek` which will output `step.re2`.

```bash
 ******************************************************
 Fluid mesh boundary info summary
 BoundaryName     BoundaryID
 inlet           1
 outlet           2
 top_wall           3
 bottom_wall           4
 front           5
 back           6
 ******************************************************
```

Make (5, 6) a period pair.

Run `rea2nbin step.re2` which will output `step.nmsh`.
