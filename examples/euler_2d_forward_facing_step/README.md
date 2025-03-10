# Wind tunnel with a forward-facing step

## Problem definition

This example simulates compressible flow over a forward-facing step in a wind tunnel, which is a classical test case in computational fluid dynamics. The case involves a Mach 3 flow entering a wind tunnel that has a step geometry, creating a complex shock wave system including:

- An attached oblique shock at the step
- Reflected shocks from the top wall
- Possible recirculation zones near the step

This benchmark problem is useful for validating numerical schemes for compressible flows, particularly their ability to capture shock waves and complex flow features.

## Physical parameters

- Inlet Mach number: 3.0
- Inlet density: 1.4
- Inlet pressure: 1.0
- Ratio of specific heats (γ): 1.4
  - Inlet (zone 1): Fixed velocity, density, and pressure
  - Top/bottom walls (zones 3, 4): Symmetry conditions
  - Front/back (zones 5, 6): Periodic pair

## References

1. Woodward, P. R., & Colella, P. (1984). The numerical simulation of two-dimensional fluid flow with strong shocks. Journal of Computational Physics, 54(1), 115-173.

2. Nazarov, M., & Larcher, A. (2017). Numerical investigation of a viscous regularization of the Euler equations by entropy viscosity. Computer Methods in Applied Mechanics and Engineering, 317, 128-152.

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

## Run the case

```bash
makeneko step.f90 && mpirun -np 16 ./neko step.case
```
