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

Generate the mesh with `gmsh` and convert it with `gmsh2nmsh`, making the front
and back surfaces a periodic pair:

```bash
gmsh step.geo -3 -o step.msh
gmsh2nmsh step.msh --periodic=front:back
```

The inlet, outlet, top wall and bottom wall become labeled zones 1 to 4.

## Run the case

```bash
makeneko step.f90 && mpirun -np 16 ./neko step.case
```
