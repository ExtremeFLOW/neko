# Inviscid Compressible Flow Over a Cylinder

## Problem definition

This example simulates compressible inviscid flow over a cylinder, which is a fundamental benchmark case in computational fluid dynamics. The simulation involves a Mach 1.1 flow approaching a cylinder, creating:

- A detached bow shock upstream of the cylinder
- Regions of supersonic flow around the cylinder shoulders
- Complex flow patterns in the wake region (will be seen when a high-order artificial viscosity is implemented)

This classic test case is valuable for validating numerical schemes for compressible flows, particularly their ability to handle curved geometries and accurately capture shock waves.

## Physical parameters

- Free-stream Mach number: 1.1
- Free-stream density: 1.4
- Free-stream pressure: 1.0
- Ratio of specific heats (γ): 1.4
- Boundary conditions:
  - Inlet/far-field (zone 1): Fixed velocity, density, and pressure
  - Cylinder surface (zone 7): No-slip wall
  - Domain boundaries (zones 3, 4, 5, 6): Symmetry conditions

## Mesh

The mesh is copied from example `cyl_boundary_layer`.

## Run the case

```bash
makeneko euler_2d_cylinder.f90 && mpirun -np 16 ./neko euler_2d_cylinder.case
```
