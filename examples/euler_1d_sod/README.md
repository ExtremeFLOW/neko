# Sod's Shock Tube Problem

This is a classical test case for compressible fluid dynamics that demonstrates the ability of numerical schemes to handle discontinuities and shock waves. The problem was first introduced by Gary Sod in 1978 and has become a standard benchmark for computational fluid dynamics codes.

## Problem Description

The problem consists of a tube divided into two regions by a membrane. At `t=0`, the membrane is removed, creating a discontinuity in density and pressure:

Left state (`x < 0.5`):

- `ρ_L = 1.0`
- `u_L = 0.0`
- `p_L = 1.0`

Right state (`x > 0.5`):

- `ρ_R = 0.125`
- `u_R = 0.0`
- `p_R = 0.1`

The simulation runs until `t = 0.2`.

## Case Setup

The case uses the compressible Euler equations with:

- `γ = 1.4` (ratio of specific heats)
- 1D domain `[0,1]`
- Zero velocity initial condition
- Transmissive boundary conditions

The initial conditions are set in the user file `sod.f90`.

## Running the case

1. Generate the mesh (100 x 1 x 1 with grid size 0.01) using:

```bash
genmeshbox 0 1 0 0.01 0 0.01 100 1 1 .false. .true. .true.
```

2. Run the simulation:

```bash
/path/to/neko_install/bin/makeneko sod.f90
mpirun -n 4 neko sod.case
```

## Expected Results

The solution should show:

1. A right-moving shock wave
2. A right-moving contact discontinuity
3. A left-moving rarefaction wave

## Reference

Sod, Gary A. (1978). "A Survey of Several Finite Difference Methods for Systems of Nonlinear Hyperbolic Conservation Laws". Journal of Computational Physics. 27: 1-31
