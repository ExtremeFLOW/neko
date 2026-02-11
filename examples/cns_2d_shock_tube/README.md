# Reflected Shock-Boundary Layer Interaction in a Shock Tube

## Problem Description

A 2D compressible Navier-Stokes simulation of a shock tube with viscous walls.
A strong incident shock (density ratio 100:1) travels rightward, reflects off the
right wall, and interacts with the developing boundary layers on the top and
bottom walls.

### Initial Conditions

Diaphragm at `x = 0.5`, fluid at rest (`u = v = 0`):

| Region | ρ | p |
|--------|-------|------------|
| Left (`x < 0.5`) | 120.0 | 120.0 / γ |
| Right (`x > 0.5`) | 1.2 | 1.2 / γ |

### Boundary Conditions

- **Left/Right walls** (`x = 0, 1`): Symmetry (adiabatic slip walls)
- **Top/Bottom walls** (`y = 0, 1`): No-slip (adiabatic viscous walls)

### Parameters

- gamma = 1.4 (ideal gas)
- Re = 200 (set in `2d_shock_tube.f90`, change to 1000 for more complex flow)
- End time: t = 0.6
- Incident shock reflects off right wall at approximately t = 0.2

## Running the Case

1. Generate a 2D mesh (50 x 50 elements):

```bash
genmeshbox 0 1 0 1 0 0.02 50 50 1 .false. .false. .true.
```

2. Compile and run:

```bash
makeneko 2d_shock_tube.f90
mpirun -np 4 ./neko 2d_shock_tube.case
```

## Expected Results

1. Right-moving shock, contact discontinuity, and left-moving rarefaction fan
2. Shock reflection off the right wall (t = 0.2)
3. Reflected shock interaction with boundary layers on top/bottom walls
4. Boundary layer separation and complex vortical structures near walls
