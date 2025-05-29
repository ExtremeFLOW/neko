# 2D Euler Smooth Wave Propagation

Test case for compressible Euler equations demonstrating advection of a smooth density wave with periodic boundary conditions.

## Problem Description

A smooth density wave advects in a periodic domain with constant background pressure:

Initial conditions:

- `ρ(x,y,0) = 1.0 + 0.2 sin(2π(x + y))`
- `u = 2.5` (constant)
- `v = -0.5` (constant)
- `w = 0.0`
- `p = 1.0` (constant)

Domain: `[0,1] × [0,1]` with periodic boundaries.

## Running the case

1. Generate 2D periodic mesh (100 x 100 x 1 with grid size 0.01):

```bash
genmeshbox 0 1 0 1 0 0.01 100 100 1 .true. .true. .true.
```

2. Run the simulation:

```bash
/path/to/neko_install/bin/makeneko euler_2d_smooth.f90
mpirun -n 4 neko euler_2d_smooth.case
```

## Reference

Nazarov, M., & Larcher, A. (2017). Numerical investigation of a viscous regularization of the Euler equations by entropy viscosity. Computer Methods in Applied Mechanics and Engineering, 317, 128-152.
