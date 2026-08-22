# Mach 800 Jet

This case exercises the Euler IDP method on a very high Mach-number jet. A
light ambient gas initially fills the domain. A narrow jet enters through the
bottom boundary and drives a bow shock through the ambient gas.

## Setup

- Computational half-domain: `[0, 0.5] x [0, 1.5] x [0, 0.01]`
- Mesh: `20 x 60 x 1` affine elements, periodic in the thin direction
- Symmetry boundary at `x = 0`
- Ambient state: `(rho, u, v, p) = (0.14, 0, 0, 1)`
- Jet state on `y = 0`, `x <= 0.05`:
  `(rho, u, v, p) = (1.4, 0, 800, 1)`
- Ratio of specific heats: `gamma = 1.4`
- Polynomial degree: 3
- Time integrator: SSPRK3

The other exterior boundaries retain the ambient state. The inlet is zone 7;
the remainder of the bottom boundary is zone 3.

## Running

```bash
makeneko euler_mach_800_jet.f90
mpirun -np 4 ./neko euler_mach_800_jet.case
```

The supplied mesh is intentionally coarse so the benchmark can be run on a
workstation. A publication-quality solution requires substantially more
resolution.
