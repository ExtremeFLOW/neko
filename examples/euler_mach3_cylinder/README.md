# Mach-3 Euler flow around a cylinder

This case reproduces the physical setup of the Ryujin Mach-3 cylinder
benchmark shown in <https://www.youtube.com/watch?v=pPP26zelb0M>:

- two-dimensional inviscid compressible Euler flow;
- channel `[0, 4] x [-1, 1]`;
- cylinder diameter `0.5`, centred at `(0.6, 0)`;
- primitive free-stream state `(rho, u, v, p) = (1.4, 3, 0, 1)`;
- `gamma = 1.4` and final time `t = 5`;
- prescribed free stream at the inlet, natural supersonic outflow, and slip
  walls;
- 25,020 structured quadrilateral spectral elements.

The exact reference parameters are in Ryujin's
[`euler-mach3-cylinder-2d.prm`](https://github.com/conservation-laws/ryujin/blob/development/prm/benchmarks/euler-mach3-cylinder-2d.prm).
The reference run used continuous Q1 finite elements with about 2.36 million
grid points per conserved component. This Neko case intentionally follows the
requested element budget instead: 25,020 elements with polynomial order 3.

Neko's opt-in invariant-domain-preserving Euler path is enabled. It is the
closest in-tree counterpart to Ryujin's convex-limiting solver, but the two
codes do not use the same spatial discretization. Positivity and internal-energy
limiting remain enabled. The local entropy lower-bound limiter is disabled
because Neko's strong curved-wall slip projection does not preserve that nodal
bound. The conservative CFL target is `0.2` rather than Ryujin's `0.9`.

## Boundary zones

| Zone | Boundary | Neko condition |
| ---: | --- | --- |
| 1 | inlet, `x = 0` | prescribed free-stream primitive state |
| 2 | outlet, `x = 4` | natural (no strong pressure condition) |
| 3 | top, `y = 1` | slip |
| 4 | bottom, `y = -1` | slip |
| 7 | cylinder | slip |

## Build and run

The supplied `mach3_cylinder.nmsh` can be used directly. To regenerate it,
first build `contrib/gmsh2nek/gmsh2nek`, then run:

```console
cd examples/euler_mach3_cylinder/mesh
./generate_mesh.sh
```

Compile the user file and run from the case directory:

```console
makeneko mach3_cylinder.f90
mpirun -n 8 ./neko mach3_cylinder.json
```

The field writer saves density, pressure, entropy, and artificial viscosity at
intervals of `0.1` simulation-time units. Density or the magnitude of its
gradient is suitable for visualizing the bow shock and wake.
