# Taylor-Green Vortex (TGV)

Classical test case for studying transition to turbulence and energy dissipation in compressible flows.

## Problem Description

Initial conditions for velocity and pressure in [0,2π]³ periodic domain:

```python
u = sin(x)cos(y)cos(z)
v = -cos(x)sin(y)cos(z) 
w = 0
p = p0 + ρ0V0²/16 * (cos(2x) + cos(2y))(2 + cos(2z))
```

Parameters:

- Mach number = 1.25
- Initial density ρ0 = 1.0
- Reference velocity V0 = 1.0
- Reference pressure p0 = ρ0V0²/(γMa²)
- γ = 1.4

## Running the case

1. Generate periodic mesh:

```bash
genmeshbox 0 6.28318530718 0 6.28318530718 0 6.28318530718 20 20 20 .true. .true. .true.
```

2. Run the simulation:

```bash
/path/to/neko_install/bin/makeneko euler_tgv.f90
mpirun -n 4 neko euler_tgv.case
```

## Analysis
Additional tools will be added later for this purpose.

Monitor:

- Kinetic energy decay
- Enstrophy evolution
- Energy spectra

Expected features:

- Initial vortex breakdown
- Development of small-scale structures
- Energy cascade to smaller scales

## Reference

Taylor, G.I. & Green, A.E. (1937) "Mechanism of the Production of Small Eddies from Large Ones" Proc. R. Soc. Lond. A 158, 499-521
