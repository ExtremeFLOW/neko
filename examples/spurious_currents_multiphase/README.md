# Spurious Currents Multiphase Test Case

Stationary circular drop (D = 0.4) in a 1x1 box with slip boundaries. Equal densities (rho = 300) and viscosities (mu = 0.1), surface tension sigma = 1.0, giving Laplace number La = 12000.

Generate the mesh with

```bash
genmeshbox 0 1 0 1 -0.1 0 10 10 1 .false. .false. .true.
```

Run with

```bash
makeneko spurious_currents.f90
mpirun -np 4 ./neko spurious_currents.case
```
