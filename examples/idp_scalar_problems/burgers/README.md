# Two-dimensional Burgers equation

This case adapts the nonlinear scalar benchmark from section 6.3 of Nazarov
(2026). The scalar equation is

```text
u_t + d_x (u^2 / 2) + d_y (u^2 / 2) = 0
```

on `(-0.25, 1.75)^2` up to `t = 0.75`. Initially, `u = 1` on
`[0.5, 1.5]^2`; otherwise `u = -0.75`. This is the state obtained by
evaluating the exact solution of Guermond and Popov at `t = 0`.

Euler density must remain positive, while this scalar changes sign. The case
therefore uses the exact change of variable `rho = u + 1`. After every complete
Forward Euler step, the user routine resets both horizontal momenta to
`u^2 / 2` and resets energy to maintain constant pressure. The next density
residual is consequently the published Burgers flux. This is an explicit test
harness for nonlinear density transport and the IDP limiter; it is not an
unforced Euler solution.

The benchmark is a Cauchy problem: the initial state is `-0.75` outside the
unit square. The checked-in 24 by 24 by 1 mesh therefore uses the larger
periodic box `[-1.3, 2.3]^2 x [-0.05, 0.05]`. Its constant exterior buffer
keeps the periodic boundary outside the domain of influence through
`t = 0.75`; plots may be cropped to the published `(-0.25, 1.75)^2` window.

Check the mesh and compile the user file with:

```sh
mesh_checker box.nmsh
makeneko burgers.f90
mpirun -np 4 ./neko burgers.case
```

The case follows the paper with final time `0.75`, target CFL `0.5`, and
entropy-viscosity coefficient `C_R = 4`. It uses polynomial degree 4, Forward
Euler, and strict local density bounds. The user wave speed drives the graph
viscosity, while the Euler internal-energy and entropy constraints are disabled
for the correction. Forward Euler is required because the
nonlinear flux is reset only after a complete timestep; an SSPRK method would
need the reset inside every IDP stage.

The user entropy pair is written in terms of `u = rho - 1`:

```text
eta(u) = u^2 / 2,
q(u) = (u^3 / 3, u^3 / 3),
lambda(u) = sqrt(2) |u|.
```

At the end, the user routine prints scalar extrema, shifted-density mass drift,
violations of the initial global scalar bounds, and accumulated IDP diagnostics.
It writes one scalar-plane CSV file per MPI rank. Run `python3 plot.py` with
NumPy and Pillow available to create `burgers.png`.
