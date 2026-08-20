# Two-dimensional Burgers equation

This case adapts the nonlinear scalar benchmark from section 6.3 of Nazarov
(2026). The scalar equation is

```text
u_t + d_x (u^2 / 2) + d_y (u^2 / 2) = 0
```

on `(-0.25, 1.75)^2` up to `t = 0.75`. Initially, `u = 1` when both
`|x - 0.5| <= 1` and `|y - 0.5| <= 1`; otherwise `u = -0.75`.

Euler density must remain positive, while this scalar changes sign. The case
therefore uses the exact change of variable `rho = u + 1`. After every complete
Forward Euler step, the user routine resets both horizontal momenta to
`u^2 / 2` and resets energy to maintain constant pressure. The next density
residual is consequently the published Burgers flux. This is an explicit test
harness for nonlinear density transport and the IDP limiter; it is not an
unforced Euler solution.

The checked-in mesh has 24 by 24 by 1 elements on the published domain. It is
non-periodic in x and y and periodic in z. The current Euler IDP implementation
does not yet treat physical boundaries, so this case is intentionally not
runnable with IDP enabled. It is retained as the target setup for the boundary
extension; Neko will stop with the fully-periodic-domain requirement until that
work is complete.

Check the mesh and compile the user file with:

```sh
mesh_checker box.nmsh
makeneko burgers.f90
```

The case follows the paper with final time `0.75`, target CFL `0.5`, and
entropy-viscosity coefficient `C_R = 4`. It uses polynomial degree 4, Forward
Euler, and relaxed local density bounds. Forward Euler is required because the
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
