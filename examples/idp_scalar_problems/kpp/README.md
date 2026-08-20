# KPP rotating wave

This case adapts the nonlinear scalar benchmark of Kurganov, Petrova, and
Popov from section 6.4 of Nazarov (2026). The scalar equation is

```text
u_t + d_x sin(u) + d_y cos(u) = 0
```

on `[-2, 2] x [-2.5, 1.5]` up to `t = 1`. The initial value is `3.5 pi`
inside the unit disk and `pi/4` outside.

The scalar is represented directly by the positive Euler density, `rho = u`.
After every complete Forward Euler step, the user routine resets momentum to
`(sin(rho), cos(rho), 0)` and resets energy to maintain constant pressure.
Consequently, the density residual at the next step has the KPP flux. This is
an explicit test harness for nonlinear density transport and the IDP limiter;
it is not an unforced Euler solution.

The local density bounds are built from Euler bar states. They preserve the
Euler invariant set but do not reduce to the scalar maximum-principle interval
`[pi/4, 3.5 pi]`. The proxy may therefore leave the initial scalar range even
when density-bound relaxation is disabled. The final diagnostic reports these
global excursions explicitly.

The standard KPP benchmark holds the exterior state at the physical boundary.
The current Euler IDP implementation instead requires a fully periodic affine
mesh. Opposing faces initially carry the same exterior state, but interaction
with the periodic seams remains a difference from the standard problem and
must be checked when interpreting the result.

The checked-in mesh has 24 by 24 by 1 elements and is periodic in all three
directions. Check the mesh, compile the user file, and run with:

```sh
mesh_checker box.nmsh
makeneko kpp.f90
mpirun -np 4 ./neko kpp.case
```

The default case uses polynomial degree 5, Forward Euler, target CFL 0.5,
entropy-viscosity coefficient `C_R = 10`, and relaxed local density bounds.
Forward Euler is required because the prescribed nonlinear flux is reset only
after a complete timestep; an SSPRK method would need the reset inside every
IDP stage.

The case selects `entropy_pair_type = user` and provides the convex scalar
entropy pair

```text
eta(u) = u^2/2,
q(u) = (u sin(u) + cos(u), u cos(u) - sin(u)),
lambda(u) = max(|cos(u)|, |sin(u)|) <= 1.
```

The implementation uses the safe constant bound `lambda = 1`. Neko forms the
BDF entropy residual and its normalization from these fields; the user routine
does not reproduce the entropy-viscosity algorithm itself.

At the end, the user routine prints density extrema, mass drift, violations of
the initial global scalar bounds, and accumulated IDP diagnostics. It writes
one density-plane CSV file per MPI rank. Run `python3 plot.py` with NumPy and
Pillow available to create `kpp.png`.
