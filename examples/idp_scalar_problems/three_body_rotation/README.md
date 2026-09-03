# Three-body rotation

This case adapts the scalar three-body rotation benchmark from section 6.2 of
Nazarov (2026) to exercise Neko's experimental Euler IDP path. The scalar
profile is stored as `rho - 1`, and the paper's rigid-rotation velocity
`beta = (-2 pi y, 2 pi x)` defines its flux. The momenta are reset to
`(rho - 1) beta` after every complete Forward-Euler timestep, so the constant
density shift is not transported.

This velocity projection is an explicit test harness. It also resets total
energy to a constant pressure, so the calculation is not an unforced Euler
solution. Because the current Euler IDP implementation supports only fully
periodic affine meshes, the unit-disk profile is embedded in the periodic box
`[-1, 1]^2 x [0, 1/12]`. The density support remains away from the periodic
seams.

The checked-in mesh has 24 by 24 by 1 elements and is periodic in all three
directions. Check the mesh, compile the user file, and run with:

```sh
mesh_checker box.nmsh
makeneko three_body_rotation.f90
mpirun -np 4 ./neko three_body_rotation.case
```

The default case uses polynomial degree 5, Forward Euler, target CFL 0.3,
entropy-viscosity coefficient `C_R = 0.5`, and relaxed local density bounds
with `C_rel(P5) = 5`. The relaxation factor follows Nazarov's benchmark;
`C_R = 0.5` gives a less diffusive solution while suppressing the visible
element-scale ripples.

The case selects `entropy_pair_type = user`. For the transported scalar
`s = rho - 1`, its callback supplies

```text
eta(s) = s^2/2,
q(s,x,y) = (-2 pi y eta, 2 pi x eta),
lambda(x,y) = 2 pi sqrt(x^2 + y^2).
```

Thus the entropy sensor follows the prescribed scalar transport equation
instead of applying the ideal-gas Euler entropy to this proxy state.

At the end, the user routine prints density bounds, one-rotation errors, mass
drift, and accumulated IDP diagnostics. It also writes one density-plane CSV
file per MPI rank. Run `python3 plot.py` with NumPy and Pillow available to
create `three_body_rotation.png`.

Forward Euler is intentional: projecting only after a complete SSPRK3 step
does not constrain its internal stage velocities. A faithful SSPRK3 version
would need to apply the prescribed-velocity projection inside every IDP stage.

The user wave speed drives the graph viscosity, and the correction is limited
only by relaxed local density bounds. The Euler internal-energy and entropy
constraints are disabled for this scalar proxy.
