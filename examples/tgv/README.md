## Simulation of the three-dimensional Taylor--Green vortex (TGV)

In this case we simulate a 3D Taylor--Green vortex in a cubical box of size $2\pi$ in each direction. During the simulation, the total kinetic energy and enstrophy are computed. The resulting flow fields may be visualised using Paraview or VisIt by opening the field0.nek5000 file.

There are three different meshes provided, $8^3$, $32^3$ and $64^3$ elements. New meshes could be constructed using a suitable Nek5000 box file, using appropriate conversion.

### Cases

Both cases are run with the same executable, built from `tgv.f90`:

```
makeneko tgv.f90
```

- `tgv.case`, the initial condition

$$u = \sin(x)\cos(y)\cos(z), \quad v = -\cos(x)\sin(y)\cos(z), \quad w = 0$$

  is prescribed in the user file, in `initial_conditions`.

- `tgv_expression.case`, the same initial condition, written as one
  mathematical expression per velocity component directly in the case file:

  ```json
  "initial_condition": {
      "type": "expression",
      "value": [
          "U_0*sin(x)*cos(y)*cos(z)",
          "-U_0*cos(x)*sin(y)*cos(z)",
          "0"
      ]
  }
  ```

  The amplitude `U_0` is declared under `case.constants`, and is referred to
  by name in the expressions. The two cases produce the same flow.

Note that the user file is still needed for `tgv_expression.case`, since it
also rescales the mesh from $0..8$ to $-\pi..\pi$ and computes the kinetic
energy and enstrophy. Only the initial condition moves into the case file, the
`initial_conditions` hook registered in `tgv.f90` is simply not called.
