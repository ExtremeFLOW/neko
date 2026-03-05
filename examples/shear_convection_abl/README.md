# Simple shear-convection atmospheric boundary layer (ABL) case with 512 and 32768 elements and N=7
Simulation of a shear-convection atmospheric boundary layer (ABL) in a domain of size $(5000.0, 5000.0, 2000.0) \rm{[m]}$. The case is conditioned with dimensional geostrophic wind $[10.0, 0.0, 0.0] \rm{[m/s]}$, air density $1.0\rm{[kg/m^3]}$ and dynamic viscosity $1e-10 \rm{[kg/(m\cdot s)]}$. This case couples the shear from the geostrophic wind and the vertical convection from the buoyancy and the interior-pointing heat flux at the surface.

The simulation is performed at an LES resolution if one uses the `genmeshbox` command below to generate the mesh, and Deardorff's model [1] is adopted, which solves a scalar equation for turbulent kinetic energy to formulate the eddy viscosity.

`./genmeshbox 0 5000 0 5000 0 2000 32 32 32 .true. .true. .false.`

Reference:

[1] Deardorff, J.W. Stratocumulus-capped mixed layers derived from a
three-dimensional model. Boundary-Layer Meteorol 18, 495–527 (1980).
https://doi.org/10.1007/BF00119502