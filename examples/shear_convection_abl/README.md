# Simple shear-convection atmospheric boundary layer (ABL) case with 512 and 32768 elements and N=7
Simulation of a shear-convection atmospheric boundary layer (ABL) in a domain of size $(5000.0, 5000.0, 2000.0) \rm{[m]}$. The case is conditioned with dimensional geostrophic wind $[10.0, 0.0, 0.0] \rm{[m/s]}$, air density $1.0\rm{[kg/m^3]}$ and dynamic viscosity $1e-10 \rm{[kg/(m\cdot s)]}$. This case couples the shear from the geostrophic wind and the vertical convection from the buoyancy and the interior-pointing heat flux at the surface.

The simulation is performed at an LES resolution if one uses the `genmeshbox` command below to generate the mesh, and Deardorff's SGS model [1] is adopted, which solves a scalar equation for turbulent kinetic energy to formulate the eddy viscosity. The wall model used in this example is MOST [2], and a sponge layer is applied at the top of the domain to prevent reflection of gravity waves down in to the domain.
More details about the setup can be found in [3].

`./genmeshbox 0 5000 0 5000 0 2000 32 32 32 .true. .true. .false.`

Reference:
[1] Deardorff, J.W. Stratocumulus-capped mixed layers derived from a
three-dimensional model. Boundary-Layer Meteorol 18, 495–527 (1980).
https://doi.org/10.1007/BF00119502
[2] R. B. Stull, An Introduction to Boundary Layer Meteorology. Springer, 1988. [Online]. Available: https://doi.org/10.1007/978-94-009-3027-8
[3] L. Huusko, T. Mukha, L. L. Donati, P. P. Sullivan, P. Schlatter, and G. Svensson, “Large Eddy Simulation of Canonical Atmospheric Boundary Layer Flows With the Spectral Element Method in Nek5000,” Journal of Advances in Modeling Earth Systems, vol. 17, no. 10, p. e2025MS005233, 2025. https://doi.org/10.1029/2025MS005233 