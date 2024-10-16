## Tollmien-Schlichting (TS) wave test in channel flow.
This test refers to [1].
The temporal instability of TS wave in a channel flow at $Re_c=5000$ is tested in this example.
Domain Size: (x, y, z) = (2pi/1.12, 2, 2pi/2.1), where x and z are two homogeneous directions.

A recommended number of elements: (x, y, z) = (36, 36, 36) for decaying case and (16, 16, 16) for transitional case.
The mesh included in the examples has 16^3=4096 elements and was generated with `genmeshbox 0 2 -1 1 0 1 16 16 16 .true. .false. .true.`, see `genmeshbox` for more detials on generating box meshes.

A 2D ($k_z=0$) TS wave and two 3D ($k_z=\pm 1$) oblique TS waves are settled up as the initial condition in a channel with 1 period on streamwise and spanwise direction. The initial modes are plotted in the jupyter notebook plot_modes.ypynb. The transition can be tracked by following the forces on the upper and lower wall of the channel.

After running the simulation, one can use contrib/map_to_equidistant_1d twice (once in x and once in z) to interpolate the fields onto a mesh whose grid points are equidistant in x- and z- directions. One can then apply the Fourier transform to get the amplitude of the TS waves as the develop in time on this equidistant grid.

# Simulation components
In this case lambda2 and the force on the top and bottom walls are calculated. 

Reference:
[1] Schlatter, P.C., 2005. Large-eddy simulation of transition and turbulence in wall-bounded shear flow (Doctoral dissertation, ETH Zurich).
