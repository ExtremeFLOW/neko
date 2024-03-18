## Tollmien-Schlichting (TS) wave test in channel flow.
This test refers to [1].
The temporal instability of TS wave in a channel flow at Re_c=5000 is tested in this example.
Domain Size: (x, y, z) = (2pi/1.12, 2, 2pi/2.1), where x and z are two homogeneous derections.
A recommended number of elements: (x, y, z) = (36, 36, 36) for decaying case and (16, 16, 16) for transitional case.
A 2D (k_z=0) TS wave and two 3D (k_z=\pm 1) oblique TS waves are settled up as the initial condition in a channel with 1 period on streamwise and spanwise direction.
One could use contrib/map_to_equidistant_1d twice to interpolate the fields onto a mesh whose grid points are equidistant in x- and z- directions. And then Fourier transform could be performed to get the amplitude of the TS wave.

Reference:
[1] Schlatter, P.C., 2005. Large-eddy simulation of transition and turbulence in wall-bounded shear flow (Doctoral dissertation, ETH Zurich).
