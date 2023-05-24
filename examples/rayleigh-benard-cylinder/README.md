#RBC in a cylindrical cell with aspect ratio 1:10
Small case for Rayleigh-Benard convection in a cylindrical cell at Ra=1e8. If there is no asymmetry in the velocity magnitude of the x/y-velocities early in the simulation the case will likely not trigger. Even if the amplitude of those are small they are a good indication of if the simulation will trigger at a later point.

From my experience, with these ics, this case triggers around T=30.

The mesh is originally a Nek5000 mesh. It was generated from a 2d gmsh mesh that was first converted to the Nek5000 format re2, then it was converted to rea with re2torea and then extruded in 3d with n2to3 with 'w' and 'v' at the top and bottom. Really it should be 'w' on both sides, but sine we set uinf = 0.0 the 'v' is really a wall as well.

To set the scalar bc at the top and bottom we set the scalar user bc in the user/f90 file which applies it on all zones with 'hardcoded' bcs.
