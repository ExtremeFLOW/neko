# ALE Simulation with Multiple Moving Bodies
This example illustrates how to set up the ALE framework in Neko when more than one moving body is present in the simulation.

The motion of one of the bodies is set using the ALE user functions in `double_oscillating_cylinders.f90`. Additionally, this example illustrates how to print the rotation angle, position of the pivot, and velocity of the pivot for each body in the log file, using the `user_check` subroutine in the user file.

## Mesh
To generate the mesh, first open the `generate_mesh.sh` script and set the correct paths for your `gmsh` and `gmsh2nmsh` executables at the top of the file. Once the paths are configured, execute the script in the `mesh` folder.
