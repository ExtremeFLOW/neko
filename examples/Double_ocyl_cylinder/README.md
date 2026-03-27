This example illustrates how to set up the ALE framework in Neko when more than one moving body is present in the simulation.

The motion of one of the bodies is set using the ALE user functions in `cylinder.f90`. Additionally, this example illustrates how to print the rotation angle, position of the pivot, and velocity of the pivot for each body in the log file, using the `user_check` subroutine in the user file.

To generate the mesh, execute `generate_mesh.sh` in the mesh folder: