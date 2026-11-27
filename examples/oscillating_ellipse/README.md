# Calculating Torque using Multiple Reference Points
In this example, the torque calculation of an inclined ellipse body is performed using 3 different possible ways: once, the torque is calculated around a point which moves rigidly with the ellipse, once around the pivot point, and once around a fixed point in the domain.

The mesh in this example is kept almost rigid for up to 0.5 units away from the ellipse wall by setting a high value for gain.

## Mesh
To generate the mesh, first open the `generate_mesh.sh` script and set the correct paths for your `gmsh` and `gmsh2nmsh` executables at the top of the file. Once the paths are configured, execute the script in the `mesh` folder.
