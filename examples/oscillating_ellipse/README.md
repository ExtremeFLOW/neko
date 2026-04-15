# Calculating Torque using Multiple Reference Points
In this example, the torque calculation of an inclined ellipse body is performed using 3 different possible ways: once, the torque is calculated around a point which moves rigidly with the ellipse, once around the pivot point, and once around a fixed point in the domain.

The mesh in this example is kept almost rigid for up to 0.5 units away from the ellipse wall by setting a high value for gain.

## Mesh
To generate the mesh, first open the `generate_mesh.sh` script and set the correct paths for your `gmsh`, `gmsh2nek`, and `rea2nbin` executables at the top of the file. Once the paths are configured, execute the script in the `mesh` folder.

- The `generate_mesh.sh` script requires the `gmsh2nek` version from Neko's `contrib/gmsh2nek/` directory. If you use a different version, executing `gmsh2nek` using the prompt values in the script may fail, and you must execute the commands manually.
