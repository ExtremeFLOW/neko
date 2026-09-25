# Setting up the ALE framework
This example illustrates how to set up the Arbitrary Lagrangian-Eulerian (ALE) framework in Neko to simulate a simple oscillating cylinder.

## Mesh
To generate the mesh, first open the `generate_mesh.sh` script and set the correct paths for your `gmsh`, `gmsh2nek`, and `rea2nbin` executables at the top of the file. Once the paths are configured, execute the script in the `mesh` folder.

- The `generate_mesh.sh` script requires the `gmsh2nek` version from Neko's `contrib/gmsh2nek/` directory. If you use a different version, executing `gmsh2nek` using the prompt values in the script may fail, and you must execute the commands manually.
