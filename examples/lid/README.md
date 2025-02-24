## Simulation of a lid-driven cavity
In this case we simulate a lid-driven cavity with smoothened lid (belt) velocity to fulfil the continuity equaiton in the corners. There is a 3D example (with 3 elements in the spanwise $z$ direction) and a pseudo-2D example (with a 2D mesh that then generates a case with one element in the spanwise direction). During the simulation, the kinetic energy and enstrophy are calculated and written out in a csv output.

The `lid.box` file needs to be converted to a .re2 file using the Nek5000 tool genbox, and then using rea2nbin into the `lid.nmsh` Neko mesh file.

The cavity flow is stable (steady) up to a Reynolds number of about 8000. The mesh size may be changed in the box file.

This example shows how to read a run-time parameter from the json case file into the user file. 
