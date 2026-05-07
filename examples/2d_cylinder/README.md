## Simulation of a pseudo-2D flow past a circular cylinder at Re=160.

While the mesh is 2D, Neko currently extrapolates it to a 3D-mesh that is one element thick with periodic BCs in the extrapolated direction (z). To avoid possible 3D instability the spanwise velocity component is set to zero after each time step.
