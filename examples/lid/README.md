## Simulation of a lid-driven cavity
In this case we simulate a lid-driven cavity with smoothened belt velocity to fulfil the continuity equaiton in the corners. There is a 3D example (with 3 elements in the spanwise $z$ direction) and a pseudo-2D example (with a 2D mesh that then generates a case with one element in the spanwise direction).

The meshes can be generated using `genmeshbox`. For the 2D example use
```
genmeshbox 0 1 0 1 0 1 6 6 1 .false. .false. .true. elem_dist.csv elem_dist.csv uniform
```
and for the 3D
```
genmeshbox 0 1 0 1 0 3 6 6 3 .false. .false. .true. elem_dist.csv elem_dist.csv uniform`
```
These commands generate a `box.nmsh`, which you can rename to `lid.nmsh` or `lid2d.nmsh`
depending on the targetted case.

The cavity flow is stable (steady) up to a Reynolds number of about 8000. The mesh size may be changed in the box file. There is also a short postprocessing script to show how to plot the simulation data.
