## Simulation of a wall-mounted cylinder
In this case simulate a cylinder mounted on a no-slip wall in a open channel. We provide 3 different case files that use the same mesh.

* cyl_bl_basic.case can be run with the usual'neko' executable
* cyl_bl_user.case uses the user file cyl_bl.f90 and the neko executable that is generated with makeneko.
* cyl_bl_rot.case also needs the the user file and changes the boundary conditions to simulate a rotating cylinder with dong outflow conditions.

The Reynolds number can probably be increased in these simulations to get some more interesting results.
