# Pipeflow at Re_b=5300 and Re_tau = 1800.

## Mesh generation

There is a pre-generated nmsh file that can be used. To regenerate the mesh with
a different resolution or domain size, one needs gmsh and gmsh2nek, the latter
available as a tool in Nek5000. To modify the mesh parameters, you should edit
the .geo file, which serves as input to gmsh.

Once you are happy generate the gmsh mesh with
```
gmsh -3 ./turb_pipe.geo -order 2
```

The next stage is to convert the mesh to the .re2 format of Nek5000. Run
gmsh2nek and the following input when prompted:

```
Enter mesh dimension: 3
Input fluid .msh file name: turb_pipe
Do you have solid mesh ? (0 for no, 1 for yes) 0
 Enter number of periodic boundary surface pairs:
1
 input surface 1 and  surface 2  Boundary ID
1 2
 please give re2 file name:
turb_pipe

```
This will generate a turb_pipe.re2 file, which can be converted to Neko's format
with

```
rea2nbin turb_pipe.re2
```

## Running

To perform a DNS of the case set the polynomial order to 7 in the case file..
should be used. We have two case files, one called pipe.case using
"flow_rate_force" to drive the flow,  and one called pipe_source.case that
instead uses a constant source term. The source term values is based on
computing the shear stress as 2 * (2 * Re_tau/Re_B)^2, which assumes that the
bulk velocity is 1.

We use quite forgiving tolerances, and this case should be possible to run
efficiently on consumer GPUs in single precision, specify '--enable-real=sp'
when configuring.
