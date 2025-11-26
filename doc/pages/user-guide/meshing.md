# Meshing {#meshing}
## General considerations {#general-considerations}

The Spectral Element Method (SEM) used in Neko is a high-order finite element
method. You generate a mesh of elements; the individual degrees of freedom
(DoFs) inside each element are then placed automatically by Neko based on the
chosen polynomial order.

In SEM, the total DoF count scales with the polynomial order `p`. Consequently,
the effective spatial resolution depends both on the element size and on `p`.
Within each element, nodal points are placed at Gauss–Lobatto–Legendre (GLL)
locations, which are clustered near element boundaries and not uniformly spaced.

Below are ASCII sketches of approximate GLL node distributions in a 1D element
for `p = 4` and `p = 6`. Each sketch spans `[−1, 1]` in the reference element.
Vertical ticks show node positions; note the boundary clustering.

```
-1                                        1 
|------|-------------|-------------|------|  p = 4
^      ^             ^             ^      ^
|      |             |             |      |
x1     x2            x3            x4     x5

-1                                        1
|--|------|----------|----------|------|--|  p = 6
^  ^      ^          ^          ^      ^  ^
|  |      |          |          |      |  |
x1 x2     x3         x4         x5     x6 x7
```

The leads to the question: What `p` to aim the element spacing for, when
constructing the mesh? A default choice is 7, this is the unofficial standard in
the SEM community. If you plan a mesh refinement study, it is useful to think if
you can make use of `p`-refinement. That is, use the same element mesh, and
increase resolution by increasing polynomial order. This saves you a lot of time
because you do not need to re-mesh. Also, projecting the solution onto a new
polynomial order is often seamless in Neko. For example, you can easily use a
field generated with `p` equal to 5 as initial conditions for a run with a
higher order. A possible strategy is then to aim the element size to get a
"good" mesh with order 7, and then treat orders 5, 9 and 11 as refinement
levels.

Neko is a 3D code. But you can provide it with a 2D mesh in the x-y plane. In
this case it will extrude it to contain 1 element in z, and periodic conditions
will be applied on the sides. This allows to generate pseudo-2d solutions. You
can use some tricks, like setting the z-component velocity to zero at each time
step in the user file. There are a few pseudo-2d examples available, one is
`lid`.

All elements of a Neko mesh must be hexahedral (or quadrilateral for the
pseudo-2d case). The topology can be fully unstructured. A practical way of
generating a hexahedral mesh of a complicated geometry is to first generate a
tetrahedral mesh (usually with triangle prisms in the boundary layer) and then
convert it by chopping each tetrahedron into hexahedra. But keep in mind that
the quality of such meshes are generally no very high.

Strong skewness, non-orthogonality and other issues will likely cause your case
to crash. The SEM is a high-order method and that comes at the price of much
higher sensitivity to mesh quality than that of, say, a finite volume solver.
Since the SEM flavour used in Neko allows for discontinuities in derivatives
across element boundaries, the size of these jumps serve as a good proxy for
mesh quality. This can be assessed by simply inspecting some derivative field
visually. Similarly the primary solution fields will start exhibiting wiggles
on coarse grids, particularly towards the element boundaries. While these visual
indicators are not quantitative, then can give you a good idea of where in the
domain things have potentially gone wrong.

Finally, an important note for those simulating canonical flows is that periodic
boundary conditions are part of the mesh definition in Neko, and not something
you specify in the case file. Baking in periodicity into the mesh is discussed
below.

## Constructing meshes
If your domain is a box, you can use Neko's built-in mesher called `genmeshbox`.
Otherwise, you have to convert your mesh into a Nek5000 format called `.re2`,
and then apply a utility in Neko called `rea2nbin`, which produces a mesh in the
native Neko format `.nmsh`. 

So, for most practical cases, the mesh generation for Neko boils down to mesh
generation for Nek5000. At this stage it is indeed helpful to have a local copy
of Nek5000 to make use of its mesh conversion tools. However, we provide one
such tool under `contrib/gmsh2nek`.  Executing the `compile.sh` script in that
folder will produce a `gmsh2nek` executable. As the name hints, this is a
convertor from the gmsh `.msh` format to `re2`. In addition to just converting,
the utility allows to define periodic boundaries.

Since gmsh is popular, a lot of meshing software supports its format. There is a
caveat though: there have been several versions of the `.msh` and `gmsh2nek`
supports only a particular flavour. Moreover, it expects the order of the mesh
to be 2, whereas most software will produce a linear mesh. The typical meshing
pipeline is therefore the following:

1. Generate your element mesh and save it to the `.msh` format. 
2. Open the mesh in gmsh, set the order to 2, and export the mesh in the legacy
   version 2 `.msh` format. 
3. Convert it to `.re2` with `gmsh2nek`, defining periodic boundaries if needed.
4. Run `rea2nbin` to get the `.nmsh`.
5. Run Neko's `mesh_checker` utility on the generated `.nmsh`. This will create
   an internal representation of the mesh, just like during a simulation, and 
   then output some statistics. If this goes smoothly, your mesh should be good
   to go.

Naturally, the element mesh can also be directly generated by `gmsh` itself! For
example, from a `.geo` file. This is done in the `turb_pipe` example. The README
there contains much the same meshing instructions as the points above, and
provides concrete commands to run to get from a `.geo` to a `.nmsh`. Finally, as
already noted, `gmsh2nek` is not the only converter provided in Nek5000, and you
may find it useful to explore alternatives, depending on your workflow.
