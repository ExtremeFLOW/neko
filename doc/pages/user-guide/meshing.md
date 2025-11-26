# Generating meshes {#meshing}
# General considerations {#general-considerations}

The Spectral Element Method (SEM) used in Neko is a high-order finite element
method. You generate a mesh of elements; the individual degrees of freedom
(DoFs) inside each element are then placed automatically by Neko based on the
chosen polynomial order.

In SEM, the total DoF count scales with the polynomial order `p`. Consequently,
the effective spatial resolution depends both on the element size and on `p`.
Within each element, nodal points are placed at Gauss–Lobatto–Legendre (GLL)
locations, which are clustered near element boundaries and not uniformly spaced.

Below are ASCII sketches of aproximate GLL node distributions in a 1D element
for `p = 4` and `p = 6`. Each sketch spans `[−1, 1]` in the reference element.
Vertical ticks show node positions; note the boundary clustering.

```
-1                                        1 
|-------|------------|------------|-------|  p = 4
^       ^            ^            ^       ^
|       |            |            |       |
x1      x2           x3           x4      x5

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
because you do not need to remesh. Also, projecting the solution onto a new
polynomial order is often seamless in Neko. For example, you can easily use a
field generated with `p` equal to 5 as initial conditions for a run with a
higher order. 

Neko is a 3D code. But you can provide it with a 2D mesh in the x-y plane. In
this case it will extrude it to contain 1 element in z, and periodic conditions
will be applied on the sides. This allows to generate pseudo-2d solutions. You
can use some tricks, like setting the z-component velocity to zero at each time
step in the user file.

All elements of a Neko mesh must be hexahedral (or quadrilateral for the
pseudo-2d case). The topology can be fully unstructured. As usual, use
high-quality structured meshes to produce the best results. A practical way of
generating a hexahedral mesh of a complicated geometry is to first generate
a tetrahedral mesh (usually with triangle prisms in the boundary layer) and then
convert it by chopping each tetrahedron into hexahedra.

# Generating meshes
If your domain is a box, you can use a built-in mesher called `genmeshbox`.
Otherwise, you have to convert your mesh into a Nek5000 format called `.re2`,
and then apply a utility in Neko called `rea2nbin`, which produces a mesh in the
native Neko format `.nmsh`. So, for most practical cases, the mesh
generation boils down to getting into the `.re2` format. 



Key takeaways:
- You mesh elements; Neko places GLL DoFs per element given `p`.
- Higher `p` increases DoFs per element and effective resolution.
- GLL nodes cluster near element boundaries; spacing is non-uniform.
