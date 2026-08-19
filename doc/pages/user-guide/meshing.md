# Meshing {#meshing}
## General considerations {#general-considerations}

The Spectral Element Method (SEM) used in Neko is a high-order finite element
method. You generate a mesh of elements; the individual degrees of freedom
(DoFs) inside each element are then placed automatically by Neko based on the
chosen polynomial order.

In SEM, the total DoF count increases with the polynomial order `p`. Consequently,
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

This leads to the question: What `p` should guide the element spacing when
constructing the mesh? A default choice is 7; this is the unofficial standard in
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
tetrahedral mesh (usually with triangular prisms in the boundary layer) and then
convert it by chopping each tetrahedron into hexahedra. But keep in mind that
the quality of such meshes are generally not very high.

Strong skewness, non-orthogonality and other issues will likely cause your case
to crash. The SEM is a high-order method and that comes at the price of much
higher sensitivity to mesh quality than that of, say, a finite volume solver.
Since the SEM flavour used in Neko allows for discontinuities in derivatives
across element boundaries, the size of these jumps serves as a good proxy for
mesh quality. This can be assessed by simply inspecting some derivative field
visually. Similarly the primary solution fields will start exhibiting wiggles
on coarse grids, particularly towards the element boundaries. While these visual
indicators are not quantitative, they can give you a good idea of where in the
domain things have potentially gone wrong.

Finally, an important note for those simulating canonical flows is that periodic
boundary conditions are part of the mesh definition in Neko, and not something
you specify in the case file. Baking in periodicity into the mesh is discussed
below.

## Constructing meshes {#constructing-meshes}

### Box meshes {#constructing-meshes-box}

If your domain is a box, you can use Neko's built-in mesher called `genmeshbox`.
The tool can distribute elements along each axis according to user input.

### Converting a Nek5000 mesh {#constructing-meshes-nek5000}

If you have a mesh in Nek5000's `.re2` format, a Neko utility called `rea2nbin`,
can convert it to the native Neko format, `.nmsh`. This is currently the most
tested way of producing meshes. Since this boils down mesh generation for Neko
to mesh generation for Nek5000, it is useful to have a local copy of Nek5000 to
make use of its mesh conversion tools.

However, we provide two such tools under `contrib`. One is `gmsh2nek`.
Executing the `compile.sh` script in `contrib/gmsh2nek` will produce a
`gmsh2nek` executable. As the name hints, this is a converter from the Gmsh
`.msh` format to `re2`. In addition to just converting, the utility allows to
define periodic boundaries.

Since Gmsh is popular, a lot of meshing software supports its format. There is a
caveat though: there have been several versions of the `.msh` and `gmsh2nek`
supports only a particular flavour. Moreover, it expects the order of the mesh
to be 2, whereas most software will produce a linear mesh. The typical meshing
pipeline is therefore the following:

1. Generate your element mesh and save it to the `.msh` format.
2. Open the mesh in Gmsh, set the order to 2, and export the mesh in the legacy
   version 2 `.msh` format.
3. Convert it to `.re2` with `gmsh2nek`, defining periodic boundaries if needed.
4. Run `rea2nbin` to get the `.nmsh`.
5. Run Neko's `mesh_checker` utility on the generated `.nmsh`. This will create
   an internal representation of the mesh, just like during a simulation, and
   then output some statistics. If this goes smoothly, your mesh should be good
   to go.

Naturally, the element mesh can also be directly generated by Gmsh itself! For
example, from a `.geo` file. This is done in the `turb_pipe` example. The README
there contains much the same meshing instructions as the points above, and
provides concrete commands to run to get from a `.geo` to a `.nmsh`. Finally, as
already noted, `gmsh2nek` is not the only converter provided in Nek5000, and you
may find it useful to explore alternatives, depending on your workflow.

The other tool is `icem2re2.py`, which can convert ICEM/Fluent meshes to `re2`.
Since it is a Python tool, it requires no installation, but there are
package dependencies. We refer to the detailed documentation under
`contrib/icem2re2` for further details.

### Experimental workflows {#constructing-meshes-experimental}

Here we list several other workflows, which have not been excessively tested,
but could be useful to users.

- A converter from OpenFOAM&reg; meshes to `.nmsh` is available in the
  [foamToNeko](https://github.com/ExtremeFLOW/foamToNeko) repository. The
  converted mesh has to be composed of hexahedrals. The natural limitation is
  that the mesh remains linear (no high-order boundary representation).
- The following [fork of gmsh](https://github.com/timofeymukha/gmsh) contains
  a native `.nmsh` writer, which makes it possible to skip the cumbersome
  conversion process. Periodic boundaries are not supported, and should be
  added using the `create_periodic_zones` utility described below.
- The following [tet to hex
  converter](https://github.com/adigh23/tet-to-hex-generator) can be used as 
  part of the meshing workflow. It supports hybrid meshes with tetrahedra in the
  bulk and prismatic layers near the wall.

## High-order boundary representation {#high-order-boundary-representation}

The polynomial order in the case file determines the number of GLL points used
inside an element.  It does not, by itself, make a boundary curved: a mesh whose
elements are specified only by their corner vertices has a linear (bilinear in
2D or trilinear in 3D) geometric map.  A high-order boundary representation is
obtained when the coordinates of the GLL points on that boundary are modified
while the geometry is constructed.

Neko provides the following routes.

### Curves stored in the mesh file

Besides its vertex and boundary-zone records, an `.nmsh` file can contain a
curve record for an element.  Each record supplies curve data and a type for
each of its edges.  Neko currently constructs geometry for circular arcs
(``'C'`` in a Nek5000 `.re2` file) and mapped curves (``'m'``): the latter
uses an edge midpoint to construct a quadratic map.  The resulting coordinates
at the GLL points are generated when the function space is initialized.

This is the appropriate route when the boundary can be described by these
Nek5000 curve types.  The usual workflow is to retain the curve information in
the `.re2` mesh and convert it with `rea2nbin`; the curve block is then carried
into the `.nmsh` file.  In particular, simply exporting a high-order mesh from
another tool is not sufficient unless its converter writes Neko curve records.
The exact record layout and the currently supported type codes are described in
the [`.nmsh` format appendix](@ref nmsh-format).

### User-defined geometry deformation

For a boundary that is not an arc or mapped curve, install an `apply_deform`
procedure on the mesh from the user-file `mesh_setup` hook.  The hook is called
after the mesh is read but before the solvers create their coordinate fields;
the deformation procedure itself is called after Neko has generated the
linear and any file-defined curved geometry at all GLL points.  It can therefore
project the nodes of selected boundary facets onto an analytic surface or apply
another high-order mapping.

The user-file setup and registration are described under [runtime mesh
deformation](@ref user-file_user-mesh-setup).

Note that changing `msh%points` directly in `mesh_setup`, as in the simple
scaling examples in the user-file guide, is useful for affine transformations
and changes the element vertices before geometry generation.  On its own it does
not create a curved, high-order boundary; use a stored curve or `apply_deform`
for that purpose.

## Adding periodicity to an existing `.nmsh` {#adding-periodicity-to-an-existing-nmsh}

Sometimes the mesh is already available as an `.nmsh`, but some pairs of labeled
boundary zones should actually be periodic. In this case, the utility
`create_periodic_zones` under `contrib` can be used to rewrite the zone
information into a new mesh file.

The utility expects the periodic surfaces to be represented by pairs of labeled
zones with a one-to-one facet correspondence and a constant translation offset
between the two surfaces. A typical invocation looks like:

```text
create_periodic_zones input.nmsh output.nmsh "(1,2),(3,4)"
```

This command converts zones `1` and `2` into one periodic pair, and zones `3`
and `4` into another. Any labeled zones not mentioned in the pair list are kept
unchanged. Existing periodic zones in the mesh are also preserved.

The utility infers the translation vector for each pair from the facet centres
and then verifies the match using the facet corner coordinates. If your mesh is
very large or your coordinates are noisy, you can override the matching
tolerance:

```text
create_periodic_zones input.nmsh output.nmsh "(1,2)" --tol=1.0e-7
```

As with any mesh conversion step, you should always run `mesh_checker` on the
resulting mesh and make sure it reports no errors.
