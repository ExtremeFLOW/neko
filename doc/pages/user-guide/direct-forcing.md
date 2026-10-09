# Direct forcing immersed boundaries {#direct-forcing}

\tableofcontents

The `direct_forcing` source term represents a solid body inside the fluid mesh
without a body-fitted grid. The surface of the body is given as a triangulated
mesh (an STL file), the mesh is covered with Lagrangian markers, and at every
time step a volume force is added to the momentum equation that brings the
fluid velocity at the markers to the velocity of the body. The body is at rest
in the current implementation.

Compared to the [Brinkman](@ref case-file_fluid-source-term) term, which
penalises the velocity in a volume around the body with a resistance that has
to be tuned to the Reynolds number, the direct forcing acts on the surface
only, has no resistance parameter and enforces the no-slip condition to within
one grid spacing. The price is a set of interpolation and spreading operators
whose properties decide the quality of the interface. This page explains the
method and how to choose and tune these operators. The case-file parameters
are listed in the [case file](@ref case-file_direct-forcing) documentation.

## The forcing step {#direct-forcing_step}

Let \f$ \mathbf{x}_m \f$ be the markers with unit normals \f$ \mathbf{n}_m \f$
taken from the surface triangles. Each time step, with the velocity
\f$ \mathbf{u}^n \f$ at the start of the step, the term

1. interpolates the velocity to the markers, \f$ \mathbf{u}_m = (I \mathbf{u}^n)_m \f$,
2. forms the marker force that cancels it over one step,
   \f$ \mathbf{F}_m = -\mathbf{u}_m / \Delta t \f$, and
3. spreads the marker forces back onto the mesh with a spreading operator
   \f$ S \f$, so that the source term is

\f[
  \mathbf{f} = -\frac{g}{\Delta t}\, S\, I\, \mathbf{u}^n ,
\f]

where \f$ g \f$ is the `spread_gain`. The result is assembled to a continuous
field, optionally filtered, and added to the right-hand side as \f$ f^u \f$
(the fluid scheme multiplies it by the density).

The term feeds the velocity back with a gain of order \f$ 1/\Delta t \f$.
Extrapolating it in time together with the other explicit terms is unstable
for any useful gain, so the term is evaluated at the start of the step and
applied as computed, without the EXTk extrapolation. This is done by the term
itself and needs no setting in the case file.

## Interpolation {#direct-forcing_interpolation}

The `interpolation` keyword selects how the velocity is evaluated at the
markers.

- `spectral` (default): the velocity is evaluated through the polynomial
  basis of the element that contains the marker, i.e. the same global
  interpolation that probes and particles use. It is exact for the
  discrete velocity and introduces no smoothing. Markers are located
  with a point search at start-up; markers that fall exactly on an element
  face or outside the mesh are moved slightly, see
  [marker placement](@ref direct-forcing_markers).
- `idw`: Shepard (inverse distance weighted) interpolation from the nodes
  within a radius `interpolation_rmax` local grid spacings of the marker,
  with the kernel of the spread below. It smooths the velocity over the
  stencil and needs no point search.
- `adjoint`: Shepard interpolation with the weights modified so that the
  interpolation is the mass-weighted transpose of the kernel spread. With
  this choice the composite operator \f$ S I \f$ is symmetric and
  positive semi-definite in the mass inner product with all eigenvalues in
  \f$ [0, 1] \f$: the forcing can never add kinetic energy to the fluid, and
  the gain of the composite is exactly `spread_gain`. The interpolation
  radius is tied to the spread radius `rmax`. Only valid with the `idw`
  spread.

## Spread {#direct-forcing_spread}

The `spread` keyword selects how the marker forces are distributed onto the
mesh.

- `idw` (default): an inverse distance kernel. With \f$ h(\mathbf{x}) \f$
  the local grid spacing at a node and \f$ r = |\mathbf{x} - \mathbf{x}_m| /
  h(\mathbf{x}) \f$,

  \f[
     K(r) = \left( \frac{r_{\max} - r}{r_{\max}\, r} \right)^p
     \quad \text{for } r < r_{\max}, \qquad K(r) = 0 \text{ otherwise},
  \f]

  with `rmax` \f$ = r_{\max} \f$ and `power_parameter` \f$ = p \f$. The
  spread is the normalised average over the markers that reach a node,

  \f[
     (S \mathbf{F})(\mathbf{x}) = \frac{\sum_m K(r_m)\, \mathbf{F}_m}
                                      {\sum_m K(r_m)} ,
  \f]

  so adding markers does not amplify the forcing, and a constant marker
  velocity is reproduced exactly. The interface has a thickness of about
  `rmax` grid spacings on each side of the surface.
- `adjoint`: the transpose of the spectral interpolation, weighted with the
  element-lumped mass. Each marker deposits its force onto the nodes of
  its element with the values of the basis functions at the marker as
  weights, the contributions are assembled, and the result is divided by
  the lumped mass. The marker forces are scaled with per-marker weights
  \f$ D_m \f$ chosen so that the eigenvalues of \f$ S I \f$ lie in
  \f$ [0, g] \f$ (an absolute row-sum bound on the marker Gram matrix).
  The footprint is exactly the element the marker sits in, which gives the
  sharpest interface the discretisation allows. Pairs with the `spectral`
  interpolation only.
- `adjoint_mass`: as `adjoint` but with the pointwise GLL mass instead of
  the element-lumped one. It is the exact \f$ L^2 \f$ projection of the point
  forces but responds with \f$ 1/B \f$ at the nodes, which makes the corner
  nodes of a cut element react far more strongly than the interior ones.
  Kept for reference; `adjoint` is the recommended variant.

## Choosing a pairing {#direct-forcing_pairings}

| Interpolation | Spread    | Gain of \f$ S I \f$       | Interface width       | Use                                                                 |
| ------------- | --------- | ------------------------- | --------------------- | ------------------------------------------------------------------- |
| `spectral`    | `idw`     | \f$ g \f$ on the constant mode | about `rmax` spacings | Default. Robust, good general-purpose choice.                      |
| `adjoint`     | `idw`     | exactly \f$ g \f$         | about `rmax` spacings | Energy-stable. High Reynolds numbers, cases where the default rings. |
| `idw`         | `idw`     | not controlled            | about `rmax` spacings | No point search; cheapest start-up. Smooths the interpolated velocity. |
| `spectral`    | `adjoint` | at most \f$ g \f$, about \f$ 0.3 g \f$ on the constant mode | one element | Sharpest interface. Needs the face nudge; see the notes below. |

Start with the default pairing. If the velocity or pressure oscillates around
the interface at high Reynolds number, switch the interpolation to `adjoint`,
reduce `spread_gain` to about 0.5 and consider an elementwise filter. Use the
adjoint spread when the interface thickness of the kernel is the limiting
factor, for instance for thin bodies or coarse meshes relative to the body.

## Gain and stability {#direct-forcing_gain}

Because the forcing is proportional to \f$ 1/\Delta t \f$, its stability does
not depend on the time step. The product of the time step and the Jacobian of
the forcing is the gain \f$ g \| S I \| \f$, nothing else, so the time step
controller cannot help or harm the forcing. What the time step sets is the
velocity increment the flow adds at the surface per step, which the forcing
then has to remove.

The gain controls two things:

- **Residual slip.** A mode with eigenvalue \f$ \lambda \f$ of \f$ S I \f$
  is reduced by the factor \f$ 1 - g \lambda \f$ each step. With \f$ g = 1 \f$
  and the `adjoint` interpolation the slip on the markers is removed in one
  step. With a lower gain a residual slip remains that scales with the
  per-step increment times \f$ (1 - g)/g \f$, i.e. with the time step.
- **Damping of interface ringing.** Removing the whole slip in one step
  excites the pressure and the high modes of the element. The sheet model
  of the interface gives a per-step amplification of these modes of roughly
  0.76 at \f$ g = 1 \f$ and 0.64 at \f$ g = 0.5 \f$, so a lower gain is the
  cheapest remedy against oscillations at the interface.

Raising the gain above one is only meaningful for the adjoint spread, whose
constant mode gain is a fraction of \f$ g \f$. Two limits apply there. The time
integrator tolerates a gain of about 4.5 on any mode, and a node whose
velocity is reversed by more than twice its value in one step grows without
bound. The latter is the operative one: the start-up log reports the unit
response of the operator and the gain bound it implies, and warns when the
chosen gain exceeds it.

With adaptive time stepping the forcing stays stable, but at a gain below one
the residual slip follows every change of the time step, which shows as drift
in the force trace. For force measurements use a fixed time step, or cap the
controller with `max_dt`.

## Two sides of the surface {#direct-forcing_one-sided}

With `one_sided` set to `true` (the default), the nodes near the surface are
split into a plus side and a minus side by the tangent plane of the nearest
marker, and interpolation and spread are carried out for each side with that
side's nodes only. The stencil never straddles the surface, so the velocity
jump across a thin body is not smeared, and both sides are driven to rest
independently. The side a node belongs to is decided by the surface normals of
the STL, but since both sides are forced the orientation of the normals does
not affect the result; it only labels the plus and minus sides in the
diagnostics. With `one_sided` set to `false` a single interpolation and spread
over all nodes is used.

Nodes within `mask_band` grid spacings of the surface belong to both sides.
The adjoint spread needs this band (default 0.1) because a surface node has to
be reachable from both sides for the deposit to be normalised; the kernel
spread has a footprint on both sides of a face by construction and uses no
band by default. For the adjoint spread the one-sided normalisation of a
marker whose stencil has little weight on one side is clamped at
`one_sided_min_weight` and switched off for a side without weight, and the
start-up log reports how many markers were affected and where.

## Marker placement {#direct-forcing_markers}

Each `boundary_mesh` object is read from an STL file, optionally transformed
into a bounding box, and then refined to the resolution of the fluid mesh
before the markers are seeded. For each triangle the local grid spacing
\f$ h \f$ is the smallest spacing of the fluid elements under its centroid and
vertices, and the triangle is split into \f$ n^2 \f$ sub-triangles with

\f[
   n = \left\lceil \frac{\text{diameter}}{\texttt{marker\_spacing} \cdot h}
       \right\rceil ,
\f]

with one marker at the centroid of each sub-triangle. Sub-triangles outside
the fluid mesh are dropped. The refinement is what keeps a coarse STL from
leaking: a triangle larger than the elements it crosses would otherwise be
represented by a single marker and force only a disc around its centre. With
`refinement` set to `local` (the default) every triangle uses its own
\f$ h \f$; `uniform` uses the smallest spacing under the whole surface, which
on a stretched mesh can produce far more markers than needed. The number of
markers is reported in the log, and the term aborts if it exceeds
`max_markers`.

The spread radius has to cover the marker spacing. With the default
`marker_spacing` of one grid spacing, `rmax` of about 1.2 or larger closes the
surface; as a rule the kernel radius should be at least \f$ 1 + s/2 \f$ for a
marker spacing of \f$ s \f$ grid spacings.

With the spectral interpolation the markers are located by a point search. A
marker that is not found, that falls outside its element in reference
coordinates, or that lies on an element face while the surface runs along
that face is moved and searched again, up to three times. Unplaced and outside
markers are moved along the surface by `nudge` times the smallest grid
spacing, then inward. Markers on a face are pushed outward along the normal by
`face_nudge` local grid spacings. A marker exactly on a face is a degenerate
point for the adjoint spread: its basis function is a Kronecker delta on the
face, the wall it builds is one node thick, and the flow leaks through it as a
jet. A displacement of a few per cent of a grid spacing removes the
degeneracy, while a larger one only thickens the wall. `face_tolerance` sets
how close to a face, in reference coordinates, a marker counts as on it;
it defaults to 0.02 for the adjoint spread and is switched off for the kernel
spread.

## Element search and padding {#direct-forcing_padding}

The markers are matched to elements through a tree of axis-aligned bounding
boxes. Each element box is grown by `padding` times the element diameter, so
that a kernel stencil centred on a marker just outside an element still finds
it. The reach of the stencils is `rmax` (and `interpolation_rmax` for the
Shepard interpolations) local grid spacings, and the start-up log prints the
padding this needs on the active elements next to the value that is set. A
warning is issued if the stencils exceed the search reach; increase `padding`
to the reported value in that case.

## Filtering {#direct-forcing_filter}

The assembled forcing can be filtered before it is added to the right-hand
side with the `filter` object, which takes the same keywords as the
[filters](@ref filter) elsewhere in Neko. An `elementwise` filter of Boyd type
on the two or three highest modes damps the element-scale oscillations the
forcing excites at high Reynolds number, at the cost of a slightly thicker
interface. It matters most for the adjoint spread, whose footprint is a single
element. Use the `transfer` object rather than a `transfer_function` list, so
that the filter does not have to be rewritten when the polynomial order
changes; its defaults keep the modes that carry the element integral of the
forcing. A `PDE` filter smooths over a chosen length instead.

## Force on the body {#direct-forcing_force}

With a `force_output` object the term writes the force on the immersed objects
to a CSV file with the columns `tstep, time, Fx, Fy, Fz`. The force is minus
the momentum the forcing adds to the fluid,

\f[
   \mathbf{F} = -\rho \int \mathbf{f} \, dV ,
\f]

integrated over the assembled and filtered forcing, times `scale`. Setting
`scale` to \f$ 2 / (\rho U^2 A) \f$ gives the force coefficients directly. The
output interval follows the usual `output_control` and `output_value` pair
(`tsteps` or `simulationtime`), and the file is appended to by default so that
a restart keeps the earlier history.

Two caveats: the rate of change of the fluid momentum inside the body is not
included, and the forcing on nodes that carry a strong velocity boundary
condition is counted although the velocity solve discards it, which
over-reports the force on objects that touch such boundaries.

## Diagnostics {#direct-forcing_diagnostics}

The start-up block of the term lists the selected schemes with only the
parameters that apply to them, followed by

- `Minimum ds`, `Maximum ds`: the range of the local grid spacing.
- Per object: the file name, the number of STL points, the number of
  markers with the largest number of sub-markers per triangle edge, and the
  bounding box of the refined triangles.
- `Tot lagpts`: the total number of markers after clipping against the
  fluid mesh.
- `Markers not placed .., outside .., on a face ..`: the result of the point
  search before and after each nudge pass, with the bounding box of the
  affected markers (spectral interpolation).
- `Shared mrks`: markers held by several ranks (Shepard interpolations).
- `Mask zeros`: the number of nodes excluded from each side.
- `Interp +/- side: min nodes, mean, empty markers`: the stencil sizes of
  the Shepard interpolations; an empty marker has no node within its radius
  on that side and is a sign that `interpolation_rmax` is too small.
- For the adjoint spread: the unit response of one forcing step and the
  gain bound it implies, the constant-mode gain (minimum, mean and where the
  minimum sits), the largest gain from a power iteration, the number of
  zero-weight markers, the momentum ratio between the physical and the
  intended forcing, and the clamped or switched-off one-sided
  normalisations.
- `Padding needed for the kernel stencils`: see
  [padding](@ref direct-forcing_padding).

With `marker_output` set to `true` every rank writes its markers to
`df_markers_<rank>.csv` with the coordinates, normals, the lumped weights and
constant-mode gains of both sides, and the one-sided weight sums, which can be
loaded into ParaView as a point cloud. With the adjoint spread the fields
`ib_response_plus`, `ib_response_minus`, `ib_mask_plus` and `ib_mask_minus`
are also registered, and the [field writer](@ref simcomp_field_writer) can
write them.

## Troubleshooting {#direct-forcing_troubleshooting}

- **Flow leaks through the body.** Check the markers line in the log: with
  the refinement active the markers sit one grid spacing apart, so a leak
  means the spread does not cover the spacing. Increase `rmax` to 1.5 or
  reduce `marker_spacing`. With the adjoint spread, check the markers on a
  face line: a jet along an element face comes from markers on that face;
  `face_nudge` of 0.02 to 0.1 closes it.
- **Oscillations around the interface at high Reynolds number.** Use the
  `adjoint` interpolation, set `spread_gain` to 0.5, and add an elementwise
  Boyd filter on the highest modes.
- **Isolated spots of zero velocity away from the body.** This is the
  \f$ 1/B \f$ response of the `adjoint_mass` spread at element corners; use
  `adjoint` instead.
- **The run diverges and the CFL explodes.** A gain-type divergence does
  not respond to a smaller time step, since the forcing scales with
  \f$ 1/\Delta t \f$. Reduce `spread_gain`, and for the adjoint spread
  compare it with the node-stable gain bound in the log.
- **Padding warning.** Set `padding` to the value the log reports.
- **Too many markers.** The term aborts above `max_markers`. Switch
  `refinement` to `local` if it was set to `uniform`, or increase
  `marker_spacing` together with `rmax`.

## Performance notes {#direct-forcing_performance}

On device backends the interpolation, the spread, the assembly and the filter
run on the device for all pairings. The only per-step host traffic is marker
sized: the partial sums down and the marker velocities up. The point search,
the marker seeding, the masks and the weights are computed once at start-up on
the host, and the start-up cost grows with the number of markers and STL
triangles. The term reuses the gather-scatter of the fluid, so it adds no
communication set-up of its own, and it takes the four fields it needs per
step (the three forcing components and a work field) from the scratch
registry instead of holding them for the whole run.

## Example {#direct-forcing_example}

A sphere in a uniform stream with the default pairing and force output of the
drag and lift coefficients:

~~~~~~~~~~~~~~~{.json}
"source_terms": [
   {
      "type": "direct_forcing",
      "rmax": 1.2,
      "power_parameter": 0.5,
      "one_sided": true,
      "objects": [
         {
            "type": "boundary_mesh",
            "name": "sphere.stl"
         }
      ],
      "force_output": {
         "output_file": "sphere_force.csv",
         "output_control": "tsteps",
         "output_value": 10,
         "scale": 2.546
      }
   }
]
~~~~~~~~~~~~~~~

The same body at a high Reynolds number with the energy-stable pairing, a
reduced gain and an order-independent Boyd filter:

~~~~~~~~~~~~~~~{.json}
"source_terms": [
   {
      "type": "direct_forcing",
      "interpolation": "adjoint",
      "spread": "idw",
      "rmax": 1.5,
      "power_parameter": 0.5,
      "spread_gain": 0.5,
      "filter": {
         "type": "elementwise",
         "elementwise_filter_type": "Boyd",
         "transfer": { "cutoff": 0.4, "strength": 3.5, "order": 1 }
      },
      "objects": [
         {
            "type": "boundary_mesh",
            "name": "sphere.stl"
         }
      ]
   }
]
~~~~~~~~~~~~~~~

And the sharp-interface pairing:

~~~~~~~~~~~~~~~{.json}
"source_terms": [
   {
      "type": "direct_forcing",
      "interpolation": "spectral",
      "spread": "adjoint",
      "spread_gain": 1.0,
      "face_nudge": 0.02,
      "marker_output": true,
      "objects": [
         {
            "type": "boundary_mesh",
            "name": "sphere.stl"
         }
      ]
   }
]
~~~~~~~~~~~~~~~
