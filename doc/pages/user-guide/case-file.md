# Case File {#case-file}

\tableofcontents

The case file defines all the parameters of a simulation.
The format of the file is JSON, making it easy to read and write case files
using the majority of the popular programming languages.
JSON is heirarchical and, and consists of parameter blocks enclosed in curly
braces.
These blocks are referred to as objects.
The case file makes use objects to separate the configuration of different parts
 of the solver.
We refer the reader to the examples shipped with the code to get a good
idea of how a case file looks.
The table below provides a complete reference for all possible configuration
choices.

## High-level structure
The current high-level structure of the case file is shown below.

~~~~~~~~~~~~~~~{.json}
{
    "version": 1.0
    "case": {
        "numerics": {}
        "fluid": {}
        "scalar": {}
        "statistics": {}
        "simulation_components" : []
        "point_zones" : []
    }
}
~~~~~~~~~~~~~~~
The `version` keywords is reserved to track changes in the format of the file.
The the subsections below we list all the configuration options for each of the high-level objects.
Some parameters will have default values, and are therefore optional.

## Output frequency control
A common scheme for controlling the output frequency is applied for various
outputs.
It is described already now in order to clarify the meaning of several
parameters found in the tables below.

The frequency is controlled by two paramters, ending with `_control` and
`_value`, respectively.
The latter name is perhaps not ideal, but it is somewhat difficult to come up
with a good one, suggestions are welcome.

The `_value` parameter is a number, that defines the output frequency, but the
interpretation of that number depends on the choice of `_control`.
The three following options are possible.
1. `simulationtime`, then `_value` is the time interval between the outputs.
2. `tsteps`, then `_value` is the number of time steps between the outputs.
3. `nsamples`, then `_value` is the total number of outputs that will be
   performed in the course of the simulation.
4. `never`, then `_value` is ignored and output is never performed.


## The case object

This object is mostly used as a high-level container for all the other objects,
but also defines several parameters that pertain to the simulation as a whole.

| Name                       | Description                                                                                           | Admissible values                               | Default value |
| -------------------------- | ----------------------------------------------------------------------------------------------------- | ----------------------------------------------- | ------------- |
| `mesh_file`                | The name of the mesh file.                                                                            | Strings ending with `.nmsh`                     | -             |
| `output_boundary`          | Whether to write a `bdry0.f0000` file with boundary labels. Can be used to check boundary conditions. | `true` or `false`                               | `false`       |
| `output_directory`         | Folder for redirecting solver output. Note that the folder has to exist!                              | Path to an existing directory                   | `.`           |
| `output_precision`         | Whether to output snapshots in single or double precision                                             | `single` or `double`                            | `single`      |
| `load_balancing`           | Whether to apply load balancing.                                                                      | `true` or `false`                               | `false`       |
| `output_partitions`        | Whether to write a `partitions.vtk` file with domain partitioning.                                    | `true` or `false`                               | `false`       |
| `output_checkpoints`       | Whether to output checkpoints, i.e. restart files.                                                    | `true` or `false`                               | `false`       |
| `checkpoint_control`       | Defines the interpretation of `checkpoint_value` to define the frequency of writing checkpoint files. | `nsamples`, `simulationtime`, `tsteps`, `never` | -             |
| `checkpoint_value`         | The frequency of sampling in terms of `checkpoint_control`.                                           | Positive real or integer                        | -             |
| `checkpoint_format`        | The file format of checkpoints                                                                        | `chkp` or `hdf5`                                | `chkp`        |
| `restart_file`             | checkpoint to use for a restart from previous data                                                    | Strings ending with `.chkp`                     | -             |
| `timestep`                 | Time-step size                                                                                        | Positive reals                                  | -             |
| `variable_timestep`        | Whether to use variable dt                                                                            | `true` or `false`                               | `false`       |
| `max_timestep`             | Maximum time-step size when variable time step is activated                                           | Positive reals                                  | -             |
| `target_cfl`               | The desired CFL number                                                                                | Positive real                                   | `0.4`         |
| `cfl_max_update_frequency` | The minimum interval between two time-step-updating steps in terms of time steps                      | Integer                                         | `0`           |
| `cfl_running_avg_coeff`    | The running average coefficient `a` where `cfl_avg_new = a * cfl_new + (1-a) * cfl_avg_old`           | Positive real between `0` and `1`               | `0.5`         |
| `max_dt_increase_factor`   | The maximum scaling factor to increase time step                                                      | Positive real greater than `1`                  | `1.2`         |
| `min_dt_decrease_factor`   | The minimum scaling factor to decrease time step                                                      | Positive real less than `1`                     | `0.5`         |
| `end_time`                 | Final time at which the simulation is stopped.                                                        | Positive reals                                  | -             |
| `job_timelimit`            | The maximum wall clock duration of the simulation.                                                    | String formatted as HH:MM:SS                    | No limit      |

<h3> Boundary type numbering in the `output_boundary` field </h3>

When the `output_boundary` setting is set to `true`, and additional `.fld` file
will be stored in the beginning of the simulation, where the recognized
boundaries will be marked with an integer number. This is a good way to debug
the simulation setup. The value of the number depends on the type of the
boundary as follows:

1. A wall boundary, i.e. the `w` label.
2. A Dirichlet boundary, i.e. the `v` label.
3. An outlet boundary, i.e. the `o` label.
4. A symmetry boundary, i.e. the `sym` label.
5. An wall-normal transpiration boundary, i.e. the `on` label.
6. A periodic boundary.

Note that the boundary conditions can be both prescribed via the labels in the
case file or built into the mesh via conversion from a `.re2` file. Both types
will be picked up and marked in the field produced by `output_boundary`.


## Numerics
Used to define the properties of the numerical discretization.

| Name                         | Description                                                                                                   | Admissible values          | Default value                   |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------- | -------------------------- | ------------------------------- |
| `polynomial_order`           | The oder of the polynomial basis.                                                                             | Integers, typically 5 to 9 | -                               |
| `time_order`                 | The order of the time integration scheme. Refer to the `time_scheme_controller` type documention for details. | 1,2, 3                     | -                               |
| `dealias`                    | Whether to apply dealiasing to advection terms.                                                               | `true` or `false`          | `false`                         |
| `dealiased_polynomial order` | The polynomial order in the higher-order space used in the dealising.                                         | Integer                    | `3/2(polynomial_order + 1) - 1` |

## Fluid

The configuration of the fluid solver and the flow problem.
Contains multiple subobjects for various parts of the setup.

### Material properties
As per the governing equations, Neko requires the value of the density and
dynamic viscosity to define the flow problem. These can be provided as `rho` and
`mu` in the case file.

Alternatively, one may opt to provide the Reynolds number, `Re`, which
corresponds to a non-dimensional formulation of the Navier-Stokes equations.
This formulation can effectively be obtained by setting \f$ \rho = 1 \f$ and \f$
\mu = 1/Re \f$. This is exactly what Neko does under the hood, when `Re` is
provided in the case file.

Note that if both `Re` and any of the dimensional material properties are
provided, the simulation will issue an error.

As an alternative to providing material properties in the case file, it is
possible to do that in a special routine in the user file. This is demonstrated
in the `rayleigh-benard-cylinder` example. Ultimately, both `rho` and `mu` have
to be set in the subroutine, but it can be based on arbitrary computations and
arbitrary parameters read from the case file. Additionally, this allows to
change the material properties in time.

### Boundary types {#case-file_boundary-types}
The optional `boundary_types` keyword can be used to specify boundary conditions.
The reason for it being optional, is that some conditions can be specified
directly inside the mesh file.
In particular, this happens when Nek5000 `.re2` files are converted to `.nmsh`.
Periodic boundary conditions are *always* defined inside the mesh file.

The value of the keyword is an array of strings, with the following possible
values:

* Standard boundary conditions
  * `w`, a no-slip wall.
  * `v`, a velocity Dirichlet boundary.
  * `sym`, a symmetry boundary.
  * `o`, outlet boundary.
  * `on`, Dirichlet for the boundary-parallel velocity and homogeneous Neumann for
   the wall-normal. The wall-parallel velocity is defined by the initial
   condition.

* Advanced boundary conditions
  * `d_vel_u`, `d_vel_v`, `d_vel_w` (or a combination of them, separated by a 
  `"/"`), a Dirichlet boundary for more complex velocity profiles. This boundary
  condition uses a [more advanced user
  interface](#user-file_field-dirichlet-update).
  * `d_pres`, a boundary for specified non-uniform pressure profiles, similar in
  essence to `d_vel_u`,`d_vel_v` and `d_vel_w`. Can be combined with other
  complex Dirichlet conditions by specifying e.g.: `"d_vel_u/d_vel_v/d_pres"`.
  * `o+dong`, outlet boundary using the Dong condition.
  * `on+dong`, an `on` boundary using the Dong condition, ensuring that the
   wall-normal velocity is directed outwards.

In some cases, only some boundary types have to be provided.
For example, when one has periodic boundaries, like in the channel flow example.
In this case, to put the specification of the boundary at the right index,
preceding boundary types can be marked with an empty string.
For example, if boundaries with index 1 and 2 are periodic, and the third one is
a wall, we can set.
```
"boundary_types": ["", "", "w"]
```

### Inflow boundary conditions {#case-file_fluid-if}
The object `inflow_condition` is used to specify velocity values at a Dirichlet
boundary.
This does not necessarily have to be an inflow boundary, so the name is not so
good, and will most likely be changed along with type changes in the code.
Since not all cases have Dirichlet boundaries (note, the special case of a
no-slip boundary is treated separately in the configuration), this object
is not obligatory.
The means of prescribing the values are controlled via the `type` keyword:

1. `user`, the values are set inside the compiled user file.
2. `uniform`, the value is a constant vector, looked up under the `value`
   keyword.
3. `blasius`, a Blasius profile is prescribed. Its properties are looked up
   in the `case.fluid.blasius` object, see below.

### Initial conditions {#case-file_fluid-ic}
The object `initial_condition` is used to provide initial conditions.
It is mandatory.
Note that this currently pertains to both the fluid, but also scalars.
The means of prescribing the values are controlled via the `type` keyword:

1. `user`, the values are set inside the compiled user file.  as explained in the 
[user defined initial condition](@ref user-file_user-ic) section of the user
file documentation.
2. `uniform`, the value is a constant vector, looked up under the `value`
   keyword.
3. `blasius`, a Blasius profile is prescribed. Its properties are looked up
   in the `case.fluid.blasius` object, see below.
4. `point_zone`, the values are set to a constant base value, supplied under the
   `base_value` keyword, and then assigned a zone value inside a point zone. The
   point zone is specified by the `name` keyword, and should be defined in the
   `case.point_zones` object. See more about point zones @ref point-zones.md.
5. `field`, where the initial condition is retrieved from a field file.
   The following keywords can be used:
   
| Name             | Description                                                                                        | Admissible values            | Default value |
| ---------------- | -------------------------------------------------------------------------------------------------- | ---------------------------- | ------------- |
| `file_name`      | Name of the field file to use (e.g. `myfield0.f00034`).                                            | Strings ending with `f*****` | -             |
| `interpolate`    | Whether to interpolate the velocity and pressure fields from the field file onto the current mesh. | `true` or `false`            | `false`       |
| `tolerance`      | Tolerance for the point search.                                                                    | Positive real.               | `1e-6`        |
| `mesh_file_name` | If interpolation is enabled, the name of the field file that contains the mesh coordinates.        | Strings ending with `f*****` | `file_name`   |

   @attention Interpolating a field from the same mesh but different 
   polynomial order is performed implicitly and does not require to enable
   interpolation.
   
   @note It is recommended to interpolate from `fld` files that were 
   written in double precision. 
   To check if your `fld` file was written in double precision, run
   the command:
   ~~~~~~~~~~~~~~~{.sh}
   head -1 field0.f00000
   ~~~~~~~~~~~~~~~
   The output `#std 4 ...` indicates single precision, 
   whereas `#std 8 ...` indicates double precision.
   Neko write single precision `fld` files by default. To write your 
   files in double precision, set `case.output_precision` to
   `"double"`. 
  
   @attention Neko does not detect wether interpolation is needed or not. 
   Interpolation will always be performed if `"interpolate"` is set 
   to `true` even if the field file matches with the current simulation. 

### Blasius profile
The `blasius` object is used to specify the Blasius profile that can be used for the
initial and inflow condition.
The boundary cannot be tilted with respect to the coordinate axes.
It requires  the following parameters:

1. `delta`, the thickness of the boundary layer.
2. `freestream_velocity`, the velocity value in the free stream.
3. `approximation`, the numerical approximation to the Blasius profile.
   - `linear`, linear approximation.
   - `quadratic`, quadratic approximation.
   - `cubic`, cubic approximation.
   - `quartic`, quartic approximation.
   - `sin`, sine function approximation.

### Source terms {#case-file_fluid-source-term}
The `source_terms` object should be used to specify the source terms in the
momentum equation. The object is not mandatory, by default no forcing term is
present. Each source term, is itself a JSON object, so `source_terms` is just an
array of them. Note that with respect to the governing equations, the source
terms define \f$ f^u \f$, meaning that the values are then multiplied by the
density.

For each source, the `type` keyword defines the kind of forcing that will be
introduced. Furthermore, the `start_time` and `end_time` keywords can be used to
set a time frame for when the source term is active. Note, however, that these
keywords have no effect on the user-defined source terms, but their execution
can, of course, be directly controlled in the user code. By default, all source
terms are active during the entire simulation.

The following types are currently implemented.

1. `constant`, constant forcing. Strength defined by the `values` array with 3
   reals corresponding to the 3 components of the forcing.
2. `boussinesq`, a source term introducing boyancy based on the Boussinesq
   approximation, \f$ \rho \beta (T - T_{ref}) \cdot \mathbf{g} \f$. Here, \f$ \rho \f$ is
   density, \f$ \beta \f$ the thermal expansion coefficient, \f$ \mathbf{g} \f$ the
   gravity vector, and \f$ T_{ref} \f$ a reference value of the scalar, typically
   temperature.

   Reads the following entries:
   - `scalar_field`: The name of the scalar that drives the source term,
     defaults to "s".
   - `reference_value`: The reference value of the scalar.
   - `g`: The gravity vector.
   - `beta`: The thermal expansion coefficient, defaults to the inverse of
      `ref_value`.
3. `coriolis`, a source term introducing a Coriolis force, defined as \f$ -2
   \Omega \times (u - U_g) \f$. Here, \f$ \Omega \f$ is the rotation vector and
   \f$ u \f$ is the velocity vector, and \f$ U_g \f$ is the geostrophic wind.
   Several ways of setting \f$ \Omega \f$ are provided via the following
   keywords.

   - `rotation_vector`: Array with 3 values. Directly assigns \f$ \Omega \f$ to
     the provided vector.
   - `omega` and `phi`: Both scalars. The latitude `phi` should be provided in
     degrees. Sets \f$ \Omega = [0, \omega \cos \phi, \omega \sin \phi ] \f$.
     Common notation when modelling the atmosphere. This assumes that the \f$ z
     \f$ axis is normal to the ground.
   - `f`: Scalar, referred to as the Coriolis parameter, \f$ f = 2 \omega \sin
     \phi \f$. Sets \f$ \Omega = [0, 0, 0.5f ] \f$. This assumes both that \f$ z
     \f$ axis is normal to the ground and that the ground-normal component of
     the Coriolis force is negligible.

   The geostrophic wind is set to 0 for all components by default. Other values
   are set via the `geostrophic_wind` keyword.
4. `user_pointwise`, the values are set inside the compiled user file, using the
   pointwise user file subroutine. Only works on CPUs!
5. `user_vector`, the values are set inside the compiled user file, using the
   non-pointwise user file subroutine. Should be used when running on the GPU.
6. `brinkman`, Brinkman permeability forcing inside a pre-defined region.

#### Brinkman
The Brinkman source term introduces regions of resistance in the fluid domain.
The volume force \f$ f_i \f$ applied in the selected regions are proportional to the
fluid velocity component \f$ u_i \f$.

\f{eqnarray*}{
   f_i(x) &=& - B(x) u_i(x), \\
   B(x) &=& \kappa_0 + (\kappa_1 - \kappa_0) \xi(x) \frac{q + 1}{q + \xi(x)},
 \f}

where, \f$ x \f$ is the current location in the domain, \f$ \xi: x \mapsto [0,1] \f$
represent an indicator function for the resistance where \f$ \xi(x) = 0 \f$ is a free
flow. \f$ \kappa_i \f$ describes the limits for the force application at \f$ \xi(x)=0 \f$
and \f$ \xi(x)=1 \f$. A penalty parameter \f$ q \f$ help us to reduce numerical problems.

The indicator function will be defined based on the object type. The following
types are currently implemented.

1. `boundary_mesh`, the indicator function for a boundary mesh is computed in
   two steps. First, the signed distance function is computed for the boundary
   mesh. Then, the indicator function is computed using the distance transform
   function specified in the case file.
2. `point_zone`, the indicator function is defined as 1 inside the point zone
   and 0 outside.

Each object are added to a common indicator field by means of a point-wise max
operator. This means that the indicator field will be the union of all the
regions defined by the objects.

To assist correct placement and scaling of objects from external sources, the
meshes can be transformed using the `mesh_transform` object. The object can be
used to apply a transformation to the boundary mesh. The following types are
currently implemented.

1. `none`, no transformation is applied.
2. `bounding_box`, the boundary mesh is transformed to fit inside a box defined
   by `box_min` and `box_max`. The box is defined by two vectors of 3 reals
   each. The `keep_aspect_ratio` keyword can be used to keep the aspect ratio of
   the boundary mesh.

After the indicator field is computed, it is filtered using a filter type
specified in the case file. The filter is used to smooth the indicator field
before computing the Brinkman force. The following types are currently
implemented.

1. `none`, no filtering is applied.

The filtering can be defined for each object separately. Additionally, the
filter can be specified for the entire source term, in which case it will be
applied to the final indicator field, after all sources have been added.

Additional keywords are available to modify the Brinkman force term.

| Name                               | Description                                                                                   | Admissible values                 | Default value |
| ---------------------------------- | --------------------------------------------------------------------------------------------- | --------------------------------- | ------------- |
| `brinkman.limits`                  | Brinkman factor at free-flow (\f$ \kappa_0 \f$) and solid domain (\f$ \kappa_1 \f$).          | Vector of 2 reals.                | -             |
| `brinkman.penalty`                 | Penalty parameter \f$ q \f$ when estimating Brinkman factor.                                  | Real                              | \f$ 1.0 \f$   |
| `objects`                          | Array of JSON objects, defining the objects to be immersed.                                   | Each object must specify a `type` | -             |
| `distance_transform.type`          | How to map from distance field to indicator field.                                            | `step`, `smooth_step`             | -             |
| `distance_transform.value`         | Values used to define the distance transform, such as cut-off distance for the step function. | Real                              | -             |
| `filter.type`                      | Type of filtering applied to the indicator field either globally or for the current object.   | `none`                            | `none`        |
| `mesh_transform.type`              | Apply a transformation to the boundary mesh.                                                  | `bounding_box`, `none`            | `none`        |
| `mesh_transform.box_min`           | Lower left front corner of the box to fit inside.                                             | Vector of 3 reals                 | -             |
| `mesh_transform.box_max`           | Upper right back corner of the box to fit inside.                                             | Vector of 3 reals                 | -             |
| `mesh_transform.keep_aspect_ratio` | Keep the aspect ratio of the boundary mesh.                                                   | `true` or `false`                 | `true`        |

Example of a Brinkman source term where a boundary mesh and a point zone are
combined to define the resistance in the fluid domain. The indicator field for
the boundary mesh is computed using a step function with a cut-off distance of
\f$ 0.1 \f$. The indicator field for the point zone is not filtered.

~~~~~~~~~~~~~~~{.json}
"source_terms": [
   {
      "type": "brinkman",
      "objects": [
         {
            "type": "boundary_mesh",
            "name": "some_mesh.stl",
            "distance_transform": {
               "type": "step",
               "value": 0.1
            },
         },
         {
            "type": "point_zone",
            "name": "cylinder_zone",
            "filter": {
               "type": "none"
            }
         }
      ],
      "brinkman": {
         "limits": [0.0, 100.0],
         "penalty": 1.0
      }
   }
]
~~~~~~~~~~~~~~~

### Gradient Jump Penalty
The optional `gradient_jump_penalty` object can be used to perform gradient jump penalty 
as an continuous interior penalty option. 
The penalty term is performed on the weak form equation of quantity \f$ T \f$ 
(could either be velocity or scalar) as a right hand side term

\f$ - < \tau |u \cdot n| h^2_{\Omega ^e} G(T) \phi_{t1} \phi_{t2} \frac{\partial \phi_{n}}{\partial n}>\f$,

where \f$ <> \f$ refers to the integral over all facets of the element, \f$ \tau \f$ is the penalty parameter, 
\f$ |u \cdot n| \f$ is the absolute velocity flux over the facet, \f$ h^2_{\Omega ^e} \f$ is the mesh size, \f$ G(T) \f$ is the gradient jump over the facet, \f$ \phi_{t1} \phi_{t2} \f$ are the polynomial on the tangential direction of the facet, and finally \f$ \frac{\partial \phi_{n}}{\partial n} \f$ is the gradient of the normal polynomial on the facet.

Here in our Neko context where hexahedral mesh is adopted, \f$ h^2_{\Omega ^e} \f$ is measured by the average distance from the vertices of the facet to the facet on the opposite side. And the distance of a vertex to another facet is defined by the average distance from the vertex to the plane constituted by 3 vertices from the other facet.

The penalty parameter  \f$ \tau \f$ could be expressed as the form \f$ \tau = a * (P + 1) ^ {-b}\f$, 
for \f$ P > 1 \f$ where \f$ P \f$ is the polynomial order while \f$ a \f$ and \f$ b \f$ are user-defined parameters.
The configuration uses the following parameters:

* `enable`, the boolean to turn on and off the gradient jump penalty option, default to be `false`.
* `tau`, the penalty parameter that can be only used for \f$ P = 1 \f$, default to be `0.02`.
* `scaling_factor`, the scaling parameter \f$ a \f$ for \f$ P > 1 \f$, default to be `0.8`.
* `scaling_exponent`, the scaling parameter \f$ b \f$ for \f$ P > 1 \f$, default to be `4.0`.


## Linear solver configuration
The mandatory `velocity_solver` and `pressure_solver` objects are used to
configure the solvers for the momentum and pressure-Poisson equation.
The following keywords are used, with the corresponding options.

* `type`, solver type.
  - `cg`, a conjugate gradient solver.
  - `pipecg`, a pipelined conjugate gradient solver.
  - `bicgstab`, a bi-conjugate gradient stabilized solver.
  - `cacg`, a communication-avoiding conjugate gradient solver.
  - `cpldcg`, a coupled conjugate gradient solver.
  - `gmres`, a GMRES solver. Typically used for pressure.
  - `fusedcg`, a conjugate gradient solver optimised for accelerators using kernel fusion.
  - `fcpldcg`, a coupled conjugate gradient solver optimised for accelerators using kernel fusion.
* `preconditioner`, preconditioner type.
  - `jacobi`, a Jacobi preconditioner. Typically used for velocity.
  - `hsmg`, a hybrid-Schwarz multigrid preconditioner. Typically used for pressure.
  - `ident`, an identity matrix (no preconditioner).
* `absolute_tolerance`, tolerance criterion for convergence.
* `max_iterations`, maximum number of iterations before giving up.
* `projection_space_size`, size of the vector space used for accelerating the
   solution procedure. If 0, then the projection space is not used.
   More important for the pressure equation.
* `projection_hold_steps`, steps for which the simulation does not use projection after starting
   or time step changes. E.g. if 5, then the projection space will start to update at the 6th
   time step and the space will be utilized at the 7th time step.
* `monitor`, monitoring of residuals. If set to true, the residuals will be printed for each iteration.

### Flow rate forcing
The optional `flow_rate_force` object can be used to force a particular flow
rate through the domain.
Useful for channel and pipe flows.
The configuration uses the following parameters:

* `direction`, the direction of the flow, defined as 0, 1, or 2, corresponding
  to x, y or z, respectively.
* `value`, the desired flow rate.
* `use_averaged_flow`, whether `value` specifies the domain-averaged (bulk)
   velocity or the volume flow rate.


### Full parameter table
All the parameters are summarized in the table below.
This includes all the subobjects discussed above, as well as keyword parameters
that can be described concisely directly in the table.

| Name                                    | Description                                                                                       | Admissible values                                | Default value |
| --------------------------------------- | ------------------------------------------------------------------------------------------------- | ------------------------------------------------ | ------------- |
| `scheme`                                | The fluid solve type.                                                                             | `pnpn`                                           | -             |
| `Re`                                    | The Reynolds number.                                                                              | Positive real                                    | -             |
| `rho`                                   | The density of the fluid.                                                                         | Positive real                                    | -             |
| `mu`                                    | The dynamic viscosity of the fluid.                                                               | Positive real                                    | -             |
| `output_control`                        | Defines the interpretation of `output_value` to define the frequency of writing checkpoint files. | `nsamples`, `simulationtime`, `tsteps`, `never`  | -             |
| `output_value`                          | The frequency of sampling in terms of `output_control`.                                           | Positive real or integer                         | -             |
| `inflow_condition.type`                 | Velocity inflow condition type.                                                                   | `user`, `uniform`, `blasius`                     | -             |
| `inflow_condition.value`                | Value of the inflow velocity.                                                                     | Vector of 3 reals                                | -             |
| `initial_condition.type`                | Initial condition type.                                                                           | `user`, `uniform`, `blasius`, `field`            | -             |
| `initial_condition.value`               | Value of the velocity initial condition.                                                          | Vector of 3 reals                                | -             |
| `initial_condition.file_name`           | If `"type" = "field"`, the path to the field file to read from.                                   | String ending with `.fld`, `.chkp`, `.nek5000` or `f*****`.  | -             |
| `initial_condition.sample_index`        | If `"type" = "field"`, and file type is `fld` or `nek5000`, the index of the file to sampled.     | Positive integer.                                | -1            |
| `initial_condition.previous_mesh`       | If `"type" = "field"`, and file type is `chkp`, the previous mesh from which to interpolate.      | String ending with `.nmsh`.                      | -             |
| `initial_condition.tolerance`           | If `"type" = "field"`, and file type is `chkp`, tolerance to use for mesh interpolation.          | Positive real.                                   | 1e-6          |
| `blasius.delta`                         | Boundary layer thickness in the Blasius profile.                                                  | Positive real                                    | -             |
| `blasius.freestream_velocity`           | Free-stream velocity in the Blasius profile.                                                      | Vector of 3 reals                                | -             |
| `blasius.approximation`                 | Numerical approximation of the Blasius profile.                                                   | `linear`, `quadratic`, `cubic`, `quartic`, `sin` | -             |
| `source_terms`                          | Array of JSON objects, defining additional source terms.                                          | See list of source terms above                   | -             |
|`gradient_jump_penalty`                  | Array of JSON objects, defining additional gradient jump penalty.                                 | See list of gradient jump penalty above                 | -  |
| `boundary_types`                        | Boundary types/conditions labels.                                                                 | Array of strings                                 | -             |
| `velocity_solver.type`                  | Linear solver for the momentum equation.                                                          | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres`      | -             |
| `velocity_solver.preconditioner`        | Linear solver preconditioner for the momentum equation.                                           | `ident`, `hsmg`, `jacobi`                        | -             |
| `velocity_solver.absolute_tolerance`    | Linear solver convergence criterion for the momentum equation.                                    | Positive real                                    | -             |
| `velocity_solver.maxiter`               | Linear solver max iteration count for the momentum equation.                                      | Positive real                                    | 800           |
| `velocity_solver.projection_space_size` | Projection space size for the momentum equation.                                                  | Positive integer                                 | 20            |
| `velocity_solver.projection_hold_steps` | Holding steps of the projection for the momentum equation.                                        | Positive integer                                 | 5             |
| `velocity_solver.monitor`               | Monitor residuals in the linear solver for the momentum equation.                                 | `true` or `false`                                | `false`       |
| `pressure_solver.type`                  | Linear solver for the pressure equation.                                                          | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres`      | -             |
| `pressure_solver.preconditioner`        | Linear solver preconditioner for the pressure equation.                                           | `ident`, `hsmg`, `jacobi`                        | -             |
| `pressure_solver.absolute_tolerance`    | Linear solver convergence criterion for the pressure equation.                                    | Positive real                                    | -             |
| `pressure_solver.maxiter`               | Linear solver max iteration count for the pressure equation.                                      | Positive real                                    | 800           |
| `pressure_solver.projection_space_size` | Projection space size for the pressure equation.                                                  | Positive integer                                 | 20            |
| `pressure_solver.projection_hold_steps` | Holding steps of the projection for the pressure equation.                                        | Positive integer                                 | 5             |
| `pressure_solver.monitor`               | Monitor residuals in the linear solver for the pressure equation.                                 | `true` or `false`                                | `false`       |
| `flow_rate_force.direction`             | Direction of the forced flow.                                                                     | 0, 1, 2                                          | -             |
| `flow_rate_force.value`                 | Bulk velocity or volumetric flow rate.                                                            | Positive real                                    | -             |
| `flow_rate_force.use_averaged_flow`     | Whether bulk velocity or volumetric flow rate is given by the `value` parameter.                  | `true` or `false`                                | -             |
| `freeze`                                | Whether to fix the velocity field at initial conditions.                                          | `true` or `false`                                | `false`       |

## Scalar {#case-file_scalar}
The scalar object allows to add a scalar transport equation to the solution.
The solution variable is called `s`, but saved as `temperature` in the fld
 files.
Some properties of the object are inherited from `fluid`: the properties of the
linear solver, the value of the density, and the output
control.

### Material properties

The scalar equation requires defining additional material properties: the
specific heat capacity and thermal conductivity. These are provided as `cp` and
`lambda`. Similarly to the fluid, one can provide the Peclet number, `Pe`, as an
alternative. In this case, `cp` is set to 1 and `lambda` to the inverse of `Pe`.

As for the fluid, turbulence modelling is enabled by setting the `nut_field` to
the name matching that set for the simulation component with the LES model.
Additionally, the turbulent Prandtl number, `Pr_t` should be set. The eddy
viscosity values will be divided by it to produce eddy diffusivity.

### Boundary types

The boundary conditions for the scalar are specified through the
`boundary_types` keyword.

The value of the keyword is an array of strings, with the following possible values:
* Standard boundary conditions
  * `d=x`, sets a uniform Dirichlet boundary of value `x` (e.g. `d=1` to set 
  `s` to `1` on the boundary, see the Rayleigh-Benard example case).
  
* Advanced boundary conditions
    * `d_s`, a Dirichlet boundary condition for more complex, non-uniform 
    and/or time-dependent profiles. This boundary condition uses a 
    [more advanced user interface](#user-file_field-dirichlet-update).

### Initial conditions

The object `initial_condition` is used to provide initial conditions.
It is mandatory.
The means of prescribing the values are controlled via the `type` keyword:

1. `user`, the values are set inside the compiled user file as explained in the 
[user defined initial condition](@ref user-file_user-ic) section of the user
file documentation.
2. `uniform`, the value is a constant scalar, looked up under the `value`
   keyword.
3. `point_zone`, the values are set to a constant base value, supplied under the
   `base_value` keyword, and then assigned a zone value inside a point zone. The
   point zone is specified by the `name` keyword, and should be defined in the
   `case.point_zones` object. See more about point zones @ref point-zones.md.
4. `field`, where the initial condition is retrieved from a field file. Works
   in the same way as for the fluid. See the 
   [fluid section](@ref case-file_fluid-ic) for detailed explanations.
  
### Source terms

The configuration of source terms is the same as for the fluid. A demonstration
of using source terms for the scalar can be found in the `scalar_mms` example.

### Full parameter table

| Name                      | Description                                              | Admissible values               | Default value |
| ------------------------- | -------------------------------------------------------- | ------------------------------- | ------------- |
| `enabled`                 | Whether to enable the scalar computation.                | `true` or `false`               | `true`        |
| `Pe`                      | The Peclet number.                                       | Positive real                   | -             |
| `cp`                      | Specific heat cpacity.                                   | Positive real                   | -             |
| `lambda`                  | Thermal conductivity.                                    | Positive real                   | -             |
| `nut_field`               | Name of the turbulent kinematic viscosity field.         | String                          | Empty string  |
| `Pr_t`                    | Turbulent Prandtl number                                 | Positive real                   | -             |
| `boundary_types`          | Boundary types/conditions labels.                        | Array of strings                | -             |
| `initial_condition.type`  | Initial condition type.                                  | `user`, `uniform`, `point_zone` | -             |
| `initial_condition.value` | Value of the velocity initial condition.                 | Real                            | -             |
| `source_terms`            | Array of JSON objects, defining additional source terms. | See list of source terms above  | -             |
|`gradient_jump_penalty`    | Array of JSON objects, defining additional gradient jump penalty. | See list of gradient jump penalty above | -  |

## Statistics

This object adds the collection of statistics for the fluid fields. For
additional details on the workflow, see the
[corresponding page](@ref statistics-guide) in the user manual.

| Name                | Description                                                          | Admissible values | Default value |
| ------------------- | -------------------------------------------------------------------- | ----------------- | ------------- |
| `enabled`           | Whether to enable the statistics computation.                        | `true` or `false` | `true`        |
| `start_time`        | Time at which to start gathering statistics.                         | Positive real     | 0             |
| `sampling_interval` | Interval, in timesteps, for sampling the flow fields for statistics. | Positive integer  | 10            |

## Simulation components
Simulation components enable the user to perform various additional operations,
which are not strictly necessary to run the solver. An example could be
computing and output of additional fields, e.g. vorticity.

A more detailed description as well as a  full list of available components and
 their setup is provided in a [separate page of the manual](simcomps.md).

## Point zones
Point zones enable the user to select GLL points in the computational domain
according to some geometric criterion. Two predefined geometric shapes are
selectable from the case file, boxes and spheres.

A point zone object defined in the case file can be retrieved from the point
zone registry, `neko_point_zone_registry`, and can be used to perform any
zone-specific operations (e.g. localized source term, probing...). User-specific
point zones can also be added manually to the point zone registry from the user
file.

A more detailed description as well as a  full list of available components and
 their setup is provided in a [separate page of the manual](point-zones.md).

## Runtime statistics

This object adds the collection of runtime statistics (timings) for identified
profiling regions. A region is defined as all functions between a call to
`profiler_start_region(name, id)` and `profiler_end_region(name, id)`. Neko
currently supports 50 regions, with id 1..25 being reserved for internal use.


| Name                | Description                                                          | Admissible values | Default value |
| ------------------- | -------------------------------------------------------------------- | ----------------- | ------------- |
| `enabled`           | Whether to enable gathering of runtime statistics                    | `true` or `false` | `false`       |
| `output_profile`    | Wheter to output all gathered profiling data as a CSV file           | `true` or `false` | `false`       |
