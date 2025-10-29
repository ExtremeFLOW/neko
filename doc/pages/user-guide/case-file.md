# Case File {#case-file}

\tableofcontents

The case file defines all the parameters of a simulation.
The format of the file is JSON, making it easy to read and write case files
using the majority of the popular programming languages.
JSON is hierarchical and, and consists of parameter blocks enclosed in curly
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
        "time": {}
        "numerics": {}
        "fluid": {}
        "scalar": {}
        "simulation_components" : []
        "point_zones" : []
    }
}
~~~~~~~~~~~~~~~
Neko also supports multiple scalar fields, using the keyword `scalars` in the
case object. Users can define multiple scalar fields, and each field can have
its own boundary conditions, source terms, and solver settings. When using
multiple scalar fields, the `name` property of each scalar field is used to
identify the scalar field in the user file, defaulted to `s_1, s_2, ...`.

The `version` keyword is reserved to track changes in the format of the file.
The subsections below we list all the configuration options for each of the high-level objects.
Some parameters will have default values, and are therefore optional.

## Output frequency control
A common scheme for controlling the output frequency is applied for various
outputs.
It is described already now in order to clarify the meaning of several
parameters found in the tables below.

The frequency is controlled by two parameters, ending with `_control` and
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

| Name                  | Description                                                                                           | Admissible values                               | Default value |
| --------------------- | ----------------------------------------------------------------------------------------------------- | ----------------------------------------------- | ------------- |
| `no_defaults`         | Prevents filling in default values to case file entries.                                              | `true` or `false`                               | `false`       |
| `mesh_file`           | The name of the mesh file.                                                                            | Strings ending with `.nmsh`                     | -             |
| `output_boundary`     | Whether to write a `bdry0.f0000` file with boundary labels. Can be used to check boundary conditions. | `true` or `false`                               | `false`       |
| `output_directory`    | Folder for redirecting solver output. Note that the folder has to exist!                              | Path to an existing directory                   | `.`           |
| `output_format`       | The file format of field data.                                                                        | `nek5000` or `adios2`                           | `nek5000`     |
| `output_precision`    | Whether to output snapshots in single or double precision                                             | `single` or `double`                            | `single`      |
| `output_layout`       | Data layout for `adios2` files. (Choose `2` or `3` for ADIOS2 supported compressors BigWhoop or ZFP.) | Positive integer `1`, `2`, `3`                  | `1`           |
| `load_balancing`      | Whether to apply load balancing.                                                                      | `true` or `false`                               | `false`       |
| `output_partitions`   | Whether to write a `partitions.vtk` file with domain partitioning.                                    | `true` or `false`                               | `false`       |
| `output_checkpoints`  | Whether to output checkpoints, i.e. restart files.                                                    | `true` or `false`                               | `false`       |
| `checkpoint_control`  | Defines the interpretation of `checkpoint_value` to define the frequency of writing checkpoint files. | `nsamples`, `simulationtime`, `tsteps`, `never` | -             |
| `checkpoint_value`    | The frequency of sampling in terms of `checkpoint_control`.                                           | Positive real or integer                        | -             |
| `checkpoint_filename` | The filename of written checkpoint.                                                                   | Strings such as `my_name`                       | `fluid`       |
| `checkpoint_format`   | The file format of checkpoints                                                                        | `chkp` or `hdf5`                                | `chkp`        |
| `restart_file`        | checkpoint to use for a restart from previous data                                                    | Strings ending with `.chkp`                     | -             |
| `restart_mesh_file`   | If the restart file is on a different mesh, specify the .nmsh file used to generate it here           | Strings ending with `.nmsh`                     | -             |
| `mesh2mesh_tolerance` | Tolerance for the restart when restarting from another mesh                                           | Positive reals                                  | 1e-6          |
| `job_timelimit`       | The maximum wall clock duration of the simulation.                                                    | String formatted as HH:MM:SS                    | No limit      |
| `output_at_end`       | Whether to always write all enabled output at the end of the run.                                     | `true` or `false`                               | `true`        |


### Time control
The time control object is used to define the time-stepping of the simulation,
including the time-step size, the start and end time, and the variables related
to the variable time-stepping algorithm.

| Name                       | Description                                                                                 | Admissible values                 | Default value |
| -------------------------- | ------------------------------------------------------------------------------------------- | --------------------------------- | ------------- |
| `start_time`               | Start time at which the simulation is initiated.                                            | Positive reals                    | `0.0`         |
| `end_time`                 | Final time after which the simulation is stopped.                                           | Positive reals                    | -             |
| `timestep`                 | Time-step size                                                                              | Positive reals                    | -             |
| `variable_timestep`        | Whether to use variable dt                                                                  | `true` or `false`                 | `false`       |
| `max_timestep`             | Maximum time-step size when variable time step is activated                                 | Positive reals                    | `huge`        |
| `min_timestep`             | Minimum time-step size when variable time step is activated                                 | Positive reals                    | `0.0`         |
| `target_cfl`               | The desired CFL number                                                                      | Positive real                     | `0.4`         |
| `max_update_frequency`     | The minimum interval between two time-step-updating steps in terms of time steps            | Integer                           | `0`           |
| `min_update_frequency`     | The maximum interval between two time-step-updating steps in terms of time steps            | Integer                           | `huge`        |
| `running_avg_coeff`        | The running average coefficient `a` where `cfl_avg_new = a * cfl_new + (1-a) * cfl_avg_old` | Positive real between `0` and `1` | `0.5`         |
| `max_dt_increase_factor`   | The maximum scaling factor to increase time step                                            | Positive real greater than `1`    | `1.2`         |
| `min_dt_decrease_factor`   | The minimum scaling factor to decrease time step                                            | Positive real less than `1`       | `0.5`         |
| `cfl_deviation_tolerance`  | The tolerance of the deviation from the target CFL number                                   | Positive real less than `1`       | `0.2`         |
| `cfl_max_update_frequency` | The minimum interval between two time-step-updating steps in terms of time steps            | Integer                           | `0`           |
| `cfl_running_avg_coeff`    | The running average coefficient `a` where `cfl_avg_new = a * cfl_new + (1-a) * cfl_avg_old` | Positive real between `0` and `1` | `0.5`         |

### Restarts and joblimit
Restarts will restart the simulation from the exact state at a given time that
the checkpoint was written. This means that the flow field and potential scalars
will be at the exact same values before as after restarts. However, derived
quantities from the flow field and any observables are not guaranteed to be
restarted. In addition, Neko does not guarantee that any files are not
overwritten. As such, it is recommended to run in different directories
if doing large scale simulations that require many restarts. Unless
`output_at_end` is disabled Neko will also ensure that all output is written to
file when reaching the `end_time` or the `job_timelimit`. In particular, unless
`output_checkpoints` and `output_at_end` are set to false a checkpoint at the
final time will be written as to avoid losing progress as far as possible.

@attention For simulations requiring restarts, it is recommended to run each
restart in a different output directory as a precaution to avoid potential overwritings of files.

### Boundary type numbering in the "output_boundary" field

When the `output_boundary` setting is set to `true`, and additional `.fld` file
will be stored in the beginning of the simulation, where the recognized boundary
conditions for the fluid  will be marked with an integer number. This is a good
way to debug the simulation setup. The value of the number depends on the type
of the boundary as follows.

| Boundary Condition              | Key |
| ------------------------------- | --- |
| no_slip                         | 1   |
| velocity_value                  | 2   |
| outflow, normal_outflow (+dong) | 3   |
| symmetry                        | 4   |
| user_velocity_pointwise         | 5   |
| periodic                        | 6   |
| user_velocity                   | 7   |
| user_pressure                   | 8   |
| shear_stress                    | 9   |
| wall_model                      | 10  |
| blasius_profile                 | 11  |

For a description of the boundary conditions themselves, see below.

## Numerics
Used to define the properties of the numerical discretization.

| Name                         | Description                                                                                                     | Admissible values          | Default value                   |
| ---------------------------- | --------------------------------------------------------------------------------------------------------------- | -------------------------- | ------------------------------- |
| `polynomial_order`           | The order of the polynomial basis.                                                                              | Integers, typically 5 to 9 | -                               |
| `time_order`                 | The order of the time integration scheme. Refer to the `time_scheme_controller` type documentation for details. | 1, 2, 3                    | -                               |
| `dealias`                    | Whether to apply dealiasing to advection terms.                                                                 | `true` or `false`          | `false`                         |
| `dealiased_polynomial order` | The polynomial order in the higher-order space used in the dealising.                                           | Integer                    | `3/2(polynomial_order + 1) - 1` |
| `oifs`                       | Whether to apply the Operator-Integration-Factor-Splitting (OIFS).                                              | `true` or `false`          | `false`                         |
| `oifs_target_cfl`            | The desired OIFS-CFL number. Requires variable_timestep = true in the time control object.                      | Positive real              | `1.9`                           |

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
in the `rayleigh_benard_cylinder` example. Ultimately, both `rho` and `mu` have
to be set in the subroutine.  Additionally, this allows to change the material
properties in time. Yet another options is to directly manipulate the case file
programmatically in the `user_startup` routine and inject the material
properties there. This is demonstrated in the `rayleigh_benard` example.

When material properties are constant or only vary in time, one can use the
simplified form of the viscous stress tensor in the governing equations.
However, when there are spatial variations, it is necessary to use the general
(full) form. The variation may come, for example,  due to a turbulence model,
the modifications in the above-mentioned user routine. The general form of the
stress tensor requires solving the 3 equations for the velocity components in a
coupled manner, which requires an appropriate linear solver. By default, Neko
will use the simplified form of the tensor, and the full one must be selected
by the user by setting `full_stress_formulation` to true.

### Turbulence modelling

Neko currently provides several LES models via the `les_model` simulation
component. The simcomp computes the viscosity and stores in the field registry
under the name selected in the case file. For more details, see the
documentation of the simcomp. To enable LES in the fluid solver, one simply has
to add the `nut_field` keyword to the configuration and set it to the name of
the field generated by the `les_model` simcomp. This will automatically change
the governing equations to feature the full viscous stress tensor, as required
for a variable viscosity field.

Note that the full viscous stress tensor requires the equations for the 3
velocity components to be solved in a coupled manner. Therefore, the `cpldcg`
solver should be used for velocity.

### Boundary conditions {#case-file_fluid-boundary-conditions}
The optional `boundary_conditions` keyword can be used to specify boundary
conditions. The reason for it being optional, is that periodic boundary
conditions are built into the definition of the  mesh, so for a periodic box
nothings needs to be added to the case file. The TGV example is such a case, for
instance.
The value of the keyword is an array of JSON objects, each specifying a single
boundary condition.

#### Specifying the boundaries

In Neko we usually refer to boundaries as "zones", which in this case are face
zones, i.e. a collection of element faces. Which zones the boundary condition is
applied to is controlled by the `zone_indices` keyword, which takes an array of
integers. It is up to the user whether to apply a single condition to multiple
zones or specify several conditions applied to one zone each. For example, if
you have two zones, which should be no-slip walls, you can either create two
`no_slip` conditions, one for each zone, or just create a single condition and
apply it to both.

The indices your boundaries have is determined by the mesh. To check them, you
can use the `mesh_checker` utility with the optional `--write_zone_indices`
argument. This will output a `zone_indices0.f00000` file that you can inspect in
Paraview, and the boundaries will be marked by their index value.

Recall that periodic conditions are built into the mesh, since they are
topological in nature. This means that you must not specify any conditions for
the corresponding zones. For example, in the `turb_pipe` example, which is a
periodic pipe simulation, two periodic zones comprise the boundary conditions in
the streamwise direction. Only one condition, corresponding to zone index 3 (the
wall) is the specified in the case file.

#### Available conditions
The conditions to apply is specified by `type` keyword inside each of the JSON
objects. The full list of possible conditions for the fluid is specified in the
table below.

| Boundary Condition      | Description                                                                                                                                            |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| symmetry                | A symmetry plain. Must be axis-aligned.                                                                                                                |
| velocity_value          | A Dirichlet condition for velocity.                                                                                                                    |
| no_slip                 | A no-slip wall.                                                                                                                                        |
| outflow                 | A pressure outlet.                                                                                                                                     |
| normal_outflow          | An Neumann condition for the surface-normal component of velocity combined with a Dirichlet for the surface-parallel components. Must be axis-aligned. |
| outflow+user            | Same as `outflow` but with user-specified pressure.                                                                                                    |
| normal_outflow+user     | Same as `normal_outflow` but with user-specified pressure.                                                                                             |
| outflow+dong            | A pressure outlet with the Dong condition applied.                                                                                                     |
| normal_outflow+dong     | The `normal_outflow` with the Dong condition applied. Must be axis-aligned.                                                                            |
| shear_stress            | Prescribed wall shear stress. Must be axis-aligned.                                                                                                    |
| wall_model              | Shear stress condition based on a wall model for large-eddy simulation.                                                                                |
| blasius_profile         | A Blasius velocity profile.                                                                                                                            |
| user_velocity           | The `field_dirichlet_vector_t` user-defined Dirichlet condition for velocity.                                                                          |
| user_pressure           | The `field_dirichlet_t` user-defined Dirichlet condition for pressure.                                                                                 |
| user_velocity_pointwise | The pointwise user-defined Dirichlet condition for velocity.                                                                                           |

A more detailed description of each boundary condition is provided below.

* `symmetry`. A symmetry plain that must be axis-aligned. Sets the
  surface-normal velocity to 0 and applies a homogenous Neumann condition to the
  surface-parallel components. Requires no additional keywords.
  ```json
  {
    "type": "symmetry",
    "zone_indices": [1, 2]
  }
  ```
* `velocity_value`. A Dirichlet condition for velocity. Suitable for velocity
  inlets, moving walls, certain freestream conditions, etc. The value is
  prescribed by the `value` keyword, that should be an array of 3 reals.

  ```json
  {
    "type": "velocity_value",
    "value": ["1, 0, 0"],
    "zone_indices": [1, 2]
  }
  ```
* `no_slip`. A standard no-slip wall, which sets velocity to zero. Requires no
  additional keywords.
  ```json
  {
    "type": "no_slip",
    "zone_indices": [1, 2]
  }
  ```
* `outflow`. A standard pressure outlet condition. Requires no additional
  keywords.
  ```json
  {
    "type": "outflow",
    "zone_indices": [1, 2]
  }
  ```
* `normal_outflow`. The condition lets the flow escape through the boundary by
  setting a homogeneous Neumann condition for the surface-normal velocity
  component, but fixes the values of the surface-parallel components. The latter
  values are not prescribed in the boundary condition's JSON, but are instead
  taken from the initial conditions. The boundary must be axis-aligned.
  ```json
  {
    "type": "normal_outflow",
    "zone_indices": [1, 2]
  }
  ```
* `outflow+user`. Same as `outflow`, but with user-specified
  pressure. The pressure is specified via the same interface as `user_pressure`,
  see the
  [relevant section](#user-file_field-dirichlet-update) for more information.

* `normal_outflow+user`. Same as `normal_outflow`, but with user-specified
  pressure. The pressure profile is specified via the same interface as
  `user_pressure`, see
  the [relevant section](#user-file_field-dirichlet-update) for more information.
  Note that, similarly to `normal_outflow`, surface-parallel velocity components
  are taken from the initial conditions.

* `outflow+dong`. Same as `outflow`, but additionally applies the Dong boundary
  condition on the pressure. This is a way to prevent backflow and therefore
  promotes numerical stability without having to create a mesh "sponge" at the
  outlet.

* `normal_outflow+dong`. Same as `normal_outflow`, but additionally applies the
  Dong boundary condition for the pressure to prevent backflow. Must be
  axis-aligned.

* `shear_stress`. Non-penetration condition combined with a set shear stress
   vector. Only works with axis-aligned boundaries. The stress value is
   specified by the `value` keyword, which should be an array of 3 reals. It is
   the responsibility of the user to set the vector in the direction parallel
   to the boundary.
  ```json
  {
    "type": "shear_stress",
    "value": ["1, 0, 0"],
    "zone_indices": [1, 2]
  }
  ```

* `wall_model`. A shear stress condition, where the values is computed by a wall
   model. Meant to be used for wall-modelled large-eddy simulation. Only works
   with axis-aligned boundaries. The model is selected using the `model`
   keyword. Additional configuration depends on the model selected.

   * The `spalding`model requires specifying `kappa` and `B`, which are the
     log-law constants. This model is suitable for smooth walls.

   * The `rough_log_law` model requires specifying `kappa` and `B`, which are
     the log-law constants, and `z0`, which is the characteristic roughness
     height.

    For all wall models, the distance to the sampling point has to be specified
    based on the off-wall index in the wall-normal direction. Thus, the sampling
    is currently from a GLL node and arbitrary distances are not yet supported.
    The index is set by the `h_index` keyword, with 1 being the minimal value, and
    the polynomial order + 1 being the maximum.

    A 3D field with the name `tau` will be registered in the field registry. At
    the boundary it will store the magnitude of the predicted stress. This can
    be used to post-process the predictions. Additionally, the sampling points
    are marked with values -1 in this field, for verification purposes.

  ```json
  {
    "type": "wall_model",
    "model": "spalding",
    "kappa": 0.41,
    "B": 5.2,
    "zone_indices": [1, 2],
    "h_index": 1
  }
  ```
* `user_pointwise`. Allows to set the velocity values using the appropriate
  routine in the user file. The routine is executed on a pointwise basis, which
  is reflected in the name of this condition. It is advisable to instead use the
  more general `user_velocity` condition. Requires no additional keywords.

  ```json
  {
    "type": "user_pointwise",
    "zone_indices": [1, 2]
  }
  ```
* `user_velocity`, a Dirichlet boundary for more complex velocity profiles. This boundary
  condition uses a [more advanced user
  interface](#user-file_field-dirichlet-update).
  ```json
  {
    "type": "user_velocity",
    "zone_indices": [1, 2]
  }
  ```
* `user_pressure`, a boundary for specified non-uniform pressure profiles, similar in
  essence to `user_velocity`.
  ```json
  {
    "type": "user_pressure",
    "zone_indices": [1, 2]
  }
  ```

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
3. `blasius`, a Blasius profile is prescribed. The boundary cannot be tilted
  with respect to the coordinate axes.
   It requires the following parameters:
   1. `delta`, the thickness of the boundary layer.
   2. `freestream_velocity`, the velocity value in the free stream.
   3. `approximation`, the numerical approximation to the Blasius profile.
      - `linear`, linear approximation.
      - `quadratic`, quadratic approximation.
      - `cubic`, cubic approximation.
      - `quartic`, quartic approximation.
      - `sin`, sine function approximation.
      - `tanh`, hyperbolic tangent approximation of Savaş (2012). In this case `delta` is the 99\% thickness.
4. `point_zone`, the values are set to a constant base value, supplied under the
   `base_value` keyword, and then assigned a zone value inside a point zone. The
   point zone is specified by the `name` keyword, and should be defined in the
   `case.point_zones` object. See more about [point zones](@ref point-zones).
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
2. `boussinesq`, a source term introducing buoyancy based on the Boussinesq
   approximation, \f$ \rho \beta (T - T_{ref}) \cdot \mathbf{g} \f$. Here, \f$
   \rho \f$ is density, \f$ \beta \f$ the thermal expansion coefficient, \f$
   \mathbf{g} \f$ the gravity vector, and \f$ T_{ref} \f$ a reference value of
   the scalar, typically temperature.

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
4. `centrifugal`, a source term introducing a centrifugal force, defined as \f$ -
   \Omega \times (\Omega \times r) \f$. Here, \f$ \Omega \f$ is the rotation
   vector and \f$ r \f$ is the position relative to the reference point, which
   is any point lying on the rotation axis. To define forcing one has to provide
   \f$ \Omega \f$ and the reference point. This is provided via the following
   keywords.

   - `rotation_vector`: Array with 3 values. Directly assigns \f$ \Omega \f$ to
     the provided vector.
   - `reference_point`: Array with 3 values. Deifines any point on the rotaion
   axis.

5. `user_pointwise`, the values are set inside the compiled user file, using the
   pointwise user file subroutine. Only works on CPUs!
6. `user_vector`, the values are set inside the compiled user file, using the
   non-pointwise user file subroutine. Should be used when running on the GPU.
7. `brinkman`, Brinkman permeability forcing inside a pre-defined region.
8. `gradient_jump_penalty`, perform gradient_jump_penalisation.

@note Notice that to perform simulation in a rotating reference frame one has to
define both `coriolis` and `centrifugal` source terms in a consistent way.

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
   function specified in the case file. This is currently not very well
   optimized, it will scale by `O(log(M)*N)`, where `M` is the number of triangles in the
   boundary mesh and `N` is the number of grid points in the simulation.
   To avoid recomputing the distance field for multiple simulations with the
   same boundary mesh and numerical discretization, the distance field can be
   cached to a file. This is controlled by the `cache` keyword. If set to
   `true`, the distance field will be saved to a file specified by the
   `cache_file` keyword. If the file already exists, it will be loaded instead
   of recomputing it. The distance field is stored in the Nek5000 `.fld` file
   format.
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
            "cache": true,
            "cache_file": "some_mesh_cache"
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

#### Gradient Jump Penalty
The optional `gradient_jump_penalty` object can be used to perform gradient jump
penalty as an continuous interior penalty option. The penalty term is performed
on the weak form equation of quantity \f$ T \f$ (could either be velocity or
scalar) as a right hand side term

\f$ - < \tau |u \cdot n| h^2_{\Omega ^e} G(T) \phi_{t1} \phi_{t2} \frac{\partial \phi_{n}}{\partial n}>\f$,

where \f$ <> \f$ refers to the integral over all facets of the element, \f$ \tau
\f$ is the penalty parameter, \f$ |u \cdot n| \f$ is the absolute velocity flux
over the facet, \f$ h^2_{\Omega ^e} \f$ is the mesh size, \f$ G(T) \f$ is the
gradient jump over the facet, \f$ \phi_{t1} \phi_{t2} \f$ are the polynomial on
the tangential direction of the facet, and finally \f$ \frac{\partial
\phi_{n}}{\partial n} \f$ is the gradient of the normal polynomial on the facet.

Here in our Neko context where hexahedral mesh is adopted, \f$ h^2_{\Omega ^e}
\f$ is measured by the average distance from the vertices of the facet to the
facet on the opposite side. And the distance of a vertex to another facet is
defined by the average distance from the vertex to the plane constituted by 3
vertices from the other facet.

The penalty parameter  \f$ \tau \f$ could be expressed as the form \f$ \tau = a
* (P + 1) ^ {-b}\f$, for \f$ P > 1 \f$ where \f$ P \f$ is the polynomial order
while \f$ a \f$ and \f$ b \f$ are user-defined parameters. The configuration
uses the following parameters:

* `tau`, the penalty parameter that can be only used for \f$ P = 1 \f$, default
  to be `0.02`.
* `scaling_factor`, the scaling parameter \f$ a \f$ for \f$ P > 1 \f$, default
  to be `0.8`.
* `scaling_exponent`, the scaling parameter \f$ b \f$ for \f$ P > 1 \f$, default
  to be `4.0`.


## Linear solver configuration
The mandatory `velocity_solver` and `pressure_solver` objects are used to
configure the solvers for the momentum and pressure-Poisson equation.
The following keywords are used, with the corresponding options.

* `type`, solver type.
  - `cg`, a conjugate gradient solver.
  - `pipecg`, a pipelined conjugate gradient solver.
  - `bicgstab`, a bi-conjugate gradient stabilized solver.
  - `cacg`, a communication-avoiding conjugate gradient solver.
  - `coupled_cg`, a coupled conjugate gradient solver. Must be used for velocity
    when viscosity varies in space.
  - `gmres`, a GMRES solver. Typically used for pressure.
  - `fused_cg`, a conjugate gradient solver optimised for accelerators using
  - `fused_coupled_cg`, a coupled conjugate gradient solver optimised for accelerators using
    kernel fusion. Must be used for velocity when viscosity varies in space and
    device backened is used.
    using kernel fusion.
* `preconditioner.type`, preconditioner type.
  - `jacobi`, a Jacobi preconditioner. Typically used for velocity.
  - `hsmg`, a hybrid-Schwarz multigrid preconditioner. Typically used for
    pressure.
  - `phmg`, a hybrid ph multigrid preconditioner. Typically used for pressure.
  - `ident`, an identity matrix (no preconditioner).
* `absolute_tolerance`, tolerance criterion for convergence.
* `max_iterations`, maximum number of iterations before giving up.
* `projection_space_size`, size of the vector space used for accelerating the
   solution procedure. If 0, then the projection space is not used.
   More important for the pressure equation.
* `projection_hold_steps`, steps for which the simulation does not use
   projection after starting or time step changes. E.g. if 5, then the
   projection space will start to update at the 6th time step and the space will
   be utilized at the 7th time step.
* `monitor`, monitoring of residuals. If set to true, the residuals will be
  printed for each iteration.

In addition to the above settings, the solvers can be configured with strict
convergence criteria. This is done by setting the
`case.fluid.strict_convergence` keyword to `true`. This will force the solver to
converge to the specified tolerance within the specified number of iterations.
If the solver does not converge, the simulation will be terminated.

### Multilevel preconditioners
The multilevel preconditioners, `hsmg` and `phmg`, come with an
additional set of parameters related to the solution of the coarse
grid problem. These parameters are specified as a `coarse_grid` block
under `preconditioner`.

For `hsmg`, the following keywords are used:

| Name                         | Description                                                                             | Admissible values       | Default value |
| ---------------------------- | --------------------------------------------------------------------------------------- | ----------------------- | ------------- |
| `coarse_grid.solver`         | Type of linear solver for the coarse grid, any of the Krylov solvers or TreeAMG  `tamg` | A solver `type`         | `cg`          |
| `coarse_grid.preconditioner` | Type of the preconditioner to use (only valid for a Krylov based `solver`)              | A preconditioner `type` | `jacobi`      |
| `coarse_grid.iterations`     | Number of linear solver iterations (only valid for a Krylov based `solver`)             | An integer              | 10            |
| `coarse_grid.monitor`        | Monitor residuals in the coarse grid (only valid for a Krylov based `solver`)           | `true` or `false`       | `false`       |
| `coarse_grid.levels`         | Number of AMG levels to construct (only valid for `solver` type `tamg`)                 | An integer              | 3             |
| `coarse_grid.iterations`     | Number of AMG iterations (only valid for `solver` type `tamg`)                          | An integer              | 1             |
| `coarse_grid.cheby_degree`   | Degree of the Chebyshev based AMG smoother                                              | An integer              | 5             |

For `phmg`, the following keywords are used:

| Name                       | Description                                                                                 | Admissible values     | Default value |
| -------------------------- | ------------------------------------------------------------------------------------------- | --------------------- | ------------- |
| `pcoarsening_schedule`     | P-multigrid coarsening schedule (polynomial order, high to low)                             | Array of integers     | `[3, 1]`      |
| `smoother_iterations`      | Number of smoother iterations in the p-multigrid parts                                      | An integer            | 10            |
| `smoother_cheby_acc`       | Type of Chebyshev acceleration (non-accelerated semi-iterative Chebyshev method if not set) | `jacobi` or `schwarz` | -             |
| `coarse_grid.levels`       | Number of AMG levels to construct (only valid for `solver` type `tamg`)                     | An integer            | 3             |
| `coarse_grid.iterations`   | Number of linear solver iterations for coarse grid solver                                   | An integer            | 1             |
| `coarse_grid.cheby_degree` | Degree of the Chebyshev based AMG smoother                                                  | An integer            | 5             |


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
All the parameters are summarized in the table below. This includes all the
subobjects discussed above, as well as keyword parameters that can be described
concisely directly in the table.

| Name                                    | Description                                                                                       | Admissible values                                           | Default value |
| --------------------------------------- | ------------------------------------------------------------------------------------------------- | ----------------------------------------------------------- | ------------- |
| `scheme`                                | The fluid solve type.                                                                             | `pnpn`                                                      | -             |
| `name`                                  | The name associated to the fluid solver.                                                          | String                                                      | `fluid`       |
| `Re`                                    | The Reynolds number.                                                                              | Positive real                                               | -             |
| `rho`                                   | The density of the fluid.                                                                         | Positive real                                               | -             |
| `mu`                                    | The dynamic viscosity of the fluid.                                                               | Positive real                                               | -             |
| `nut_field`                             | The name of the turbulent viscosity field.                                                        | String                                                      | -             |
| `output_control`                        | Defines the interpretation of `output_value` to define the frequency of writing checkpoint files. | `nsamples`, `simulationtime`, `tsteps`, `never`             | -             |
| `output_value`                          | The frequency of sampling in terms of `output_control`.                                           | Positive real or integer                                    | -             |
| `output_filename`                       | The output filename.                                                                              | String                                                      | `field`       |
| `inflow_condition.type`                 | Velocity inflow condition type.                                                                   | `user`, `uniform`, `blasius`                                | -             |
| `inflow_condition.value`                | Value of the inflow velocity.                                                                     | Vector of 3 reals                                           | -             |
| `initial_condition.type`                | Initial condition type.                                                                           | `user`, `uniform`, `blasius`, `field`                       | -             |
| `initial_condition.value`               | Value of the velocity initial condition.                                                          | Vector of 3 reals                                           | -             |
| `initial_condition.file_name`           | If `"type" = "field"`, the path to the field file to read from.                                   | String ending with `.fld`, `.chkp`, `.nek5000` or `f*****`. | -             |
| `initial_condition.sample_index`        | If `"type" = "field"`, and file type is `fld` or `nek5000`, the index of the file to sampled.     | Positive integer.                                           | -1            |
| `initial_condition.previous_mesh`       | If `"type" = "field"`, and file type is `chkp`, the previous mesh from which to interpolate.      | String ending with `.nmsh`.                                 | -             |
| `initial_condition.tolerance`           | If `"type" = "field"`, and file type is `chkp`, tolerance to use for mesh interpolation.          | Positive real.                                              | 1e-6          |
| `blasius.delta`                         | Boundary layer thickness in the Blasius profile.                                                  | Positive real                                               | -             |
| `blasius.freestream_velocity`           | Free-stream velocity in the Blasius profile.                                                      | Vector of 3 reals                                           | -             |
| `blasius.approximation`                 | Numerical approximation of the Blasius profile.                                                   | `linear`, `quadratic`, `cubic`, `quartic`, `sin`, `tanh`    | -             |
| `shear_stress.value`                    | The shear stress vector value for `sh` boundaries                                                 | Vector of 3 reals                                           | `[0, 0, 0]`   |
| `wall_modelling.type`                   | The wall model type for `wm` boundaries. See documentation for additional config parameters.      | `rough_log_law`, `spalding`                                 | -             |
| `source_terms`                          | Array of JSON objects, defining additional source terms.                                          | See list of source terms above                              | -             |
| `gradient_jump_penalty`                 | Array of JSON objects, defining additional gradient jump penalty.                                 | See list of gradient jump penalty above                     | -             |
| `boundary_types`                        | Boundary types/conditions labels.                                                                 | Array of strings                                            | -             |
| `velocity_solver.type`                  | Linear solver for the momentum equation.                                                          | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres`                 | -             |
| `velocity_solver.preconditioner.type`   | Linear solver preconditioner for the momentum equation.                                           | `ident`, `hsmg`, `jacobi`                                   | -             |
| `velocity_solver.absolute_tolerance`    | Linear solver convergence criterion for the momentum equation.                                    | Positive real                                               | -             |
| `velocity_solver.maxiter`               | Linear solver max iteration count for the momentum equation.                                      | Positive real                                               | 800           |
| `velocity_solver.projection_space_size` | Projection space size for the momentum equation.                                                  | Positive integer                                            | 0             |
| `velocity_solver.projection_hold_steps` | Holding steps of the projection for the momentum equation.                                        | Positive integer                                            | 5             |
| `velocity_solver.monitor`               | Monitor residuals in the linear solver for the momentum equation.                                 | `true` or `false`                                           | `false`       |
| `pressure_solver.type`                  | Linear solver for the pressure equation.                                                          | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres`                 | -             |
| `pressure_solver.preconditioner.type`   | Linear solver preconditioner for the pressure equation.                                           | `ident`, `hsmg`, `jacobi`                                   | -             |
| `pressure_solver.absolute_tolerance`    | Linear solver convergence criterion for the pressure equation.                                    | Positive real                                               | -             |
| `pressure_solver.maxiter`               | Linear solver max iteration count for the pressure equation.                                      | Positive real                                               | 800           |
| `pressure_solver.projection_space_size` | Projection space size for the pressure equation.                                                  | Positive integer                                            | 0             |
| `pressure_solver.projection_hold_steps` | Holding steps of the projection for the pressure equation.                                        | Positive integer                                            | 5             |
| `pressure_solver.monitor`               | Monitor residuals in the linear solver for the pressure equation.                                 | `true` or `false`                                           | `false`       |
| `flow_rate_force.direction`             | Direction of the forced flow.                                                                     | 0, 1, 2                                                     | -             |
| `flow_rate_force.value`                 | Bulk velocity or volumetric flow rate.                                                            | Positive real                                               | -             |
| `flow_rate_force.use_averaged_flow`     | Whether bulk velocity or volumetric flow rate is given by the `value` parameter.                  | `true` or `false`                                           | -             |
| `freeze`                                | Whether to fix the velocity field at initial conditions.                                          | `true` or `false`                                           | `false`       |
| `advection`                             | Whether to compute the advection term.                                                            | `true` or `false`                                           | `true`        |
| `full_stress_formulation`               | Whether to use the full form of the visous stress tensor term.                                    | `true` or `false`                                           | `false`       |

## Scalar {#case-file_scalar}
The scalar object allows to add a scalar transport equation to the solution. The
solution variable is called `s` by default, but can be controlled by the
 `field_name` entry in the case file. In the fld files, it is saved as
`temperature`. Some properties of the object are inherited from `fluid`: the
value of the density, and the output control.

### Material properties

The scalar equation requires defining additional material properties: the
specific heat capacity and thermal conductivity. These are provided as `cp` and
`lambda`. Similarly to the fluid, one can provide the Peclet number, `Pe`, as an
alternative. In this case, `cp` is set to 1 and `lambda` to the inverse of `Pe`.

As for the fluid, turbulence modelling is enabled by setting the `nut_field` to
the name matching that set for the simulation component with the LES model.
Additionally, the turbulent Prandtl number, `Pr_t` should be set. The eddy
viscosity values will be divided by it to produce eddy diffusivity.

### Turbulence modelling

The configuration is identical to the Fluid, however, one additionally has to
provide the value of the turbulent Prandl number via the `Pr_t` keyword.

### Boundary conditions

The boundary conditions for the scalar are specified through the
`boundary_conditions` keyword, which follows the same format as the fluid, for
specifying the type of the condition and where it is applied.
Four types of conditions are available for the scalar:

* `dirichlet`. Sets the value of the scalar, controlled by the `value` keyword.
  ```json
  {
    "type": "dirichlet",
    "value": 1,
    "zone_indices": [1, 2]
  }
  ```
* `neumann`. Sets the flux of the scalar, controlled by the `flux` keyword.
  ```json
  {
    "type": "neumann",
    "flux": 1,
    "zone_indices": [1, 2]
  }
  ```
* `user_pointwise`. Sets the scalar in the pointwise user interface routine.
  ```json
  {
    "type": "user_poinwise",
    "zone_indices": [1, 2]
  }
  ```
* `user`. User boundary condition, see [further documentation](#user-file_field-dirichlet-update).
  ```json
  {
    "type": "user",
    "zone_indices": [1, 2]
  }
  ```

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
   `case.point_zones` object. See more about [point zones](@ref point-zones).
4. `field`, where the initial condition is retrieved from a field file. Works
   in the same way as for the fluid. See the
   [fluid section](@ref case-file_fluid-ic) for detailed explanations.

### Source terms

The configuration of source terms is the same as for the fluid. A demonstration
of using source terms for the scalar can be found in the `scalar_mms` example.

### Linear solver configuration

Should be provided as an object under the `solver` keyword. For available
configuration options, see the corresponding documentation for the fliud. A
standard choice would be `"type": "cg"` and `"preconditioner": "jacobi"`.

### Full parameter table

| Name                           | Description                                                       | Admissible values                           | Default value |
| ------------------------------ | ----------------------------------------------------------------- | ------------------------------------------- | ------------- |
| `enabled`                      | Whether to enable the scalar computation.                         | `true` or `false`                           | `true`        |
| `name`                         | The name associated to the scalar solver.                         | String                                      | `scalar`      |
| `field_name`                   | The name of the solution in the field registry.                   | A string                                    | `s`           |
| `Pe`                           | The Peclet number.                                                | Positive real                               | -             |
| `cp`                           | Specific heat capacity.                                           | Positive real                               | -             |
| `lambda`                       | Thermal conductivity.                                             | Positive real                               | -             |
| `nut_field`                    | Name of the turbulent kinematic viscosity field.                  | String                                      | Empty string  |
| `Pr_t`                         | Turbulent Prandtl number                                          | Positive real                               | -             |
| `boundary_types`               | Boundary types/conditions labels.                                 | Array of strings                            | -             |
| `initial_condition.type`       | Initial condition type.                                           | `user`, `uniform`, `point_zone`             | -             |
| `initial_condition.value`      | Value of the velocity initial condition.                          | Real                                        | -             |
| `source_terms`                 | Array of JSON objects, defining additional source terms.          | See list of source terms above              | -             |
| `gradient_jump_penalty`        | Array of JSON objects, defining additional gradient jump penalty. | See list of gradient jump penalty above     | -             |
| `advection`                    | Whether to compute the advetion term.                             | `true` or `false`                           | `true`        |
| `solver.type`                  | Linear solver for scalar equation.                                | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres` | -             |
| `solver.preconditioner.type`   | Linear solver preconditioner for the momentum equation.           | `ident`, `hsmg`, `jacobi`                   | -             |
| `solver.absolute_tolerance`    | Linear solver convergence criterion for the momentum equation.    | Positive real                               | -             |
| `solver.maxiter`               | Linear solver max iteration count for the momentum equation.      | Positive real                               | 800           |
| `solver.projection_space_size` | Projection space size for the scalar equation.                    | Positive integer                            | 0             |
| `solver.projection_hold_steps` | Holding steps of the projection for the scalar equation.          | Positive integer                            | 5             |


## Simulation components
Simulation components enable the user to perform various additional operations,
which are not strictly necessary to run the solver. An example could be
computing and output of additional fields, e.g. vorticity.

A more detailed description as well as a  full list of available components and
 their setup is provided in a [separate page of the manual](@ref simcomps).

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
 their setup is provided in a [separate page of the manual](@ref point-zones).

## Runtime statistics

This object adds the collection of runtime statistics (timings) for identified
profiling regions. A region is defined as all functions between a call to
`profiler_start_region(name, id)` and `profiler_end_region(name, id)`. Neko
currently supports 50 regions, with id 1..25 being reserved for internal use.


| Name             | Description                                                 | Admissible values | Default value |
| ---------------- | ----------------------------------------------------------- | ----------------- | ------------- |
| `enabled`        | Whether to enable gathering of runtime statistics           | `true` or `false` | `false`       |
| `output_profile` | Whether to output all gathered profiling data as a CSV file | `true` or `false` | `false`       |
