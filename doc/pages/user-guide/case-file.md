# Case File {#case-file}

\tableofcontents

The case file defines all the parameters of a simulation.
The format of the file is JSON, making it easy to read and write case files
using the majority of the popular programming languages.
JSON is hierarchical and consists of parameter blocks enclosed in curly
braces.
These blocks are referred to as objects.
The case file makes use of objects to separate the configuration of different parts
 of the solver.
We refer the reader to the examples shipped with the code to get a good
idea of how a case file looks.
The table below provides a complete reference for all possible configuration
choices.

## High-level structure
The current high-level structure of the case file is shown below.

~~~~~~~~~~~~~~~{.json}
{
    "version": 1.0,
    "case": {
        "constants": [],
        "time": {},
        "numerics": {},
        "fluid": {},
        "scalar": {},
        "simulation_components" : [],
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
The subsections below list all the configuration options for each of the
high-level objects. Some parameters will have default values, and are therefore
optional.

## Output frequency control
A common scheme for controlling the output frequency is applied for various
outputs.
It is described already now in order to clarify the meaning of several
parameters found in the tables below.

The frequency is controlled by two parameters, ending with `_control` and
`_value`, respectively.
The latter name is perhaps not ideal, but it is somewhat difficult to come up
with a good one, suggestions are welcome.

The `_value` parameter is a *real* number, that defines the output frequency,
but the interpretation of that number depends on the choice of `_control`. The
three following options are possible.
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
| `output_format`       | The file format of field data.                                                                        | `nek5000`, `adios2`, or `vtkhdf`                | `nek5000`     |
| `output_precision`    | Whether to output snapshots in single or double precision                                             | `single` or `double`                            | `single`      |
| `output_layout`       | Data layout for `adios2` files. (Choose `2` or `3` for ADIOS2 supported compressors BigWhoop or ZFP.) | Positive integer `1`, `2`, `3`                  | `1`           |
| `load_balancing`      | Whether to apply load balancing.                                                                      | `true` or `false`                               | `false`       |
| `output_partitions`   | Whether to write a `partitions.vtk` file with domain partitioning.                                    | `true` or `false`                               | `false`       |
| `output_checkpoints`  | Whether to output checkpoints, i.e. restart files.                                                    | `true` or `false`                               | -             |
| `checkpoint_control`  | Defines the interpretation of `checkpoint_value` to define the frequency of writing checkpoint files. | `nsamples`, `simulationtime`, `tsteps`, `never` | -             |
| `checkpoint_value`    | The frequency of sampling in terms of `checkpoint_control`.                                           | Positive real or integer                        | -             |
| `checkpoint_filename` | The filename of written checkpoint.                                                                   | Strings such as `my_name`                       | `fluid`       |
| `checkpoint_format`   | The file format of checkpoints                                                                        | `chkp` or `hdf5`                                | `chkp`        |
| `restart_file`        | checkpoint to use for a restart from previous data                                                    | Strings ending with `.chkp`                     | -             |
| `restart_mesh_file`   | If the restart file is on a different mesh, specify the .nmsh file used to generate it here           | Strings ending with `.nmsh`                     | -             |
| `mesh2mesh_tolerance` | Tolerance for the restart when restarting from another mesh                                           | Positive reals                                  | 1e-6          |
| `job_timelimit`       | The maximum wall clock duration of the simulation.                                                    | String formatted as HH:MM:SS                    | No limit      |
| `output_at_end`       | Whether to always write all enabled output at the end of the run.                                     | `true` or `false`                               | `true`        |

Some additional practical comments are provided regarding the output triggered
by `job_timelimit` and `output_at_end` keywords.

If `output_at_end` is set to `true`, an additional write is performed after the
execution of the simulation time-loop is finished. This triggers most outputs,
like the fluid solvers, the checkpoint, etc. Note that if your case settings are
such that a particular output is written at the last time step regardless of
`output_at_end` (e.g. `end_time: 5`, `checkpoint_value: 5`,
 `checkpoint_control: simulationtime` ) you will get two outputs with the same
values: one from your ordinary write and one triggered by `output_at_end`.

@note This has a rather detrimental effect on outputs from various
statistics-related [simulation components](@ref simcomps). Since the collected
statistics are reset on write, the data written by `output_at_end` will be just
zeroes.

The purpose of `job_timelimit` is to gracefully stop the simulation in a typical
supercomputer environment, where your runtime is limited. When Neko detects that
the time of the run exceeds the `job_timelimit`, it exits the time-loop. At this
point, if one sets `output_at_end` to `true`, this will trigger a write as per
usual. However, if `output_at_end` is `false`, Neko will still write a special
checkpoint file, with the filename called `joblimit#####.chkp`. This is done so
that the user is at least provided a restart file, and none of the computer time
spent on the simulation is wasted. Generally, however, it is recommended to
have `output_at_end` set to `true` in tandem with `job_timelimit`, so that what
exactly gets written is controlled by the case file settings.

### Constants
The `constants` array allows the user to define parameters that are global to
the case file, and can be referred to when setting the values of other
parameters. Two types of parameters can be defined: scalars and arrays. Each is
represented as a subobject inside the `constants` object and should containt two
entries: `name` and `value`. Here is an example:

```json
"constants":
[
  {
    "name": "const1",
    "value": 3.5
  },
  {
    "name": "vector1",
    "value": [1, 0, 1]
  }
]
```

Other parameters in the case file that require a scalar or array entry, can
instead be defined as a string, pointing to the name of the corresponding
parameter in the `constants` object. As an example, recall that output frequency
is controlled by the keyword `output_value`. It is a plausible scenario that the
frequency is the same for multiple solvers, simulation components, etc. Assuming
a simulation with both [fluid](@ref case-file_fluid) and [scalar](@ref
case-file_scalar) solvers active, the following could be used.

```json
"constants":
[
  {
    "name": "common_output_value",
    "value": 10
  }
],
"fluid":
{
  "output_value": "common_output_value"
},
"scalar":
{
  "output_value": "common_output_value"
}
```
The advantage is that this guarantees that the fluid and scalar output will be
in sync, and if one wants to change the frequency only does that in one place in
the case file. Another use case is demonstrated in the `hemi` example, where the
freestream velocity is defined under `constants` and then used to setup both
initial and boundary conditions.

Under the hood, Neko stores the constants in an object called
`neko_const_registry`, which is of the type `registry_t` (same as
`neko_registry`). The object is accessible in the [user file](@ref user-file).

### Time control
The `time` object is used to define the time-stepping of the simulation,
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
| no_slip (stationary wall)                        | 1   |
| velocity_value                  | 2   |
| outflow, normal_outflow (+dong) | 3   |
| symmetry                        | 4   |
| periodic                        | 6   |
| user_velocity                   | 7   |
| user_pressure                   | 8   |
| shear_stress                    | 9   |
| wall_model                      | 10  |
| blasius_profile                 | 11  |
| no_slip (moving wall)           | 12  |

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

## Fluid {#case-file_fluid}

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

### Compressible flows

Neko supports compressible flow simulations via the compressible solver.
To enable compressible flow, set `"scheme": "compressible"` in the fluid
configuration. This solver integrates the compressible Euler equations (full
Navier-Stokes will be enabled in upcoming updates) using a Runge-Kutta time
integration scheme with artificial viscosity for stability.

The compressible solver requires the following parameters:

| Name    | Description                           | Admissible values | Default value |
| ------- | ------------------------------------- | ----------------- | ------------- |
| `gamma` | Ratio of specific heats for ideal gas | Positive reals    | `1.4`         |

Additional numerics parameters specific to compressible flows:

| Name              | Description                                        | Admissible values | Default value |
| ----------------- | -------------------------------------------------- | ----------------- | ------------- |
| `c_avisc_low`     | Coefficient for low-order artificial viscosity     | Positive reals    | `0.5`         |
| `c_avisc_entropy` | Coefficient for entropy-based artificial viscosity | Positive reals    | `1.0`         |

The compressible solver uses variable time-stepping controlled by the CFL
number. Set `variable_timestep` to `true` and specify `target_cfl` in the time
control object.

Example configuration:
~~~~~~~~~~~~~~~{.json}
{
  "fluid": {
    "scheme": "compressible",
    "gamma": 1.4,
    "initial_condition": {
      "type": "user"
    },
    "boundary_conditions": [
      {
        "type": "velocity_value",
        "zone_indices": [1],
        "value": [3, 0, 0]
      },
      {
        "type": "density_value",
        "zone_indices": [1],
        "value": 1.4
      },
      {
        "type": "pressure_value",
        "zone_indices": [1],
        "value": 1
      }
    ],
    "output_control": "nsamples",
    "output_value": 20
  },
  "numerics": {
    "time_order": 3,
    "polynomial_order": 5,
    "c_avisc_low": 0.5,
    "c_avisc_entropy": 0.5
  }
}
~~~~~~~~~~~~~~~

#### Compressible boundary conditions

The compressible solver supports the following boundary conditions:

| Boundary Condition | Description                               |
| ------------------ | ----------------------------------------- |
| velocity_value     | Dirichlet condition for velocity (inflow) |
| density_value      | Dirichlet condition for density           |
| pressure_value     | Dirichlet condition for pressure          |
| no_slip            | Zero velocity wall                        |
| symmetry           | Symmetry plane                            |
| outflow            | Pressure outlet (zero gradient)           |
| normal_outflow     | Normal outflow condition                  |

For examples of compressible flow setups, see the `euler_1d_sod`,
`euler_2d_forward_facing_step`, and `euler_tgv` examples.

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
velocity components to be solved in a coupled manner. Therefore, the `coupled_cg`
(or `fused_coupled_cg`) solver should be used for velocity.

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

It is possible to assign specific names to the boundary conditions through the
`name` keyword. Boundary conditions can then be retireved in the code by using
the name or the `zone_index` where it is applied.

The default name of the boundary conditions is given by the `<variable>_bc_<zone_index>`
pattern. i.e., the pressure boundary condition that applies in zone index 5 can be
retrieved by the `pressure_bc_5` name.

#### Available conditions
The conditions to apply is specified by `type` keyword inside each of the JSON
objects. The full list of possible conditions for the fluid is specified in the
table below.

| Boundary Condition  | Description                                                                                                                                            |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| symmetry            | A symmetry plane. Must be axis-aligned.                                                                                                                |
| velocity_value      | A Dirichlet condition for velocity.                                                                                                                    |
| no_slip             | A no-slip wall. It can be stationary or moving.                                                                                                                                         |
| outflow             | A pressure outlet.                                                                                                                                     |
| normal_outflow      | An Neumann condition for the surface-normal component of velocity combined with a Dirichlet for the surface-parallel components. Must be axis-aligned. |
| outflow+user        | Same as `outflow` but with user-specified pressure.                                                                                                    |
| normal_outflow+user | Same as `normal_outflow` but with user-specified pressure.                                                                                             |
| outflow+dong        | A pressure outlet with the Dong condition applied.                                                                                                     |
| normal_outflow+dong | The `normal_outflow` with the Dong condition applied. Must be axis-aligned.                                                                            |
| shear_stress        | Prescribed wall shear stress. Must be axis-aligned.                                                                                                    |
| wall_model          | Shear stress condition based on a wall model for large-eddy simulation.                                                                                |
| blasius_profile     | A Blasius velocity profile.                                                                                                                            |
| user_velocity       | The `field_dirichlet_vector_t` user-defined Dirichlet condition for velocity.                                                                          |
| user_pressure       | The `field_dirichlet_t` user-defined Dirichlet condition for pressure.                                                                                 |
| overset_interface   | A Dirichlet condition that prescribes values from another neko simulation running concurrently.                                                        |

A more detailed description of each boundary condition is provided below.

* `symmetry`. A symmetry plane that must be axis-aligned. Sets the
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
* `no_slip`. A standard no-slip wall, which sets velocity to zero relative to the wall. For moving walls, setting the optional argument `"moving": true` is required. This also requires setting up the ALE module separately. further details can be found in [ALE user guide](#case-file_fluid-ale). For stationary walls, no additional keyword is needed.
  ```json
  {
    "type": "no_slip",
    "zone_indices": [1, 2],
    "moving": false
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

   * The `most` model is a version of the `rough_log_law` adapted for flows with temperature stratification, such as atmospheric boundary layer (ABL) flows. The model uses Monin-Obukhov stability theory (MOST) to account for the local temperature gradient. More details and required keywords are given [below](#most-wall-model).

   * The `richardson` model is similar to the `most` model, but it assesses the stability dependence based on the Richardson number instead of the Obukhov length. More details and required keywords are given [below](#richardson-wall-model).

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
* `overset_interface`, a Dirichlet boundary condition that retrieves values
  from another simulation with an overlapping domain. For this case, it is
  recommended that all zone indices that need to be considered as an overset
  interface are included in one boundary. This avoids repeated calls to
  interpolation routines.

  As a note, both meshes must be overlapping with at least one element.

  Since this requires another concurrent simulation, you must execute neko in
  [multiple-program-multiple-data (MPMD)](#user-file_tips_mpmd) mode. Note that both simulations are otherwise
  independent, therefore the rest of the user case can be modified as seen fit.

  The *time-step* must be the same between the simulations. A variable time-step is not
  supported at the moment.

  ```json
  {
    "type": "overset_interface",
    "zone_indices": [1, 2]
  }
  ```

#### MOST wall model
The `most` model is based on Monin-Obukhov similarity theory (Monin and Obukhov, 1954) and adds a correction to the rough log law according to

\f{eqnarray*}{
   \frac{\partial{V}}{\partial z} &=& \frac{u_*}{\kappa z}\phi_m\left(\frac{z}{L}\right), \\
   \frac{\partial{\theta}}{\partial z} &=& \frac{\overline{(w'\theta')}}{u_* \kappa z}\phi_h\left(\frac{z}{L}\right),
 \f}

 where \f$V\f$ is the horizontal wind speed (given that \f$z\f$ is the wall-normal direction) and \f$\theta\f$ is the potential temperature.

 The formulations of the correction functions \f$\phi_m\f$ and \f$\phi_h\f$ are taken from Dyer 1974 for the convective regime, and from Holstlag and De Bruin 1988 for the stable regime.

 The keywords for this model are:
 - `kappa`: The von Kàrmàn constant, defaults to 0.4 (as is the standard in the ABL literature).

 - `Pr`: The turbulent Prandtl number, defaults to 1.0.
 - `z0`: The characteristic roughness length for momentum.
 - `z0h`: The characteristic roughness length for heat. If a negative value is given, the roughness length for heat is computed using the formula of Zilitinkevich 1995, with the provided value acting as the constant \f$-A_0\f$ in the Zilitinkevich formula. Defaults to be the same as `z0`.
 - `type_of_temp_bc`: Accepted values are the lowercase strings `neumann` or `dirichlet`. If `neumann`, the provided value of `bottom_bc_flux_or_temp` is used directly as the surface heat flux in the computation of the wall stress. If `dirichlet`, the value of `bottom_bc_flux_or_temp` is interpreted as a surface temperature, which is then used to compute a heat flux using the MOST relationship.
 - `bottom_bc_flux_or_temp`: Value of the surface heat flux if `type_of_temp_bc` is `neumann`, or value of the surface temperature if `type_of_temp_bc` is `dirichlet`.
 - `scalar_field`: The name of the scalar field to be used as the potential temperature in the equations.
 - `time_dependent_temp_bc`: Boolean. If `false` the value of `bottom_bc_flux_or_temp` will be kept constant throughout the simulation. If `true`, the wall model will look for `bc_value` in `neko_const_registry` and assign that value at each time step. The value of `bc_value` can then be updated in the user file, for example in `user_check`.
 <details>
  <summary><b><u>Example of user file implementation</u></b></summary>

```fortran
   subroutine user_check(time)
      type(time_state_t), intent(in) :: time
      real(kind=rp), pointer :: bc_value

      bc_value => neko_const_registry%get_real_scalar("bc_value")

      bc_value = scalar_bc

   end subroutine user_check

  subroutine dirichlet_update(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    integer i

      if (fields%items(1)%ptr%name .eq. "temperature") then

       associate(s => fields%items(1)%ptr)
            do i = 1, bc%msk(0)
               s%x(bc%msk(i), 1, 1, 1) = scalar_bc(time)
            end do
            if (neko_bcknd_device .eq. 1) then
               call device_memcpy(s%x, s%x_d, s%size(), &
                     host_to_device, sync=.false.)
            end if
         end associate
      end if
   end subroutine dirichlet_update

   function scalar_bc(time) result(bc)
      type(time_state_t), intent(in) :: time
      real(kind=rp) :: bc

      bc = 265.0_rp - 0.25_rp/3600.0_rp*time%t

   end function scalar_bc
```

</details>

 @attention This wall model uses a `neumann` or `dirichlet` value for the scalar field to compute the surface shear stress, but it does not set the boundary condition for the scalar. The same boundary condition should be set separately for the scalar (see [Boundary conditions](#boundary-conditions)).

  <details>
  <summary><b><u>Example code snippet</u></b></summary>

  ```json
  {
    "type": "wall_model",
    "model": "most",
    "kappa": 0.4,
    "Pr": 1.0,
    "z0": 0.1,
    "z0h": 0.1,
    "type_of_temp_bc": "neumann",
    "bottom_bc_flux_or_temp": 0.05,
    "scalar_field": "temperature",
    "time_dependent_temp_bc": "false",
    "zone_indices": [5],
    "h_index": 1
  }
  ```

  </details>

   <details>
   <summary><b><u>References</u></b></summary>

 Dyer, A. J. (1974). A review of flux-profile relationships. Boundary-Layer Meteorology, 7(3), 363–372. https://doi.org/10.1007/BF00240838

  Holtslag, A. A. M., & De Bruin, H. A. R. (1988). Applied Modeling of the Nighttime Surface Energy Balance over Land. Journal of Applied Meteorology, 27(6), 689–704. https://doi.org/10.1175/1520-0450(1988)027%253C0689:AMOTNS%253E2.0.CO;2

  Monin, A. S., & Obukhov, A. M. (1954). Basic laws of turbulent mixing in the surface layer of the atmosphere. Tr Akad Nauk SSSR Geofiz Inst, 24(151), 163–187.

  Zilitinkevich, S. S., 1995: Non-local turbulent transport: Pollution dispersion aspects of coherent structure of convective flows. Air Pollution III, H. Power, N. Moussiopoulos, and C. A. Brebbia, Eds., Vol. 1, Air Pollution Theory and Simulation, Computational Mechanics Publications, 53–60.
</details>


#### MOST wall model
The `most` model is based on Monin-Obukhov similarity theory (Monin and Obukhov, 1954) and adds a correction to the rough log law based on the Obukhov length \f$L\f$, according to

\f{eqnarray*}{
   \frac{\partial{V}}{\partial z} &=& \frac{u_*}{\kappa z}\phi_m\left(\frac{z}{L}\right), \\
   \frac{\partial{\theta}}{\partial z} &=& \frac{\overline{(w'\theta')}}{u_* \kappa z}\phi_h\left(\frac{z}{L}\right),
 \f}

 where \f$V\f$ is the horizontal wind speed (given that \f$z\f$ is the wall-normal direction) and \f$\theta\f$ is the potential temperature.

 The formulations of the correction functions \f$\phi_m\f$ and \f$\phi_h\f$ are taken from Dyer 1974 for the convective regime, and from Holstlag and De Bruin 1988 for the stable regime.

 The keywords for this model are:
 - `kappa`: The von Kàrmàn constant, defaults to 0.4 (as is the standard in the ABL literature).

 - `Pr`: The turbulent Prandtl number, defaults to 1.0.
 - `z0`: The characteristic roughness length for momentum.
 - `z0h`: The characteristic roughness length for heat. If a negative value is given, the roughness length for heat is computed using the formula of Zilitinkevich 1995, with the provided value acting as the constant \f$-A_0\f$ in the Zilitinkevich formula. Defaults to be the same as `z0`.
 - `type_of_temp_bc`: Accepted values are the lowercase strings `neumann` or `dirichlet`. If `neumann`, the provided value of `bottom_bc_flux_or_temp` is used directly as the surface heat flux in the computation of the wall stress. If `dirichlet`, the value of `bottom_bc_flux_or_temp` is interpreted as a surface temperature, which is then used to compute a heat flux using the MOST relationship.
 - `bottom_bc_flux_or_temp`: Value of the surface heat flux if `type_of_temp_bc` is `neumann`, or value of the surface temperature if `type_of_temp_bc` is `dirichlet`.
 - `scalar_field`: The name of the scalar field to be used as the potential temperature in the equations.
 - `time_dependent_temp_bc`: Boolean. If `false` the value of `bottom_bc_flux_or_temp` will be kept constant throughout the simulation. If `true`, the wall model will look for `bc_value` in `neko_const_registry` and assign that value at each time step. The value of `bc_value` can then be updated in the user file, for example in `user_check`.
 <details>
  <summary><b><u>Example of user file implementation</u></b></summary>

```fortran
   subroutine user_check(time)
      type(time_state_t), intent(in) :: time
      real(kind=rp), pointer :: bc_value

      bc_value => neko_const_registry%get_real_scalar("bc_value")

      bc_value = scalar_bc

   end subroutine user_check

  subroutine dirichlet_update(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    integer i

      if (fields%items(1)%ptr%name .eq. "temperature") then

       associate(s => fields%items(1)%ptr)
            do i = 1, bc%msk(0)
               s%x(bc%msk(i), 1, 1, 1) = scalar_bc(time)
            end do
            if (neko_bcknd_device .eq. 1) then
               call device_memcpy(s%x, s%x_d, s%size(), &
                     host_to_device, sync=.false.)
            end if
         end associate
      end if
   end subroutine dirichlet_update

   function scalar_bc(time) result(bc)
      type(time_state_t), intent(in) :: time
      real(kind=rp) :: bc

      bc = 265.0_rp - 0.25_rp/3600.0_rp*time%t

   end function scalar_bc
```

</details>

 @attention This wall model uses a `neumann` or `dirichlet` value for the scalar field to compute the surface shear stress, but it does not set the boundary condition for the scalar. The same boundary condition should be set separately for the scalar (see [Boundary conditions](#boundary-conditions)).

  <details>
  <summary><b><u>Example code snippet</u></b></summary>

  ```json
  {
    "type": "wall_model",
    "model": "most",
    "kappa": 0.4,
    "Pr": 1.0,
    "z0": 0.1,
    "z0h": 0.1,
    "type_of_temp_bc": "neumann",
    "bottom_bc_flux_or_temp": 0.05,
    "scalar_field": "temperature",
    "time_dependent_temp_bc": false,
    "zone_indices": [5],
    "h_index": 1
  }
  ```

  </details>

   <details>
   <summary><b><u>References</u></b></summary>

 Dyer, A. J. (1974). A review of flux-profile relationships. Boundary-Layer Meteorology, 7(3), 363–372. https://doi.org/10.1007/BF00240838

  Holtslag, A. A. M., & De Bruin, H. A. R. (1988). Applied Modeling of the Nighttime Surface Energy Balance over Land. Journal of Applied Meteorology, 27(6), 689–704. https://doi.org/10.1175/1520-0450(1988)027%253C0689:AMOTNS%253E2.0.CO;2

  Monin, A. S., & Obukhov, A. M. (1954). Basic laws of turbulent mixing in the surface layer of the atmosphere. Tr Akad Nauk SSSR Geofiz Inst, 24(151), 163–187.

  Zilitinkevich, S. S., 1995: Non-local turbulent transport: Pollution dispersion aspects of coherent structure of convective flows. Air Pollution III, H. Power, N. Moussiopoulos, and C. A. Brebbia, Eds., Vol. 1, Air Pollution Theory and Simulation, Computational Mechanics Publications, 53–60.
</details>

### Richardson wall model {#richardson-wall-model}
This Richardson-number based wall model is conceptually similar to the more well-known MOST-based wall model, but it computes the effect of the temperature stratification based on the bulk Richardson number instead of the Obukhov length.

In the convective regime, the surface shear stress, \f$\tau\f$, and surface heat flux, \f$\overline{u'\theta'}\f$ are computed using the formulations of Louis 1979:
\f{eqnarray*}{
\tau &=& a^2 V^2 F_m\left(\frac{z}{z_0}, \mathrm{Ri}_b\right), \\
\overline{u'\theta'} &=& \frac{a^2}{R}\, V\, \Delta\theta \, F_h\left(\frac{z}{z_{0h}}, \mathrm{Ri}_b\right).
\f}

Here, \f$V\f$ is the horizontal wind speed (given that \f$z\f$ is the wall-normal direction); \f$\theta\f$ is the potential temperature; \f$\mathrm{Ri}_b\f$ is the bulk Richardson number; \f$z_0\f$ and \f$z_{0h}\f$ are the roughness lengths for momentum and heat, respectively; \f$F_m\f$ and \f$F_h\f$ are stability functions as defined in Louis 1979; and \f$a\f$ and \f$R\f$ are constants, also as defined in Louis 1979.

In the stable regime, the surface shear stress and surface heat flux are computed based on Mauritsen et al. 2007:
\f{eqnarray*}{
\tau &=& \frac{V^2}{\left[\ln\left(\dfrac{z}{z_0}\right)\right]^2} \,\frac{f_{\tau}(\mathrm{Ri}_b)}{f_{\tau}(0)} \left(\frac{\ell}{z}\right)^2, \\
\overline{u'\theta'} &=& \frac{\Delta\theta}{\ln\left(\dfrac{z}{z_{0h}}\right)} \,\frac{f_{\theta}(\mathrm{Ri}_b)}{\left|f_{\theta}(0)\right|} \left(\frac{\ell}{z}\right) \frac{u_*}{\mathrm{Pr}}.
\f}

Here, \f$V, \theta, \mathrm{Ri}_b, z_0\f$, and \f$z_{0h}\f$ are the same as above; \f$f_{\tau}\f$ and \f$f_{\theta}\f$ are defined in Mauritsen et al. 2007; \f$l\f$ is a lengthscale (we use \f$l = \kappa z\f$, where \f$\kappa=0.4\f$ is the von Kàrmàn constant); \f$u_*\f$ is the friction velocity, and \f$\mathrm{Pr}\f$ is the turbulent Prandtl number.

The keywords for this wall model are the same as for the [MOST model](#most-wall-model), and a time-varying temperature boundary condition can be applied in the same way as described for the MOST model.


   <details>
   <summary><b><u>References</u></b></summary>
  Louis, J.-F. (1979). A parametric model of vertical eddy fluxes in the atmosphere. Boundary-Layer Meteorology, 17(2), 187–202. https://doi.org/10.1007/BF00117978.

  Mauritsen, T., Svensson, G., Zilitinkevich, S. S., Esau, I., Enger, L., & Grisogono, B. (2007). A Total Turbulent Energy Closure Model for Neutrally and Stably Stratified Atmospheric Boundary Layers. Journal of the Atmospheric Sciences, 64(11), 4113–4126. https://doi.org/10.1175/2007JAS2294.1.
  </details>

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
   | Name             | Description                                                                                        | Admissible values            | Default value  |
   | ---------------- | -------------------------------------------------------------------------------------------------- | ---------------------------- | -------------- |
   | `file_name`      | Name of the field file to use (e.g. `myfield0.f00034`).                                            | Strings ending with `f*****` | -              |
   | `interpolate`    | Whether to interpolate the velocity and pressure fields from the field file onto the current mesh. | `true` or `false`            | `false`        |
   | `mesh_file_name` | If interpolation is enabled, the name of the field file that contains the mesh coordinates.        | Strings ending with `f*****` | `file_name`    |
   | `interpolation.tolerance`| Tolerance for the point search.                                                            | Positive real.               | `NEKO_EPS*1e3` |
   | `interpolation.padding`  | Padding for the point search.                                                              | Positive real.               | `1e-2`         |

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
   Neko writes single precision `fld` files by default. To write your
   files in double precision, set `case.output_precision` to
   `"double"`.

   @attention Neko does not automatically detect if interpolation is needed.
   Interpolation will always be performed if `"interpolate"` is set
   to `true`, even if the field file matches with the current simulation.


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

@note Notice that to perform simulation in a rotating reference frame one has to
define both `coriolis` and `centrifugal` source terms in a consistent way.

5. `user`, the values are set inside the compiled user file, using the
   corresponding user file subroutine.
6. `brinkman`, Brinkman permeability forcing inside a pre-defined region.
7. `gradient_jump_penalty`, perform gradient_jump_penalisation.
8. `sponge`, adds a sponge term based on a reference velocity field, which is
   applied in a user-specified region of the domain.
9. `field`, uses fields in the `neko_registry` as values of the source term. The
   fields are selected with the `field_names` keyword.

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

#### Sponge

The sponge source term adds a term to each of the momentum equations of the form

\f$ \mathbf{\lambda} f(\mathbf{x}) ( \mathbf{u}^{bf} - \mathbf{u}) \f$

where:

- \f$ \mathbf{\lambda} \f$ is a 3-element vector of amplitudes of the sponge forcing in each Cartesian direction,
- \f$ \mathbf{u}^{\text{bf}} \f$ is a reference (baseflow) velocity field,
- \f$ f(\mathbf{x}) \f$ is a user-defined sponge mask field, defining where the sponge is active.

Amplitudes are specified using the `amplitudes` keyword with an array of
3 reals. Any of those values can be set to 0 to suppress the forcing in that
particular direction. For example `[1.0, 1.0, 0.0]` will multiply the fringe
field by 1, 1, and 0 in the `x`, `y` and `z` directions respectively,
effectively removing the forcing in the `z` direction.

The reference velocity field, or `baseflow` can be set from three methods:
1. `constant`, applies constant values according to the `values` keyword:

   <details>
   <summary><b><u>Example code snippet</u></b></summary>
   ```json
   "source_terms": [
      {
         "type": "sponge",
         "amplitudes": [1.0, 1.0, 1.0],
         "baseflow": {
             "method": "constant",
             "value": [2.0, 0.0, 0.0]
         }
      }
   ]
   ```
   </details>

2. `field`, where the velocity fields are retrieved from an `fld` file.
   Uses the same parameters as the field initial condition.
   @note The same parameters as the `field` initial condition apply here.

   <details>
   <summary><b><u>Example code snippet</u></b></summary>
   ```json
   "source_terms": [
      {
         "type": "sponge",
         "amplitudes": [1.0, 1.0, 1.0],
         "baseflow": {
             "method": "field",
             "file_name": "my_field0.f00016",
             "mesh_file_name": "my_field0.f00000",
             "interpolate": true
         }
      }
   ]
   ```
   </details>

3. `user`, where the velocity field is set according to what
   is defined in the user file. Useful for setting
   velocity fields manually. In this case, the base flow fields must be
   created and added to the `neko_registry` (see fortran code snippet
   below).
   <details>
   <summary><b><u>Example code snippet</u></b></summary>
   ```json
   "source_terms": [
      {
         "type": "sponge",
         "amplitudes": [1.0, 1.0, 1.0],
         "baseflow": {
             "method": "user"
         }
      }
   ]
   ```
   </details>

Finally, the fringe function field must be filled by the user. This must be
done through the user file by adding the fringe field to the
`neko_registry` in either `initialize` or `initial_conditions` (more
specifically, before the first call to compute the sponge source term). Note that `initial_conditions` is not called when doing a restart from a checkpoint file, so if restarts will be done the sponge should be implemented in `user_init_modules`.

The fringe field must be set by adding a field to the `neko_registry`
under a specific name that can be retrieved internally. By default, Neko will
search for the field `"sponge_fringe"` in the registry, but this can be changed
by setting the parameter `fringe_registry_name`, which is important when using
more than one sponge source term.

The same principle applies for the base flow fields (if `"method": "user"`).
By default, neko will search for the base flow fields in the registry using
the prefix `"sponge_bf_"`, meaning that `u` will be in `sponge_bf_u`, etc.
This prefix can be changed by setting the parameter `bf_registry_prefix`.

<details>
<summary><b><u>Example using `initialize`</u></b></summary>

```fortran
module user
  use neko
  implicit none

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initialize => user_initialize
  end subroutine user_setup

  ! User-defined initialization called just before time loop starts
  subroutine user_initialize(t)
    type(time_state_t), intent(in) :: t

    type(field_t), pointer :: u, fringe, ubf, vbf, wbf
    real(kind=rp) :: x, y, xmin1, delta_rise1, xmin2, delta_rise2
    integer :: i, imask

    !
    ! 1. Add the "sponge_field" to the field registry.
    !    NOTE: The name of the fringe field in the registry
    !    can be changed with the parameter `fringe_registry_name`.
    !
    !
    u => neko_registry%get_field("u")
    call neko_registry%add_field(u%dof,"sponge_fringe")
    fringe => neko_registry%get_field("sponge_fringe")

    ! Initialize the base flows
    call neko_registry%add_field(u%dof,"sponge_bf_u")
    ubf => neko_registry%get_field("sponge_bf_u")
    call neko_registry%add_field(u%dof,"sponge_bf_v")
    vbf => neko_registry%get_field("sponge_bf_v")
    call neko_registry%add_field(u%dof,"sponge_bf_w")
    wbf => neko_registry%get_field("sponge_bf_w")

    !
    ! 2. Set the function f(x,y,z) from 0 to 1. in two zones of the mesh,
    !    a top region in x \in [xmin1, +\infty[, y \in [0, +\infty[
    !  and a bottom region x \in [xmin2, +\infty[, y \in ]-\infty, 0[
    !
    !    A smoothing function S(x) is applied at the beginning of each zone,
    !    with a rising distance of delta_rise1 and delta_rise2
    !

    ! Bottom boundary
    xmin1 = 3.0_rp
    delta_rise1 = 3.0_rp

    ! Top boundary
    xmin2 = 20.0_rp
    delta_rise2 = 7.0_rp

    fringe = 0.0_rp
    do i = 1, fringe%size()
        x = fringe%dof%x(i,1,1,1)
        y = fringe%dof%y(i,1,1,1)

        ! Bottom boundary
        if ( (y .lt. 0.0_rp) .and. (x .gt. xmin1)) then
           fringe%x(i,1,1,1) = S( (x - xmin1)/delta_rise1 )

        ! Top boundary
        else if ( (y .gt. 0.0_rp) .and. (x .gt. xmin2)) then
           fringe%x(i,1,1,1) = S( (x - xmin2)/delta_rise2 )
        end if

       ! Set ubf,vbf to something random
       ubf%x(i,1,1,1) = sin(3.1415926_rp*2.0_rp/10.0_rp * x)
       vbf%x(i,1,1,1) = cos(3.1415926_rp*2.0_rp/10.0_rp * y)

    end do

    wbf = 0.0_rp
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(ubf%x, ubf%x_d, ubf%size(), &
            HOST_TO_DEVICE, .false.)
       call device_memcpy(vbf%x, vbf%x_d, vbf%size(), &
            HOST_TO_DEVICE, .false.)
       call device_memcpy(fringe%x, fringe%x_d, fringe%size(), &
            HOST_TO_DEVICE, .false.)
    end if

    ! NOTE: You can dump the fringe field to file using the `dump_fields`
    ! parameter. The fringe field will be stored under `pressure`.

    nullify(fringe)
    nullify(u)
    nullify(ubf)
    nullify(vbf)
    nullify(wbf)

  end subroutine user_initialize

  ! Smooth step function, 0 if x <= 0, 1 if x >= 1, 1/erp(1/(x-1) + 1/x) between 0 and 1
  function S(x) result(y)
    real(kind=rp), intent(in) :: x
    real(kind=rp)             :: y

    if ( x.le.0._rp ) then
       y = 0._rp
    else if ( x.ge.1._rp ) then
       y = 1._rp
    else
       y = 1._rp / (1._rp + exp( 1._rp/(x-1._rp) + 1._rp/x))
    end if

  end function S

end module user
```

</details>

In order to visualize your baseflow and fringe field, you may set
`dump_fields` to `true`. An `fld` file will be written to disk as
`spng_fields.fld`(note, not in `output_directory`) with the fringe field
stored as `pressure`. You may change the name of the field file by setting
`dump_file_name` (must have the extension `fld`).

The parameters for the sponge source term are summarized in the table below:

| Name                       | Description                                                             | Admissible values                 | Default value     |
| -------------------------- | ----------------------------------------------------------------------- | --------------------------------- | ----------------- |
| `amplitudes`               | Sponge forcing strength in each Cartesian direction                     | Array of 3 reals                  | -                 |
| `baseflow.method`          | Method to define the reference (baseflow) velocity                      | `"constant"`, `"field"`, `"user"` | -                 |
| `baseflow.value`           | Velocity vector for constant baseflow                                   | Array of 3 reals                  | -                 |
| `baseflow.file_name`       | File containing baseflow velocity field                                 | String                            | -                 |
| `baseflow.mesh_file_name`  | Mesh file corresponding to the baseflow field                           | String                            | -                 |
| `baseflow.interpolate`     | Whether to interpolate field values to current mesh                     | Boolean                           | `false`           |
| `baseflow.tolerance`       | Tolerance for interpolation convergence                                 | Real                              | -                 |
| `fringe_registry_name`     | Name of the fringe mask field in `neko_registry`                        | String                            | `"sponge_fringe"` |
| `baseflow_registry_prefix` | Prefix of the base flow fields in `neko_registry`                       | String                            | `"sponge_bf"`     |
| `dump_fields`              | If `true`, dumps the fringe and baseflow fields for visualization       | Boolean                           | `false`           |
| `dump_file_name`           | Name of the `fld` file in which to dump the base flow and fringe fields | String ending with `fld`          | `spng_fields.fld` |


### Arbitrary Lagrangian-Eulerian Framework {#case-file_fluid-ale}
Neko supports the simulation of moving walls through the Arbitrary Lagrangian-Eulerian (ALE) framework. The current implementation allows for an arbitrary number of individually moving or deformable walls, collectively referred to as bodies.

@note Currently, only the CPU backend of the ALE framework is supported. GPU acceleration for ALE computations will be available in future updates.

The `"ale"` block in case file is part of the `"fluid"` object, and has the following high-level structure:

~~~~~~~~~~{.json}
{
  "fluid": {

    "ale": {
      "enabled": true,

      "solver": {
        // Linear solver configuration to obtain the base shape for smooth blending function used for mesh deformation
      },

      "mesh_preview": {
        // Settings to preview mesh motion prior to full simulation
      },

      "bodies": [
        // Array of registered moving bodies
      ]
    },
  },

}
~~~~~~~~~~
To run an ALE simulation, the framework must be set up as follows:

* **Boundary Conditions:** Under the `"boundary_condition"` block, any moving wall must be set to `"no_slip"` with `"moving": true`.
* **Enable ALE:** Under the `"ale"` block, the keyword `"enabled"` must be set to `true`.

@attention If the `"ale"` block is present, the `"enabled"` keyword is mandatory.

* **Body Registration:** There must be a strict mapping between moving boundaries and ALE bodies:
  * Any `"no_slip"` boundary marked as `"moving": true` **must** be registered as an ALE body under the `"bodies"` object.
  * Conversely, any boundary registered as an ALE body **must** be defined as a `"no_slip"` boundary with `"moving": true`.
  * A single boundary zone can only be assigned to a maximum of one ALE body.
  * If any of these rules are violated, Neko will not start the simulation loop and will print an informative error message.

In the following, the main blocks of `"ale"` object are explained.

#### Solver {#case-file_fluid-ale-solver}

To smoothly move the mesh during an ALE simulation, we use a global smooth blending function, \f$ \phi_{total} \f$. This function is found by solving the following variable stiffness Laplace equation:

\f{eqnarray*}{
 \nabla \cdot (h(\mathbf{x}) \nabla \phi_{total}) = 0,
\f}

where \f$ h(\mathbf{x}) \f$ is the spatial mesh stiffness. Due to the linearity of the Laplace operator, the total blending function can be decomposed into the sum of individual base shapes \f$ \phi_i \f$ associated with each registered body \f$ i \f$:

\f{eqnarray*}{
 \phi_{total}(\mathbf{x}) = \sum_{i} \phi_i(\mathbf{x}).
\f}

Consequently, the problem reduces to solving a separate Laplace equation for each individual body:

\f{eqnarray*}{
 \nabla \cdot (h(\mathbf{x}) \nabla \phi_i) = 0.
\f}

Currently, the base shapes \f$ \phi_i \f$ are computed only once during initialization, and these same blending functions are used throughout the entire simulation.

During the simulation, the mesh velocity \f$ \mathbf{w}_{mesh}(\mathbf{x}, t) \f$ at any given grid point \f$ \mathbf{x} \f$ is calculated by multiplying the prescribed velocity \f$ \mathbf{v}_i(\mathbf{x}, t) \f$ of each body by its local base shape value \f$ \phi_i(\mathbf{x}) \f$, and summing the contributions across all registered bodies:

\f{eqnarray*}{
 \mathbf{w}_{mesh}(\mathbf{x}, t) = \sum_{i} \phi_i(\mathbf{x}) \mathbf{v}_i(\mathbf{x}, t)
\f}

@note A separate Laplace equation is solved for each registered body. During the solve for body \f$ i \f$, the Dirichlet boundary conditions are set such that \f$ \phi_i = 1 \f$ on its own moving boundary, and \f$ \phi_i = 0 \f$ on all other boundaries (including fixed walls and other moving bodies), with the exception of periodic boundaries. This guarantees that body \f$ i \f$ deforms or moves as intended, while other ALE bodies remain unaffected by this motion, and all other boundaries which define the simulation domain remain fixed. The value of \f$ \phi_i \f$ on periodic boundaries is determined naturally by solving the Laplace equation.

Details regarding the configuration of the mesh stiffness \f$ h(\mathbf{x}) \f$ can be found in the [Mesh Stiffness](#case-file_fluid_ale_stiff_geom) section.

Within the `"solver"` block, the parameters of the linear solver used to solve the Laplace equation are set. This block accepts the following keywords:

| Name | Description | Admissible values | Default value |
| :--- | :--- | :--- | :--- |
| `type` | Type of linear solver for the Laplace equation |  `"cg"`, `"gmres"` | `"cg"` |
| `preconditioner.type` | Type of preconditioner to use |  `"jacobi"`, `"hsmg"`, `"phmg"` | `"jacobi"` |
| `absolute_tolerance` | Absolute tolerance for solver convergence | Positive real | `1.0e-10` |
| `max_iterations` | Maximum number of linear solver iterations | Positive integer | `10000` |
| `monitor` | Monitor residuals in the linear solver | `true` or `false` | `false` |
| `output_base_shape` | Enables output of the base shape field \f$ \phi \f$ | `true` or `false` | `true` |
| `output_stiffness` | Enables output of the computed mesh stiffness field \f$ h(\mathbf{x}) \f$ | `true` or `false` | `false` |

##### Output Files and Diagnostics
If the output flags are enabled, Neko will generate `.fld` files during the initialization phase. These files are highly useful for verifying that the mesh deformation fields and stiffness regions are configured correctly before running the simulation:

* `phi_<body_name>0.f00000`: Generated if `"output_base_shape": true`. Contains the computed base shape \f$ \phi_i \f$ for a specific body.
* `phi_total0.f00000`: Generated if `"output_base_shape": true` **and** there is more than one body registered. Contains the sum of all base shapes (\f$ \phi_{total} = \sum \phi_i \f$).
* `stiffness0.f00000`: Generated if `"output_stiffness": true`. Contains the global spatial mesh stiffness field \f$ h(\mathbf{x}) \f$.

@attention Due to the linearity and the maximum principle of the Laplace equation, the combined base shape field \f$ \phi_{total} \f$ is guaranteed to be strictly bounded between 0 and 1 everywhere in the domain, provided that the solver's `absolute_tolerance` is set appropriately.


@note It is also possible to provide a custom base shape \f$ \phi \f$ using a `user_ale_base_shapes` user subroutine. In this case, the internal Laplace solver is bypassed entirely, even if the custom subroutine is only used for one of the ALE bodies. It is thus up to the user to ensure the validity of the base shape. Setting `"output_base_shape": true` will still write your custom user shapes to `.fld` files, allowing you to easily visualize and debug your custom implementations. More details about implementing this user subroutine can be found [here](#user-file_ale-base-shapes).


#### Mesh preview

One of the available features in the Neko ALE module is the Mesh Preview. This allows users to investigate mesh quality over time **without** running an expensive full fluid simulation. The preview applies the exact same prescribed motion configured by the user, ensuring the visualized mesh motion perfectly mirrors what will happen during the actual simulation. Enabling this feature generates `.fld` files named `mesh_preview0.f*****`, which contain the deformed mesh and the mass matrix.

The `"mesh_preview"` block accepts the following keywords:

| Name | Description | Admissible values | Default value |
| :--- | :--- | :--- | :--- |
| `enabled` | Toggles the mesh preview feature on or off | `true` or `false` | `false` |
| `start_time` | Start time for the mesh preview simulation | Positive real | 0.0 |
| `end_time` | End time for the mesh preview simulation | Positive real | - |
| `output_freq` | Number of timesteps between each generated output file | Positive integer | - |
| `dt` | Constant time step size used for the mesh preview | Positive real | - |


@attention The `"mesh_preview"` feature is strictly a pre-processing step. Once the preview completes, whether successfully or due to a failure, the Neko run will terminate and output a corresponding success or failure message. To proceed with the actual fluid simulation, the user **must** set `"mesh_preview.enabled: false"` in order for the actual simulation to run. Additionally, if the solver detects an inverted mesh element during the preview phase, it will save the exact time step at which the Jacobian becomes negative to assist with debugging.

@note The `"mesh_preview"` feature uses the same time integration order as defined in `"case.numerics.time_order"`. However, the `start_time`, `end_time`, constant `dt`, and `output_freq` parameters must be explicitly defined within the `"mesh_preview"` block itself.

@attention When the ALE module is `enabled`, Neko will **always** save the mesh in the `.fld` files at every output step.

#### Bodies

The `"bodies"` block defines an array of objects, where each object represents an individually controlled ALE body. Each body must map to corresponding physical boundary zone using `"bodies.zone_indices"`.

Each individual body object accepts the following general keywords and base kinematics:

| Name | Description | Admissible values | Default value |
| :--- | :--- | :--- | :--- |
| `name` | The name identifier of the body | String | `"body_<body_ID>"`  |
| `zone_indices` | The physical boundary zones associated with this body | Array of positive integers | -  |
| `oscillation` | Sub-object defining the translational oscillation kinematics | JSON object | - |
| `rotation` | Sub-object defining the rotational kinematics applied to the body | JSON object | - |
| `pivot` | Sub-object defining the center point for rotational kinematics | JSON object | - |
| `stiff_geom` | Sub-object defining the mesh stiffness region | JSON object | - |

@attention The body_ID for ALE bodies is defined based on the order in which they are added to the `"bodies"` array, not based on their `"zone_indices"`.

@note If multiple moving `no_slip` zone IDs are assigned to `"zone_indices"` of a single ALE body, the code will treat all those boundaries as a unified rigid body.

##### Oscillation

The `"oscillation"` sub-object defines the harmonic translational motion of the body. It takes the following mandatory keywords:

| Name | Description | Admissible values | Default value |
| :--- | :--- | :--- | :--- |
| `oscillation.amplitude` | Amplitude of translational oscillation in \f$ [x, y, z] \f$ | Array of 3 reals | - |
| `oscillation.frequency` | Frequency of translational oscillation in \f$ [x, y, z] \f$ | Array of 3 reals | - |

@warning If the `"oscillation"` block is included in the case file, both `"amplitude"` and `"frequency"` become **mandatory**.

When translational oscillation is configured, the displacement \f$ x_i(t) \f$ and velocity \f$ v_i(t) \f$ of the body follow a simple harmonic motion for each active directional component \f$ i \in \{x, y, z\} \f$:

\f{eqnarray*}{
  x_i(t) &=& A_i \sin(2\pi f_i t), \\
  v_i(t) &=& 2\pi f_i A_i \cos(2\pi f_i t),
\f}

where \f$ A_i \f$ is the `oscillation.amplitude` and \f$ f_i \f$ is the `oscillation.frequency` for that specific direction.


##### Rotation

If the body undergoes rotational motion, the `"rotation"` sub-object can be configured. Depending on the `rotation.type`, different parameters become applicable:

| Name | Description | Admissible values | Default value |
| :--- | :--- | :--- | :--- |
| `rotation.type` | The type of rotational kinematics applied | `"harmonic"`, `"ramp"`, `"smooth_step"` | - |
| `rotation.amplitude_deg` | Rotational amplitude in **degrees** <i>(only for </i>`harmonic`<i>)</i> | Array of 3 reals | - |
| `rotation.frequency` | Rotational frequency <i>(only for </i>`harmonic`<i>)</i> | Array of 3 reals | - |
| `rotation.ramp_omega0` | Target angular velocity <i>(only for </i>`ramp`<i>)</i> | Array of 3 reals | - |
| `rotation.ramp_t0` | Time constant for the ramp <i>(only for </i>`ramp`<i>)</i> | Array of 3 reals | - |
| `rotation.axis` | Axis of rotation <i>(only for </i>`smooth_step`<i>)</i> | `1` (x), `2` (y), `3` (z) | `3` |
| `rotation.target_angle_deg`| Target rotation angle in **degrees** <i>(only for </i>`smooth_step`<i>)</i> | Real | - |
| `rotation.step_control_times`| Control times \f$ [t_0, t_1, t_2, t_3] \f$ <i>(only for </i>`smooth_step`<i>)</i>| Array of 4 reals | - |

@warning If the `"rotation"` block is included in the case file, a valid `"pivot"` block to specify the center of rotation must be defined. The `"pivot"` object is explained [here](#case-file_fluid-ale-pivot). Additionally, the specific parameters corresponding to the chosen `rotation.type` become **mandatory**.

@attention Positive rotation is defined counter-clockwise in a right-handed coordinate system.

The angular velocity vector \f$ \mathbf{\omega}(t) \f$ and angular position vector \f$ \mathbf{\theta}(t) \f$ are computed based on the selected `rotation.type` as follows:

**1. Harmonic Rotation** (`"harmonic"`)
Applies a simple harmonic oscillation to the angular velocity and position for each axis \f$ i \in \{x, y, z\} \f$:

\f{eqnarray*}{
 \omega_i(t) &=& A_{rad,i} (2\pi f_i) \cos(2\pi f_i t), \\
 \theta_i(t) &=& A_{rad,i} \sin(2\pi f_i t),
\f}

where \f$ A_{rad,i} \f$ is the `amplitude_deg` converted to radians, and \f$ f_i \f$ is the `frequency`.

**2. Ramp Rotation** (`"ramp"`)
Gradually ramps up the angular velocity to a target value for each axis \f$ i \in \{x, y, z\} \f$:

\f{eqnarray*}{
 \omega_i(t) &=& \Omega_{0,i} \left( 1 - \exp\left(-4.6 \frac{t}{t_{0,i}}\right) \right), \\
 \theta_i(t) &=& \Omega_{0,i} \left[ t - \frac{t_{0,i}}{4.6} \left( 1 - \exp\left(-4.6 \frac{t}{t_{0,i}}\right) \right) \right],
\f}

where \f$ \Omega_{0,i} \f$ is the target angular velocity (`ramp_omega0`) and \f$ t_{0,i}\f$  is the time constant (`ramp_t0`). The factor of 4.6 ensures the velocity reaches approximately 99% of its target at \f$ t = t_{0,i} \f$. At times beyond the ramp parameter (\f$ t \gt t_{0,i} \f$), the body achieves a steady rotation rate \f$ \Omega_{0,i} \f$.

**3. Smooth Step Rotation** (`"smooth_step"`)
Applies a smooth rotation around a single specified `axis` using a derivative step function, \f$ \text{dstep}(\tau) \f$, to reach a target angle and return to zero. The motion along the chosen axis is divided into four phases defined by `step_control_times` \f$ [t_0, t_1, t_2, t_3] \f$, prescribing both the angular velocity \f$ \omega(t) \f$ and angular position \f$ \theta(t) \f$:

- **Rise Phase** (\f$ t_0 \le t < t_1 \f$):

    \f{eqnarray*}{
     \omega(t) &=& \frac{\theta_{rad}}{t_1 - t_0} \text{dstep}\left(\frac{t - t_0}{t_1 - t_0}\right), \\
     \theta(t) &=& \theta_{rad} \, S\left(\frac{t - t_0}{t_1 - t_0}\right)
    \f}

- **Hold Phase** (\f$ t_1 \le t < t_2 \f$):

    \f{eqnarray*}{
     \omega(t) &=& 0, \\
     \theta(t) &=& \theta_{rad}
    \f}

- **Fall Phase** (\f$ t_2 \le t < t_3 \f$):

    \f{eqnarray*}{
     \omega(t) &=& -\frac{\theta_{rad}}{t_3 - t_2} \text{dstep}\left(\frac{t - t_2}{t_3 - t_2}\right), \\
     \theta(t) &=& \theta_{rad} \left[ 1 - S\left(\frac{t - t_2}{t_3 - t_2}\right) \right]
    \f}

- **Rest Phase** (\f$ t < t_0\f$  or \f$ t \ge t_3 \f$):
    \f{eqnarray*}{
     \omega(t) &= 0, \\
     \theta(t) &= 0,
    \f}
where \f$ \theta_{rad} \f$ is the `target_angle_deg` converted to radians.

The base smooth step function \f$ S(\tau) \f$ and its analytical derivative \f$ \text{dstep}(\tau) \f$ are defined for \f$ \tau \in (0, 1) \f$ as:

\f{eqnarray*}{
 S(\tau) &=& \frac{1}{1 + \exp(g(\tau))}, \quad \text{where} \quad g(\tau) = \frac{1}{\tau - 1} + \frac{1}{\tau}, \\
 \text{dstep}(\tau) &=& -S(\tau)(1 - S(\tau))g'(\tau), \quad \text{where} \quad g'(\tau) = -\frac{1}{(\tau - 1)^2} - \frac{1}{\tau^2}.
\f}

For bounds where \f$ \tau \le 0 \f$, \f$ S(\tau) = 0 \f$ and \f$ \text{dstep}(\tau) = 0 \f$.
For bounds where \f$ \tau \ge 1 \f$, \f$ S(\tau) = 1 \f$ and \f$ \text{dstep}(\tau) = 0 \f$.

@attention Within the case file, both translational oscillation and rotational motion can be applied simultaneously to a body. However, only one rotation mode can be active at a time.

@note A custom motion logic can always be applied to a body via custom `user_ale_rigid_kinematics` or `user_ale_mesh_velocity` subroutines in the user file (see [here](#user-file_ale-rigid_motion) for more information). The example `"Double_ocyl_cylinder"` shows an example for this usage. More details regarding the supported ALE interfaces and subroutines can be found [here](#user-file_ale).

@attention There are several ways to calculate torque on a moving ALE body. The user is referred to [this part of the documentation](#simcomp_force_torque) for more information.

##### Pivot {#case-file_fluid-ale-pivot}

The `"pivot"` sub-object defines the center point around which the body rotates.

| Name | Description | Admissible values | Default value |
| :--- | :--- | :--- | :--- |
| `pivot.type` | Type of pivot definition | `"relative"` | `"relative"` |
| `pivot.value` | The spatial coordinates of the rotation center \f$ [x, y, z] \f$ | Array of 3 reals | - |


@note The rotation center (i.e., the `"pivot"`) moves rigidly with the body. This means that if a body undergoes both translational oscillation and rotation, `"pivot.value"` defines the initial position of the rotation center. Throughout the simulation, the location of the pivot point is numerically updated using the translational velocity of the body, even if a custom rigid motion is applied using a `user_ale_rigid_kinematics` subroutine (see [here](#user-file_ale-rigid_motion) for more information).

##### Mesh Stiffness {#case-file_fluid_ale_stiff_geom}

The `"stiff_geom"` sub-object defines a local geometric region where the mesh is kept highly rigid to maintain the original mesh quality in regions of interest, e.g., within the boundary layer.

The global mesh stiffness field \f$ h(\mathbf{x}) \f$ is constructed by taking a base stiffness of \f$ 1.0 \f$ and adding the maximum local stiffness contribution from all registered bodies:

\f{eqnarray*}{
 h(\mathbf{x}) = 1.0 + \max_{b \in \text{bodies}} (\text{Stiffness}_b(\mathbf{x})).
\f}

@attention A valid `"stiff_geom"` definition is **mandatory** for every registered body.

| Name | Description | Admissible values | Default value |
| :--- | :--- | :--- | :--- |
| `stiff_geom.type` | The shape of the stiffness region | `"cylinder"`, `"sphere"`, `"cheap_dist"` | -  |
| `stiff_geom.decay_profile`| How the stiffness decays away from the body | `"gaussian"`, `"tanh"` | -  |
| `stiff_geom.gain` | The gain multiplier for the stiffness field | Positive real | - |
| `stiff_geom.cutoff_coef` | Controls the steepness of the spatial decay | Positive real | `9.0` (gaussian), `3.5` (tanh) |
| `stiff_geom.center` | Center coordinates \f$ [x, y, z] \f$ <i>(only for </i>`cylinder`<i>, </i>`sphere`<i>)</i>| Array of 3 reals | - |
| `stiff_geom.radius` | Radius of the geometry <i>(only for </i>`cylinder`<i>, </i>`sphere`<i>)</i> | Positive real | - |
| `stiff_geom.stiff_dist` | Distance to maintain stiffness <i>(only for </i>`cheap_dist`<i>)</i> | Positive real | - |

###### Local stiffness
The local stiffness contribution from a body, \f$ \text{Stiffness}_b(\mathbf{x}) \f$, is evaluated based on the chosen `decay_profile`. These profiles depend on a raw distance \f$ r \f$ (detailed in the next section) and a characteristic decay length \f$ d \f$, which is defined as:

- \f$ d \f$ = `"radius"` if `"type"` is `"cylinder"` or `"sphere"`.
- \f$ d \f$ = `"stiff_dist"` if `"type"` is `"cheap_dist"`.

Based on the `decay_profile`, the stiffness is calculated as follows:

**1. Gaussian Profile** (`"gaussian"`)
Applies an exponential decay based on the squared normalized distance:

\f{eqnarray*}{
 \text{Stiffness}_b(\mathbf{x}) = \text{gain} \cdot \exp\left( -\left(\frac{r}{d}\right)^2 \cdot \text{cutoff_coef} \right).
\f}

**2. Tanh Profile** (`"tanh"`)
Applies a smooth hyperbolic tangent transition:

\f{eqnarray*}{
 \text{Stiffness}_b(\mathbf{x}) = \text{gain} \cdot \left[ 1 - \tanh\left( \frac{r}{d} \cdot \text{cutoff_coef} \right) \right].
\f}

###### Distance Calculation
For a given coordinate \f$ \mathbf{x} = (x, y, z) \f$, the raw distance \f$ r \f$ used in the stiffness formulas above is calculated based on the selected `stiff_geom.type` and the body's stiffness center \f$ C = (c_x, c_y, c_z) \f$:

* **Sphere** (`"sphere"`):

  \f{eqnarray*}{
   r = \sqrt{(x - c_x)^2 + (y - c_y)^2 + (z - c_z)^2}.
  \f}
* **Cylinder** (`"cylinder"`):
  Calculates the distance to the Z-axis passing through the center \f$ (c_x, c_y) \f$.

  \f{eqnarray*}{
   r = \sqrt{(x - c_x)^2 + (y - c_y)^2}.
  \f}

* **Wall Distance** (`"cheap_dist"`):
  \f$ r \f$ is assigned from a precomputed pseudo distance field based on the boundary `zone_indices`.

@attention Within the region defined by `radius` (from the center) or `stiff_dist` (from the boundary), the mesh stiffness is at its highest. If the `gain` parameter is set large enough, the mesh within this region moves rigidly with the body, preserving its original element quality without deformation. Users are encouraged to check the `ocyl_cylinder3D`, `ocyl_ellipse3D`, and `Double_ocyl_cylinder` examples to get a better idea of how these parameters are configured in practice.

@note Setting a very large value for `gain` (e.g., `1.0e6`) is recommended if the mesh immediately surrounding the body is intended to be fully rigid (i.e., \f$ \phi_i \approx 1 \f$). Conversely, if two moving objects are in close proximity, the `gain` should be kept low enough to ensure the mesh in the gap region remains soft and deformable. Visualizing the generated `phi_total0.f00000` file, and checking the mesh quality using `mesh_preview` are highly recommended.

#### Restarting ALE simulations

Neko supports checkpointing and restarting for ALE simulations. No additional parameters need to be set apart from the usual configuration for saving `.chkp` files.

@attention A `.chkp` file generated from a standard static simulation (i.e., `"ale.enabled": false`) cannot be used to restart an ALE simulation. However, if you run a static simulation to establish a base flow, that output field can be loaded as an `initial_condition` for a subsequent ALE simulation. In this case, saving the file in `double precision` is recommended.

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
* `projection_reorthogonalize_basis`, logical flag to update and re-orthogonalize
   the projection basis at each time. This option works only for pressure projection.
* `monitor`, monitoring of residuals. If set to true, the residuals will be
  printed for each iteration.

In addition to the above settings, the solvers can be configured with strict
convergence criteria. This is done by setting the
`case.fluid.strict_convergence` keyword to `true`. This will force the solver to
converge to the specified tolerance within the specified number of iterations.
If the solver does not converge, the simulation will be terminated.
This can in some situations cause issues if the initial condition is far from a
valid solution. Therefore a user can allow an initial stabilization phase by
setting the `case.fluid.allow_stabilization` keyword to `true`. In this case,
the strict convergence will be ignored untill all components of the velocity
field converge within the specified tolerance. After this initial stabilization
phase, strict convergence will be enforced for the rest of the simulation.

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
| `coarse_grid.cheby_degree`   | Degree of the Chebyshev based AMG smoother                                              | An integer              | 4             |

For `phmg`, the following keywords are used:

| Name                       | Description                                                             | Admissible values             | Default value |
| -------------------------- | ----------------------------------------------------------------------- | ----------------------------- | ------------- |
| `pcoarsening_schedule`     | P-multigrid coarsening schedule (polynomial order, high to low)         | Array of integers             | `[3, 1]`      |
| `smoother_iterations`      | Number of smoother iterations in the p-multigrid parts                  | An integer                    | 3             |
| `smoother_cheby_acc`       | Type of Chebyshev acceleration                                          | `none`, `jacobi` or `schwarz` | `jacobi`      |
| `coarse_grid.levels`       | Number of AMG levels to construct (only valid for `solver` type `tamg`) | An integer                    | 3             |
| `coarse_grid.iterations`   | Number of linear solver iterations for coarse grid solver               | An integer                    | 1             |
| `coarse_grid.cheby_degree` | Degree of the Chebyshev based AMG smoother                              | An integer                    | 4             |


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
* `log`, whether to print the flow-rate forcing log message each time the
  forcing is adjusted. Defaults to `true`.


### Full parameter table
All the parameters are summarized in the table below. This includes all the
subobjects discussed above, as well as keyword parameters that can be described
concisely directly in the table.

| Name                                               | Description                                                                                       | Admissible values                                           | Default value |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------- | ----------------------------------------------------------- | ------------- |
| `scheme`                                           | The fluid solve type.                                                                             | `pnpn`                                                      | -             |
| `name`                                             | The name associated to the fluid solver.                                                          | String                                                      | `fluid`       |
| `Re`                                               | The Reynolds number.                                                                              | Positive real                                               | -             |
| `rho`                                              | The density of the fluid.                                                                         | Positive real                                               | -             |
| `mu`                                               | The dynamic viscosity of the fluid.                                                               | Positive real                                               | -             |
| `nut_field`                                        | The name of the turbulent viscosity field.                                                        | String                                                      | -             |
| `output_control`                                   | Defines the interpretation of `output_value` to define the frequency of writing checkpoint files. | `nsamples`, `simulationtime`, `tsteps`, `never`             | -             |
| `output_value`                                     | The frequency of sampling in terms of `output_control`.                                           | Positive real or integer                                    | -             |
| `output_mesh_in_all_files`                         | Indicates if the mesh should be written in every output fld file.                                 | `true` or `false`                                           | `false`       |
| `output_filename`                                  | The output filename.                                                                              | String                                                      | `field`       |
| `output_subdivide`                                 | Whether to subdivide spectral elements into linear sub-cells for VTKHDF output.                   | `true` or `false`                                           | `false`       |
| `inflow_condition.type`                            | Velocity inflow condition type.                                                                   | `user`, `uniform`, `blasius`                                | -             |
| `inflow_condition.value`                           | Value of the inflow velocity.                                                                     | Vector of 3 reals                                           | -             |
| `initial_condition.type`                           | Initial condition type.                                                                           | `user`, `uniform`, `blasius`, `field`                       | -             |
| `initial_condition.value`                          | Value of the velocity initial condition.                                                          | Vector of 3 reals                                           | -             |
| `initial_condition.file_name`                      | If `"type" = "field"`, the path to the field file to read from.                                   | String ending with `.fld`, `.chkp`, `.nek5000` or `f*****`. | -             |
| `initial_condition.sample_index`                   | If `"type" = "field"`, and file type is `fld` or `nek5000`, the index of the file to sampled.     | Positive integer.                                           | -1            |
| `initial_condition.previous_mesh`                  | If `"type" = "field"`, and file type is `chkp`, the previous mesh from which to interpolate.      | String ending with `.nmsh`.                                 | -             |
| `initial_condition.tolerance`                      | If `"type" = "field"`, and file type is `chkp`, tolerance to use for mesh interpolation.          | Positive real.                                              | 1e-6          |
| `blasius.delta`                                    | Boundary layer thickness in the Blasius profile.                                                  | Positive real                                               | -             |
| `blasius.freestream_velocity`                      | Free-stream velocity in the Blasius profile.                                                      | Vector of 3 reals                                           | -             |
| `blasius.approximation`                            | Numerical approximation of the Blasius profile.                                                   | `linear`, `quadratic`, `cubic`, `quartic`, `sin`, `tanh`    | -             |
| `shear_stress.value`                               | The shear stress vector value for `sh` boundaries                                                 | Vector of 3 reals                                           | `[0, 0, 0]`   |
| `wall_modelling.type`                              | The wall model type for `wm` boundaries. See documentation for additional config parameters.      | `rough_log_law`, `spalding`                                 | -             |
| `source_terms`                                     | Array of JSON objects, defining additional source terms.                                          | See list of source terms above                              | -             |
| `gradient_jump_penalty`                            | Array of JSON objects, defining additional gradient jump penalty.                                 | See list of gradient jump penalty above                     | -             |
| `boundary_types`                                   | Boundary types/conditions labels.                                                                 | Array of strings                                            | -             |
| `velocity_solver.type`                             | Linear solver for the momentum equation.                                                          | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres`                 | -             |
| `velocity_solver.preconditioner.type`              | Linear solver preconditioner for the momentum equation.                                           | `ident`, `hsmg`, `jacobi`                                   | -             |
| `velocity_solver.absolute_tolerance`               | Linear solver convergence criterion for the momentum equation.                                    | Positive real                                               | -             |
| `velocity_solver.maxiter`                          | Linear solver max iteration count for the momentum equation.                                      | Positive real                                               | 800           |
| `velocity_solver.projection_space_size`            | Projection space size for the momentum equation.                                                  | Positive integer                                            | 0             |
| `velocity_solver.projection_hold_steps`            | Holding steps of the projection for the momentum equation.                                        | Positive integer                                            | 5             |
| `velocity_solver.monitor`                          | Monitor residuals in the linear solver for the momentum equation.                                 | `true` or `false`                                           | `false`       |
| `pressure_solver.type`                             | Linear solver for the pressure equation.                                                          | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres`                 | -             |
| `pressure_solver.preconditioner.type`              | Linear solver preconditioner for the pressure equation.                                           | `ident`, `hsmg`, `jacobi`                                   | -             |
| `pressure_solver.absolute_tolerance`               | Linear solver convergence criterion for the pressure equation.                                    | Positive real                                               | -             |
| `pressure_solver.maxiter`                          | Linear solver max iteration count for the pressure equation.                                      | Positive real                                               | 800           |
| `pressure_solver.projection_space_size`            | Projection space size for the pressure equation.                                                  | Positive integer                                            | 0             |
| `pressure_solver.projection_hold_steps`            | Holding steps of the projection for the pressure equation.                                        | Positive integer                                            | 5             |
| `pressure_solver.projection_reorthogonalize_basis` | Whether to enable pressure projection basis reorthogonalization.                                  | `true` or `false`                                           | `false`       |
| `pressure_solver.monitor`                          | Monitor residuals in the linear solver for the pressure equation.                                 | `true` or `false`                                           | `false`       |
| `flow_rate_force.direction`                        | Direction of the forced flow.                                                                     | 0, 1, 2                                                     | -             |
| `flow_rate_force.value`                            | Bulk velocity or volumetric flow rate.                                                            | Positive real                                               | -             |
| `flow_rate_force.use_averaged_flow`                | Whether bulk velocity or volumetric flow rate is given by the `value` parameter.                  | `true` or `false`                                           | -             |
| `flow_rate_force.log`                              | Whether to print the flow-rate forcing log message during the volume-flow adjustment.             | `true` or `false`                                           | `true`        |
| `freeze`                                           | Whether to fix the velocity field at initial conditions.                                          | `true` or `false`                                           | `false`       |
| `strict_convergence`                               | Whether to enforce strict convergence in the linear solvers.                                      | `true` or `false`                                           | `false`       |
| `allow_stabilization`                              | Whether to allow an initial stabilization phase before enforcing strict convergence.              | `true` or `false`                                           | `false`       |
| `advection`                                        | Whether to compute the advection term.                                                            | `true` or `false`                                           | `true`        |
| `full_stress_formulation`                          | Whether to use the full form of the visous stress tensor term.                                    | `true` or `false`                                           | `false`       |

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

Different from the setup in the fluid, turbulence modelling is enabled by setting the `alphat` json entry.

### Turbulence modelling

The user could choose to either relate the eddy diffusivity field to the eddy viscosity
field in the Fluid, or model the eddy diffusivity field by some particular SGS models,
by setting up the `nut_dependency` entry.
If the eddy diffusivity field is associated to the eddy viscosity field by a coefficient
`Pr_t`, the eddy viscosity values will be divided by it to produce eddy diffusivity.
And the corresponding setting could be done by the following:

```json
"alphat":{
    "nut_dependency": true,
    "nut_field": "nut",
    "Pr_t": 0.7
},
```

Otherwise one could have some SGS models providing an eddy diffusivity field, and the
user could set it up by the following manner to include an eddy diffusivity field called
`temperature_alphat`:

```json
"alphat":{
    "nut_dependency": false,
    "alphat_field": "temperature_alphat"
},
```

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

| Name                           | Description                                                           | Admissible values                           | Default value |
| ------------------------------ | --------------------------------------------------------------------- | ------------------------------------------- | ------------- |
| `enabled`                      | Whether to enable the scalar computation.                             | `true` or `false`                           | `true`        |
| `name`                         | The name associated to the scalar solver.                             | String                                      | `scalar`      |
| `field_name`                   | The name of the solution in the field registry.                       | A string                                    | `s`           |
| `Pe`                           | The Peclet number.                                                    | Positive real                               | -             |
| `cp`                           | Specific heat capacity.                                               | Positive real                               | -             |
| `lambda`                       | Thermal conductivity.                                                 | Positive real                               | -             |
| `alphat.nut_dependency`        | Whether the eddy diffusivity depends on the eddy kinematic viscosity. | `true` or `false`                           | -             |
| `alphat.alphat_field`          | Name of the turbulent diffusivity field.                              | String                                      | Empty string  |
| `alphat.nut_field`             | Name of the turbulent kinematic viscosity field.                      | String                                      | Empty string  |
| `alphat.Pr_t`                  | Turbulent Prandtl number                                              | Positive real                               | -             |
| `boundary_types`               | Boundary types/conditions labels.                                     | Array of strings                            | -             |
| `initial_condition.type`       | Initial condition type.                                               | `user`, `uniform`, `point_zone`             | -             |
| `initial_condition.value`      | Value of the velocity initial condition.                              | Real                                        | -             |
| `source_terms`                 | Array of JSON objects, defining additional source terms.              | See list of source terms above              | -             |
| `gradient_jump_penalty`        | Array of JSON objects, defining additional gradient jump penalty.     | See list of gradient jump penalty above     | -             |
| `advection`                    | Whether to compute the advetion term.                                 | `true` or `false`                           | `true`        |
| `solver.type`                  | Linear solver for scalar equation.                                    | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres` | -             |
| `solver.preconditioner.type`   | Linear solver preconditioner for the momentum equation.               | `ident`, `hsmg`, `jacobi`                   | -             |
| `solver.absolute_tolerance`    | Linear solver convergence criterion for the momentum equation.        | Positive real                               | -             |
| `solver.maxiter`               | Linear solver max iteration count for the momentum equation.          | Positive real                               | 800           |
| `solver.projection_space_size` | Projection space size for the scalar equation.                        | Positive integer                            | 0             |
| `solver.projection_hold_steps` | Holding steps of the projection for the scalar equation.              | Positive integer                            | 5             |


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
