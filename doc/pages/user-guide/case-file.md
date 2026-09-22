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

@note We use the [json-fortran](https://github.com/jacobwilliams/json-fortran)
module to handle the JSON files. In addition to standard JSON, it permits
using `/`, `!` or `#` to make line comments and tolerates trailing commas
before a closing `}` or `]`.

A JSON Schema for the case file is available in `doc/schemas/`. The helper
`contrib/validate_case_schema.py` can be used to validate case files against
that schema. This is a convenient way to quickly check that your case file is 
valid without having to run a simulation.

@warning Neko distinguishes JSON integers such as `10` from JSON
reals such as `10.0`. Make sure that you use the correct number type for the
parameter you are setting.

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
        "point_zones" : [],
        "runtime_statistics": {}
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

The `_value` parameter defines the output frequency, but both its JSON type and
interpretation depend on the choice of `_control`. Neko requires an integer
literal for `tsteps` and `nsamples`, and a real literal for `simulationtime`.
For example, use `10` for ten time steps but `10.0` for ten simulation-time
units. The following options are possible.
1. `simulationtime`, then `_value` is the time interval between the outputs.
2. `tsteps`, then `_value` is the number of time steps between the outputs.
3. `nsamples`, then `_value` is the total number of outputs that will be
   performed in the course of the simulation.
4. `never`, then `_value` is ignored and output is never performed.

The times at which an output is written are a property of the case, not of
the run. With `simulationtime` they are the whole multiples of
\f$ \Delta t_{out} \f$, so a case asking for an output every 1.0 writes at
11.0, 12.0, 13.0 whether it started from 0, from 10.5, or from a checkpoint
in between. With `nsamples` the samples divide the interval over which the
output is active instead, \f$ t_k = t_{start} + k (t_{end} - t_{start}) / n
\f$, where \f$ t_{start} \f$ is the start of the simulation or, for the
statistics, the start of the averaging, and the last sample falls on
`end_time`. Only the times inside the interval where the output is active
are written, and since the schedule does not depend on when the run was
started, restarting neither skips nor repeats an output and the file numbering
carries on.

Each of those times is written at the first step to reach it, at most one file
per output per step, so a time step longer than the interval writes once and
moves on rather than building up a backlog. The write scheduled for `end_time`
still happens when the last step overshoots it.

The initial state is written in addition to the scheduled times, unless
`output_at_start` is `false`. Checkpoints and statistics never write it: a
checkpoint holds nothing the initial condition does not, and an average over
an interval of zero length is empty.

@note A statistics file covers the interval between two writes, and the first
one covers only the interval from the start of the averaging to the first
scheduled time, which is shorter than the rest unless the two coincide. One
more reason to use the weighted averages when post-processing.


## The case object

This object is mostly used as a high-level container for all the other objects,
but also defines several parameters that pertain to the simulation as a whole.

| Name                  | Description                                                                                           | Admissible values                               | Default value |
| --------------------- | ----------------------------------------------------------------------------------------------------- | ----------------------------------------------- | ------------- |
| `no_defaults`         | Prevents filling in default values to case file entries.                                              | `true` or `false`                               | `false`       |
| `mesh_file`           | The name of the mesh file.                                                                            | Strings ending with `.nmsh`                     | -             |
| `output_boundary`     | Whether to write a `bdry0.f0000` file with boundary labels. Can be used to check boundary conditions. | `true` or `false`                               | `false`       |
| `output_directory`    | Folder for redirecting solver output. Note that the folder has to exist!                              | Path to an existing directory                   | `.`           |
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
| `job_timelimit`       | The maximum wall clock duration of the simulation.                                                    | String formatted as [[[DD-]HH:]MM:]SS           | No limit      |
| `output_at_end`       | Whether to always write all enabled output at the end of the run.                                     | `true` or `false`                               | `true`        |
| `output_at_start`     | Whether to write the initial state of the simulation. Checkpoints and statistics never do.            | `true` or `false`                               | `true`        |

Some additional practical comments are provided regarding the output triggered
by `job_timelimit` and `output_at_end` keywords.

If `output_at_end` is set to `true`, a write is performed after the execution
of the simulation time-loop is finished. This triggers most outputs, like the
fluid solvers, the checkpoint, etc. An output that has already been written at
the last time step of the run is not written a second time, so case settings
under which an output lands on the last step anyway (e.g. `end_time: 5`,
`checkpoint_value: 5`, `checkpoint_control: simulationtime`) produce one file
at the end of the run and not two. `output_at_end` does override a `never`
control, which is the way to ask for an output that is written once, at the
end of the run, and never in between.

The file numbering follows the schedule: the k-th file is the k-th scheduled
write, and a restart resumes where the schedule says the run has got to. A run
that repeats an interval it has already covered, or that changes the output
frequency at the restart, therefore writes over the files of the earlier run,
and Neko warns when the first file it is about to write already exists.

The purpose of `job_timelimit` is to gracefully stop the simulation in a typical
supercomputer environment, where your runtime is limited. When Neko detects that
the time of the run exceeds the `job_timelimit`, it exits the time-loop. At this
point, if one sets `output_at_end` to `true`, this will trigger a write as per
usual. However, if `output_at_end` is `false`, Neko will still write a special
checkpoint file, with the filename called `joblimit#####.chkp`. This is done so
that the user is at least provided a restart file, and none of the computer time
spent on the simulation is wasted. Generally, however, it is recommended to
have `output_at_end` set to `true` in tandem with `job_timelimit`, so that what
exactly gets written is controlled by the case file settings. Note that the time
format is flexible. For example, `1-01:00:00` and `25:00:00` are both valid ways
to specify a 25-hour time limit.

### Constants
The `constants` array allows the user to define parameters that are global to
the case file, and can be referred to when setting the values of other
parameters. Two types of parameters can be defined: scalars and arrays. Each is
represented as a subobject inside the `constants` object and should containt two
entries: `name` and `value`. Here is an example:

```json
{
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
}
```

Other parameters in the case file that require a scalar or array entry, can
instead be defined as a string, pointing to the name of the corresponding
parameter in the `constants` object. As an example, recall that output frequency
is controlled by the keyword `output_value`. It is a plausible scenario that the
frequency is the same for multiple solvers, simulation components, etc. Assuming
a simulation with both [fluid](@ref case-file_fluid) and [scalar](@ref
case-file_scalar) solvers active, the following could be used.

```json
{
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

### Expressions {#case-file_expressions}
Some parameters accept a mathematical expression, written as a string, instead
of a number. The expression is evaluated in every point where the parameter is
needed, which makes it possible to prescribe a spatially varying profile
without writing a [user file](@ref user-file).

```json
{
  "constants":
  [
    { "name": "U_b", "value": 1.0 },
    { "name": "H", "value": 2.0 }
  ],
  "fluid":
  {
    "initial_condition":
    {
      "type": "expression",
      "value": ["6*U_b*y*(H - y)/H^2", "0", "0"]
    }
  }
}
```

The following can be used in an expression.

| Kind        | Available                                                                                                                     |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------- |
| Coordinates | `x`, `y`, `z`                                                                                                                 |
| Time        | `t`, `dt`, where the parameter is evaluated during the simulation                                                             |
| Constants   | `pi`, and any scalar declared under `case.constants`                                                                          |
| Operators   | `+`, `-`, `*`, `/`, `^` (also written `**`), parentheses                                                                      |
| Functions   | `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `exp`, `log`, `log10`, `sqrt`, `abs`, `erf`, `erfc`, `step` |
|             | `atan2(y, x)`, `min(a, b)`, `max(a, b)`, `mod(a, b)`                                                                          |

`step(a)` is `0` for `a < 0` and `1` otherwise. Precedence and associativity
follow Fortran, so `-x^2` is `-(x^2)` and `2^3^2` is `2^(3^2)`.

Constants are read when the expression is compiled, at setup, and multiplication
must always be written out, so `2*x` is valid whereas `2x` is not. An expression
that refers to an unknown name is reported as an error at setup, rather than
during the simulation. An expression that does not evaluate to a finite value
everywhere, typically because of a division by zero or the square root of a
negative number, is also reported as an error, at setup if it does not depend on
time and otherwise every time it is evaluated.

### Time control
The `time` object is used to define the time-stepping of the simulation,
including the time-step size, the start and end time, and the variables related
to the variable time-stepping algorithm. For the variable timestep, one can
specify a `timestep`, such that the first timestep will be assigned to the
smallest of `timestep` and the value calculated from the target CFL number.

| Name                       | Description                                                                                 | Admissible values                 | Default value |
| -------------------------- | ------------------------------------------------------------------------------------------- | --------------------------------- | ------------- |
| `start_time`               | Start time at which the simulation is initiated.                                            | Non-negative reals                | `0.0`         |
| `end_time`                 | Final time after which the simulation is stopped.                                           | Positive reals                    | -             |
| `timestep`                 | Time-step size                                                                              | Positive reals                    | -             |
| `variable_timestep`        | Whether to use variable dt                                                                  | `true` or `false`                 | `false`       |
| `max_timestep`             | Maximum time-step size when variable time step is activated                                 | Positive reals                    | `huge`        |
| `min_timestep`             | Minimum time-step size when variable time step is activated                                 | Non-negative reals                | `0.0`         |
| `target_cfl`               | The desired CFL number                                                                      | Positive real                     | `0.4`         |
| `max_update_frequency`     | Minimum number of time steps between two time-step updates triggered by CFL deviation       | Non-negative integer              | `0`           |
| `min_update_frequency`     | Maximum number of time steps before forcing a time-step update                              | Non-negative integer              | `huge`        |
| `max_dt_increase_factor`   | The maximum scaling factor to increase time step                                            | Positive real greater than `1`    | `1.2`         |
| `min_dt_decrease_factor`   | The minimum scaling factor to decrease time step                                            | Positive real less than `1`       | `0.5`         |
| `cfl_deviation_tolerance`  | The tolerance of the deviation from the target CFL number                                   | Positive real less than `1`       | `0.2`         |
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
| no_slip (stationary wall)       | 1   |
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
| `polynomial_order`           | The order of the polynomial basis.                                                                              | Integers >= 3 for `pnpn`; integers >= 1 for `compressible` (typically 5 to 9) | -                               |
| `time_order`                 | The order of the time integration scheme. Refer to the `time_scheme_controller` type documentation for details. | 1, 2, 3                    | -                               |
| `dealias`                    | Whether to apply dealiasing to advection terms.                                                                 | `true` or `false`          | `false`                         |
| `dealiased_polynomial order` | The polynomial order in the higher-order space used in the dealising.                                           | Integer                    | `3/2(polynomial_order + 1) - 1` |
| `oifs`                       | Whether to apply the Operator-Integration-Factor-Splitting (OIFS).                                              | `true` or `false`          | `false`                         |
| `oifs_target_cfl`            | The desired OIFS-CFL number. Requires variable_timestep = true in the time control object.                      | Positive real              | `1.9`                           |

## Solvers {#case-file_solvers}

The configuration of the fluid and scalar solvers is described on separate
pages:

- \subpage case-file_fluid
- \subpage case-file_fluid-compressible
- \subpage case-file_scalar

## Simulation components
Simulation components enable the user to perform various additional operations,
which are not strictly necessary to run the solver. An example could be
computing and output of additional fields, e.g. vorticity.

A more detailed description as well as a  full list of available components and
 their setup is provided in a [separate page of the manual](@ref simcomps).

## Point zones
Point zones enable the user to select GLL points in the computational domain
according to some geometric criterion. Three predefined geometric shapes are
selectable from the case file: boxes, spheres, and cylinders.

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
