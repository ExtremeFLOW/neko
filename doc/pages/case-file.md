# Case File {#case-file}

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

Name                 | Description                                                    | Admissable values                         | Default value 
-----                | -----------                                                    | -----------------                         | -------------
`mesh_file`          | The name of the mesh file.                                                                            | Strings ending with `.nmsh`               | -  
`output_boundary`    | Whether to write a `bdry0.f0000` file with boundary labels. Can be used to check boundary conditions. | `true` or `false`                         | `false`       
`output_directory`   | Folder for redirecting solver output. Note that the folder has to exist!                              | Path to an existing directory             | `.` 
`load_balancing`     | Whether to apply load balancing.                                                                      | `true` or `false`                         | `false` 
`output_partitions`  | Whether to write a `partitions.vtk` file with domain partitioning.                                    | `true` or `false`                         | `false` 
`output_checkpoints` | Whether to output checkpoints, i.e. restart files.                                                    | `true` or `false`                         | `false` 
`checkpoint_control` | Defines the interpretation of `checkpoint_value` to define the frequency of writing checkpoint files. | `nsamples`, `simulationtime`, `tsteps`, `never` | -  
`checkpoint_value` | The frequency of sampling in terms of `checkpoint_control`. | Positive real or integer | -  
`restart_file` | checkpoint to use for a restart from previous data | Strings ending with `.chkp` | -  
`time_step`          | Time-step size.                                                                                       | Positive reals                            | -  
`end_time`           | Final time at which the simulation is stopped.                                                        | Positive reals                            | -   
`job_timelimit`      | The maximum wall clock duration of the simulation.                                                    | String formatted as HH:MM:SS              | No limit 

## Numerics
Used to define the properties of the numerical discretization.

Name                 | Description                                                        | Admissable values                         | Default value  
----                 | -----------                                                        | -----------------                         | -------------
`polynomial_order`   | The oder of the polynomial basis.                           | Integers, typically 5 to 9               | -             |
`time_order`    | The order of the time integration scheme. Refer to the `time_scheme_controller` type documention for details. | 1,2, 3                         | -    
`dealias`    | Whether to apply dealiasing to advection terms. | `true` or `false`                         | `false`     
`dealiased_polynomial order`    | The polynomial order in the higher-order space used in the dealising. | Integer                         | `3/2(polynomial_order + 1) - 1`    

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


### Inflow boundary conditions
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

### Initial conditions
The object `initial_condition` is used to provide initial conditions.
It is mandatory.
Note that this currently pertains to both the fluid, but also scalars.
The means of prescribing the values are controlled via the `type` keyword:

1. `user`, the values are set inside the compiled user file. The only way to
   initialize scalars.
2. `uniform`, the value is a constant vector, looked up under the `value`
   keyword.
3. `blasius`, a Blasius profile is prescribed. Its properties are looked up
   in the `case.fluid.blasius` object, see below.

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

### Source terms
The `source_terms` object should be used to specify the source terms in the
momentum equation. The object is not mandatory, by default no forcing term is
present. Each source term, is itself a JSON object, so `source_terms` is just an
array of them. Note that with respect to the governing equations, the source
terms define \f$ f^u \f$, meaning that the values are then multiplied by the
density.

For each source, the `type` keyword defines the kind of forcing that will be
introduced. The following types are currently implemented.

1. `constant`, constant forcing. Strength defined by the `values` array with 3
   reals corresponding to the 3 components of the forcing. 
2. `user_pointwise`, the values are set inside the compiled user file, using the
   pointwise user file subroutine. Only works on CPUs!
2. `user_vector`, the values are set inside the compiled user file, using the 
   non-pointwise user file subroutine. Should be used when running on the GPU.
   

### Boundary types
The optional `boundary_types` keyword can be used to specify boundary conditions.
The reason for it being optional, is that some conditions can be specified
directly inside the mesh file.
In particular, this happens when Nek5000 `.re2` files are converted to `.nmsh`.
Periodic boundary conditions are *always* defined inside the mesh file.

The value of the keyword is an array of strings, with the following possible
values:
* `w`, a no-slip wall.
* `v`, a Dirichlet boundary.
* `sym`, a symmetry boundary.
* `on`, Dirichlet for the boundary-parallel velocity and homogeneous Neumann for
   the wall-normal. The wall-parallel velocity is defined by the initial
   condition. 
* `o`, outlet boundary. 
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

## Linear solver configuration
The mandatory `velocity_solver` and `pressure_solver` objects are used to 
configure the solvers for the momentum and pressure-Poisson equation.
The following keywords are used, with the corresponding options.

* `type`, solver type.
  - `cg`, a conjugate gradient solver.
  - `pipecg`, a pipelined conjugate gradient solver.
  - `bicgstab`, a bi-conjugate gradient stabilized solver.
  - `cacg`, a communication-avoiding conjugate gradient solver.
  - `gmres`, a GMRES solver. Typically used for pressure.
* `preconditioner`, preconditioner type.
  - `jacobi`, a Jacobi preconditioner. Typically used for velocity.
  - `hsmg`, a hybrid-Schwarz multigrid preconditioner. Typically used for pressure.
  - `ident`, an identity matrix (no preconditioner).
* `absolute_tolerance`, tolerance criterion for convergence.
* `max_iterations`, maximum number of iterations before giving up.
* `projection_space_size`, size of the vector space used for accelerating the
   solution procedure. If 0, then the projection space is not used.
   More important for the pressure equation.

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

Name                                    | Description                                                                      | Admissable values                                | Default value
----------------------------------------|----------------------------------------------------------------------------------|--------------------------------------------------|--------------
`scheme`                                | The fluid solve type.                                                            | `pnpn`                                           | -
`Re`                                    | The Reynolds number.                                                             | Positive real                                    | -
`rho`                                   | The density of the fluid.                                                        | Positive real                                    | -
`mu`                                    | The dynamic viscosity of the fluid.                                              | Positive real                                    | -
`output_control` | Defines the interpretation of `output_value` to define the frequency of writing checkpoint files. | `nsamples`, `simulationtime`, `tsteps`, `never` | -  
`output_value` | The frequency of sampling in terms of `output_control`. | Positive real or integer | -  
`inflow_condition.type`                 | Velocity inflow condition type.                                                  | `user`, `uniform`, `blasius`                     | -
`inflow_condition.value`                | Value of the inflow velocity.                                                    | Vector of 3 reals                                | -
`initial_condition.type`                | Initial condition type.                                                          | `user`, `uniform`, `blasius`                     | -
`initial_condition.value`               | Value of the velocity initial condition.                                         | Vector of 3 reals                                | -
`blasius.delta`                         | Boundary layer thickness in the Blasius profile.                                 | Positive real                                    | -
`blasius.freestream_velocity`           | Freestream velocity in the Blasius profile.                                      | Vector of 3 reals                                | -
`blasius.approximation`                 | Numerical approximation of the Blasius profile.                                  | `linear`, `quadratic`, `cubic`, `quartic`, `sin` | -
`source_term.type`                      | Source term in the momentum equation.                                            | `noforce`, `user`, `user_vector`                 | -
`boundary_types`                        | Boundary types/conditions labels.                                                | Array of strings                                 | -
`velocity_solver.type`                  | Linear solver for the momentum equation.                                         | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres`      | -
`velocity_solver.preconditioner`        | Linear solver preconditioner for the momentum equation.                          | `ident`, `hsmg`, `jacobi`                        | -
`velocity_solver.absolute_tolerance`    | Linear solver convergence criterion for the momentum equation.                   | Positive real                                    | -
`velocity_solver.maxiter`               | Linear solver max iteration count for the momentum equation.                     | Positive real                                    | 800
`velocity_solver.projection_space_size` | Projection space size for the momentum equation.                                 | Positive integer                                 | 0
`pressure_solver.type`                  | Linear solver for the momentum equation.                                         | `cg`, `pipecg`, `bicgstab`, `cacg`, `gmres`      | -
`pressure_solver.preconditioner`        | Linear solver preconditioner for the momentum equation.                          | `ident`, `hsmg`, `jacobi`                        | -
`pressure_solver.absolute_tolerance`    | Linear solver convergence criterion for the momentum equation.                   | Positive real                                    | -
`pressure_solver.maxiter`               | Linear solver max iteration count for the momentum equation.                     | Positive real                                    | 800
`pressure_solver.projection_space_size` | Projection space size for the momentum equation.                                 | Positive integer                                 | 0
`flow_rate_force.direction`             | Direction of the forced flow.                                                    | 0, 1, 2                                          | -
`flow_rate_force.value`                 | Bulk velocity or volumetric flow rate.                                           | Positive real                                    | -
`flow_rate_force.use_averaged_flow`     | Whether bulk velocity or volumetric flow rate is given by the `value` parameter. | `true` or `false`                                | -            
`freeze`                                | Whether to fix the velocity field at initial conditions.                          | `true` or `false`                                | `false`            

## Scalar
The scalar object allows to add a scalar transport equation to the solution.
The solution variable is called `s`, but saved as `temperature` in the fld
 files. 
Some properties of the object are inherited from `fluid`: the properties of the
linear solver, the value of the density, and the output
control.

The scalar equation requires defining additional material properties: the
specific heat capacity and thermal conductivity. These are provided as `cp` and
`lambda`. Similarly to the fluid, one can provide the Peclet number, `Pe`, as an
alternative. In this case, `cp` is set to 1 and `lambda` to the inverse of `Pe`. 

The boundary conditions for the scalar are specified through the
`boundary_types` keyword.
It is possible to directly specify a uniform value for a Dirichlet boundary.
The syntax is, e.g. `d=1`, to set the value to 1, see the Ryleigh-Benard
example case.

Note that the source term configuration for the scalar currently differs from
that of the fluid. This will be addressed in a future release. For now, the
`source_term` keyword should be used, set to either `noforce` (no forcing),
`user` (same as `user_pointwise` for the fluid), and `user_vector` (same as for
the fluid).

Name               | Description                               | Admissable values                | Default value
-------------------|-------------------------------------------|----------------------------------|--------------
`enabled`          | Whether to enable the scalar computation. | `true` or `false`                | `true`
`Pe`               | The Peclet number.                        | Positive real                    | -
`cp`               | Specific heat cpacity.                    | Positive real                    | -
`lambda`           | Thermal conductivity.                     | Positive real                    | -
`source_term.type` | Source term in the momentum equation.     | `noforce`, `user`, `user_vector` | -
`boundary_types`   | Boundary types/conditions labels.         | Array of strings                 | -

## Statistics

This object adds the collection of statistics for the flud fields.
For additional details on the workflow, see the corresponding page in the
user manual.

Name                | Description                                                          | Admissable values | Default value
--------------------|----------------------------------------------------------------------|-------------------|--------------
`enabled`           | Whether to enable the statistics computation.                        | `true` or `false` | `true`
`start_time`        | Time at which to start gathering statistics.                         | Positive real     | 0
`sampling_interval` | Interval, in timesteps, for sampling the flow fields for statistics. | Positive integer  | 10

## Simulation components
Simulation components enable the user to perform various additional operations,
which are not strictly necessary to run the solver. An example could be
computing and output of additional fields, e.g. vorticity.

A more detailed description as well as a  full list of available components and
 their setup is provided in a [separate page of the manual](simcomps.md).

## Point zones
Point zones enable the user to select GLL points in the computational domain according to some geometric criterion. Two predefined geometric shapes are selectable from the case file, boxes and spheres.

A point zone object defined in the case file can be retrieved from the point zone registry, `neko_point_zone_registry`, and can be used to perform any zone-specific operations (e.g. localized source term, probing...). User-specific point zones can also be added manually to the point zone registry from the user file.

A more detailed description as well as a  full list of available components and
 their setup is provided in a [separate page of the manual](point-zones.md).
