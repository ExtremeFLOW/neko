# Simulation components {#simcomps}

\tableofcontents

## What are simulation components?
Simulation components, or simcomps for short, encapsulate additional
functionalities that may be useful for certain cases but not necessary to run the
solver.
This can include computation and output of additional fields, in-situ
post-processing operations, data sampling, etc.

By design, simulation components can tap into every aspect of the simulation,
so they can be quite powerful.
As the code grows, we expect to add more and more simcomps to the code.

## Adding simulation components to the case
Each simcomp is defined as a single JSON object at are added to an array
of objects called `simulation_components`, which resides directly under the
`case` object.

All simcomps support a `name` keyword in their JSON object. When ommitted, a
default `name`, which as a rule coincided with its `type`. For completness,
these default names are provided for each simcomp in the documentation below.
However, each simcomp's name must be unique, so if two or more simcomps of the
same type are present in case, a unique `name` for each must be provided
manually.

## List of simulation components

The following is a list of simulation components that are currently available
in Neko. The list will be updated as new simcomps are added.

- Differential operators
  - Computation of the curl of a vector field \ref simcomp_curl
  - Computation of the gradient of a scalar field \ref simcomp_gradient
  - Computation of the weak gradient of a field \ref simcomp_weak_gradient
  - Computation of the derivative of a scalar field \ref simcomp_derivative
  - Computation of the divergence of a vector field \ref simcomp_divergence
- Computation of \f$ \lambda_2 \f$ \ref simcomp_lambda2
- Probing of fields at selected points \ref simcomp_probes
- Output of registered fields to a file \ref simcomp_field_writer
- Computation of forces and torque on a surface \ref simcomp_force_torque
- Computation of subgrid-scale (SGS) eddy viscosity via a SGS model \ref
  simcomp_les_model
- User defined components \ref user-file_simcomps
- Fluid statistics simcomp, "fluid_stats", for more details see the
  [statistics guide](@ref statistics-guide)
- Fluid SGS statistics simcomp, "fluid_sgs_stats", for more details see the
  [statistics guide](@ref statistics-guide)
- Scalar statistics simcomp, "scalar_stats", for more details see the
  [statistics guide](@ref statistics-guide). If a `field` is specified in the
  case file without `name` being specified, the field name is appended to the 
  default simcomp name as `scalar_stats_{field}`.
- Scalar SGS statistics simcomp, "scalar_sgs_stats", for more details see the
  [statistics guide](@ref statistics-guide)
- User statistics simcomp, "user_stats" \ref user_stats
- Computation of the spectral error indicator \ref simcomp_speri
- Streaming of data for in-situ field manipulation \ref simcomp_data_streamer
- Sub-sampling of fields by changing polynomial order and masking by point
  zones \ref simcomp_field_subsampler

## Controlling execution and file output
Each simulation component is, by default, executed once per time step to perform
associated computations and output. However, this can be modified by using the
`compute_control` and `compute_value` parameters for the computation and the
`output_control` and `output_value` for the output to disk. The parameters
for the `_control` values are the same as for the fluid and checkpointing.
Additionally, one can set `output_control` to `global` and `never`. The former
will sync the `output_` parameter to that of the fluid. Choosing `never` will
suppress output all together. If no parameters for the `output_` parameters are
provided, they are set to be the same as for `compute_`. In order to simplify
the configuration, the `compute_control` can be set to `fluid_output` to sync
the computation to the fluid output. All `_value` keywords can also be strings
pointing to an entry under `constants.scalars`.

For simcomps that compute 3D fields, the output can be either added to the main
`.fld` file, containing velocity and pressure, or saved to a separate file. For
the latter, the `output_filename` keyword should be provided. One can
additionally provide the `precision` keyword, which can be set to either
`single` or `double` to control the precision of the written data.

For example, in the `tgv` example case the `curl` component is executed
once per 50 time steps. The `output_` parameters are synced to that, and the
vorticity fields will be added to the main `.fld` file.
~~~~~~~~~~~~~~~{.json}
{
  "type": "curl",
  "name": "curl",
  "field_names": ["u", "v", "w"],
  "computed_field": "vorticity"
  "compute_control": "tsteps",
  "compute_value": 50
}
~~~~~~~~~~~~~~~

### Differential operators

There is a set of simcomps that allow one to apply various differential
operators on registered fields. They share common configuration traits. The
field or fields that the operators is applied to us controlled, respectively, by
the `field` or `feilds` keyword. The produced fields are also added to the
registry, and each operator provides a default name. However, it can be
overriden using the `computed_field` keyword. Operators that output a vector
field will register three fields, adding `_x`, `_y`, and `_z` to the base of the
name.

All of these simcomps also support saving the result to `.fld` files. The \ref
simcomp_field_writer simcomp is used for that under the hood, so the associated
JSON keywords can be found in its documentation (`output_filename`,
`precision`).

#### derivative {#simcomp_derivative}
Computes the derivative of field along a chosen direction (x, y, or z). The
field to derivate is controlled by the `field` keyword and the direction by the
`direction` keyword. The simcomp will, by default, register the computed
derivative in the registry as `d[field]_d[direction]`, where the values in the
brackets correspond to the choice of the user keywords.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "derivative",
   "name": "derivative",
   "field": "u",
   "direction": "y"
   "computed_field": "dudy"
 }
 ~~~~~~~~~~~~~~~

#### curl {#simcomp_curl}
Takes a list of three field names from the `fields` keyword, and computes
the curl.  By default, registers the result in `curl_x`, `curl_y` and `curl_z`.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "curl"
   "name": "curl"
   "fields": ["u", "v", "w"],
   "computed_field": "vorticity"
 }
 ~~~~~~~~~~~~~~~

#### divergence {#simcomp_divergence}
Takes a list of three field names from the `fields` keyword, and computes
the divergence.  By default, registers the result in `div`.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "divergence"
   "name": "divergence"
   "fields": ["u", "v", "w"],
   "computed_field": "continuity"
 }
 ~~~~~~~~~~~~~~~

### gradient {#simcomp_gradient}
Computes the gradient of a field.
The field to derivate is controlled by the `field` keyword. The simcomp will, by
default, register the computed components of the gradients in the registry as
`grad_[field]_x`, `grad_[field]_y`, `grad_[field]_z` where the
value in the brackets corresponds to the choice of the user keyword.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "gradient"
   "name": "gradient"
   "field": "u",
 }
 ~~~~~~~~~~~~~~~

### weak_gradient {#simcomp_weak_gradient}
Computes the weak gradient of a field. The weak gradient is value of the
gradient multiplied by the local value of the mass matrix. This is how a
gradient term appears in the weak formulation of the governing equations. The
field to derivate is controlled by the `field` keyword. The simcomp will, by
default, register the computed components of the gradients in the registry as
`weak_grad_[field]_x`, `weak_grad_[field]_y`, `weak_grad_[field]_z` where the
value in the brackets corresponds to the choice of the user keyword.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "weak_gradient"
   "name": "weak_gradient"
   "field": "u",
 }
 ~~~~~~~~~~~~~~~

### lambda2 {#simcomp_lambda2}
Computes \f$ \lambda_2 \f$ for the velocity field and stores it in the normal
output files as the first unused field. This means that \f$ \lambda_2 \f$ can be
found in the temperature field in then fld files if running without a scalar and
s1 if neko is run with one scalar. To output in a different `fld` series, use
the `"output_filename"` parameter.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "lambda2"
   "name": "lambda2"
 }
 ~~~~~~~~~~~~~~~

### probes {#simcomp_probes}
Probes selected solution fields at a list of points. This list of points can be
generated in a variety of ways, but the most common is to use the `csv` type.

Mandatory fields for this simcomp are:
- `fields`: a list of fields to probe. Should be a list of field names that
  exist in the registry. Example: `"fields": ["u", "v", "p", "s"]`.
- `output_file`: Name of the file in which to output the probed fields. Must be
  `.csv`.

Optional arguments:
- Interpolation parameters can be provided as a JSON sub-dictionary, 
  `interpolation`. If not provided, default values defined in 
  the `global_interpolation` module will be used.
  ~~~~~~~~~~~~~~~{.json}
  "interpolation": {
    "tolerance": 1e-8,
    "padding": 1e-3
  }
  ~~~~~~~~~~~~~~~

It is also possible to set a `start_time` before which the probes will not be
executed (same behavior as the statistics).

#### Supported types

- `file`: Reads a list of points from a CSV file. The name of the file is
   provided with the `file_name` keyword. The CSV file should have the
   following format:
   ~~~~~~~~~~~~~~~{.csv}
   x_0, y_0, z_0
   x_1, y_1, z_1
   ...
   x_N, y_N, z_N
   ~~~~~~~~~~~~~~~
   The points are assumed to be in the same units as the simulation.
- `points`: Reads a list of points from a JSON file. The points are specified
  based in the `coordinates` keyword and should be a list of x,y,z values.
  The file should have the following format:
  ~~~~~~~~~~~~~~~{.json}
  {
    "type": "points",
    "coordinates": [
      0.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0,
      ...
    ]
  }
  ~~~~~~~~~~~~~~~
  The points are assumed to be in the same units as the simulation.
- `line`: Generates a list of points along a line. The line is defined by two
  points, `start` and `end`, and the number of points to generate, `amount`.
  The points are generated by linearly interpolating between `start` and `end`.
  The line is defined as:
  ~~~~~~~~~~~~~~~{.json}
  {
    "type": "line",
    "start": [0.0, 0.0, 0.0],
    "end": [1.0, 1.0, 1.0],
    "amount": 10
  }
  ~~~~~~~~~~~~~~~
- `circle`: Generates a list of points along a circle. The circle is defined by
  a center, `center`, a radius, `radius` and the normal, `normal`, the number of
  points to generate is controlled by `amount`. The points are generated by
  rotating a point around the center starting from the specified axis projected
  onto the circle. The circle is defined:
  ~~~~~~~~~~~~~~~{.json}
  {
    "type": "circle",
    "center": [0.0, 0.0, 0.0],
    "radius": 1.0,
    "normal": [0.0, 0.0, 1.0],
    "axis": "x",
    "amount": 4
  }
  ~~~~~~~~~~~~~~~
  Leads to the following points:
  ~~~~~~~~~~~~~~~{.csv}
  1.0, 0.0, 0.0
  0.0, 1.0, 0.0
  -1.0, 0.0, 0.0
  0.0, -1.0, 0.0
  ~~~~~~~~~~~~~~~

 #### Example usage

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "probes",
   "name": "probes",
   "compute_control": "simulationtime",
   "compute_value"    : 1.0,
   "fields": ["w", "s"],
   "output_file":  "output.csv",
   "points": [
      {
        "type": "file",
        "file_name": "points.csv"
      }
    ],
 }
 ~~~~~~~~~~~~~~~
This probes the fields 'w', and 's' in the points described by points.csv and
outputs into output.csv every 1 time units.

The probed information will be saved in the output file in the following format:

~~~~~~~~~~~~~~~{.csv}
N_p, N_f, fields[0], fields[1], ..., fields[N_f-1]
p_0_x, p_0_y, p_0_z
p_1_x, p_1_y, p_1_z
...
p_N_p_x, p_N_p_y, p_N_p_z
time_0, p_0_field_0, p_0_field_1, ..., p_0_field_N_f-1
time_0, p_1_field_0, p_1_field_1, ..., p_1_field_N_f-1
...
time_0, p_N_p_field_0, p_N_p_field_1, ..., p_N_p_field_N_f-1
time_1, p_0_field_0, p_0_field_1, ..., p_0_field_N_f-1
time_1, p_1_field_0, p_1_field_1, ..., p_1_field_N_f-1
...
time_N_p, p_N_p_field_0, p_N_p_field_1, ..., p_N_p_field_N_f-1
~~~~~~~~~~~~~~~

### field_writer {#simcomp_field_writer}
Outputs registered 3D fields to a file. Requires a list of field names
in the `fields` keyword. Primarily to be used for outputting new fields defined
in the user file. The fields are added to the `neko_registry` object and
are expected to be updated in the user file, or, perhaps, by other simcomps.
Since this simcomp does not compute anything, `compute_` configuration is
irrelevant.

- The output format is controlled by the `output_format` keyword, which can be
  set to `nek5000` (default), `vtkhdf`, or `adios2`. 
- The `output_precision` keyword controls the precision of the written data 
  and can be set to `single` (default) or `double`.
- When using the `vtkhdf` format, the `output_subdivide` keyword can be set to
`true` to subdivide spectral elements into linear sub-cells instead of
writing high-order Lagrange cells. See the [cell representation](@ref
vtkhdf-cell-representation) section for more details.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "field_writer",
   "name": "field_writer",
   "fields": ["my_field1", "my_field2"],
   "output_filename": "myfields",
   "output_precision": "double",
   "output_format": "nek5000",
   "output_subdivide": false,
   "output_control" : "simulationtime",
   "output_value" : 1.0
 }
 ~~~~~~~~~~~~~~~


The `field_writer` may be used in conjunction with a `point_zone` to sample
the corresponding subsection of the domain. At the moment, this capability
can only be used with `nek5000` files.
@attention When using `point_zone` with the `nek5000` format, an 
`output_filename` must be provided.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "field_writer",
   "name": "field_writer",
   "fields": ["my_field1", "my_field2"],
   "output_format": "nek5000",
   "output_control" : "simulationtime",
   "output_value" : 1.0,
   "output_filename": "my_point_zone_field",
   "point_zone": "my_point_zone"
 }
 ~~~~~~~~~~~~~~~

### force_torque {#simcomp_force_torque}
Computes the force on a specified zone and the corresponding torque around a
center point. The compute control specifies how often they are computed and
printed into the log. Scale specifies a scale for the computed force/torque.
Conventient if one wants to scale with the area or similar. long_print is
default false and can be set to true to print all digits in the calculation.
Subroutines used in the simcomp can be found in src/qoi/drag_torque.f90

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "force_torque",
   "name": "force_torque",
   "zone_id": 1,
   "center_type": "fixed",
   "center": [0.0, 0.0, 0.0],
   "zone_name": "some chosen name, optional",
   "scale": 1.0
   "long_print" : false
   "compute_control" : "tsteps",
   "compute_value" : 10
 }
 ~~~~~~~~~~~~~~~

#### Torque calculation for moving bodies

When an object undergoes translational or rotational movement, it is often necessary to calculate the torque around its dynamic center of rotation, or around another specific reference point that moves rigidly with the body. The `center_type` parameter enables accurate torque computation for these scenarios by dictating how the tracking point behaves:

* `"fixed"` <i>(Default):</i> The torque is calculated around the static coordinates provided in the `center` array, regardless of how the body moves.
* `"pivot"` <i>(ALE only):</i> The torque is calculated directly around the ALE body's dynamic pivot point. If this is selected, the `center` array in the JSON is ignored, and the pivot coordinate at each time step is used automatically.
* `"body_attached"` <i>(ALE only):</i> The torque is calculated around a custom point that translates and rotates *with* the rigid movement of the ALE body. The initial position of this point is defined by the `center` array.
  > <i>Example use case:</i> If you are simulating a pitching and heaving airfoil, you might want to calculate the torque acting on a trailing-edge flap. By using `"body_attached"`, you simply define the initial coordinates of the hinge in the `center` array, and the code will automatically track its dynamic position as the main airfoil moves.

@attention For static simulations, the `center_type` parameter is completely optional. If omitted from the case file, the code will automatically default to `"fixed"`.

@note If `center_type` is set to `"pivot"` or `"body_attached"` but the specified `zone_id` is not registered as an ALE body (or ALE is globally inactive), the code will print a warning and automatically revert back to `"fixed"` using the provided `"center"` in the case file.

@attention For ALE simulations, the wall normal vectors are re-calculated at every time step to account for body movement and deformation. If the ALE module is not enabled, this calculation is performed only once during initialization.

@note **Restarting Simulations:** When restarting an ALE simulation, the code automatically calculates the correct current position of the torque center at the restart time. Therefore, if the intended torque calculation point remains the same, the `center` array in the JSON file should **not** be modified between restarts. If you wish to calculate torque around a *new* point upon restart, the `center` array must specify the coordinates of that new point in the **original, undeformed mesh** (at \f$ t=0 \f$), not its current spatial location.

### les_model {#simcomp_les_model}
Computes a subgrid eddy viscosity field using an SGS model. **Note*:* The simcomp
*only* computes the eddy viscosity field. You have to select the corresponding
`nut_field` in the fluid and/or scalar JSON object to actually enable LES, see
corresponding documentation. The simcomp is controlled by the following
keywords:

- `model`: Selects the SGS model. Currently available models are:
  - `smagorinsky`: The standard Smagorinsky model. Configured by the
    following additional keyword:
    - `c_s`: The Smagorinsky constant, defaults to 0.17.
  - `dynamic_smagorinsky`: The dynamic Smagorinsky model.
    - `test_filter`: The test filter for the dynamic Smagorinsky model
  - `vreman`: The Vreman model. Configured by the following additional keywords:
    - `c`: The model constant, defaults to 0.07.
    - `buoyancy_correction`: Whether or not to apply a correction to the eddy
      viscosity field based on the local Richardson number as described by Moeng
      and Sullivan 2015 (http://dx.doi.org/10.1016/B978-0-12-382225-3.00201-2).
      Defaults to `false`.
      - `true`: Add a buoyancy correction according to the following parameters:
        - `scalar_field`: Name of the scalar field based on which the buoyancy
          effect is computed.
        - `Ri_c`: The critical Richardson number.
        - `reference_temperature`: The reference temperature for computation of
          the Richardson number.
        - `g`: The gravity vector.
      - `false`: Compute the standard Vreman eddy viscosity.
  - `sigma`: The Sigma model. Configured by the following additional keyword:
    - `c`: The model constant, defaults to 1.35.
  - `wale`: The WALE model. Configured by the following additional keyword:
    - `c_w`: The WALE constant, defaults to 0.55.
  - `deardorff`: The Deardorff model dedicated for atmospheric boundary-layer
    applications. Please find the usage in `examples/shear_convection_ABL`.
    Configured by the following additional keyword:
    - `c_k`: The model constant, defaults to 0.1.
    - `T0`: The reference temperature.
    - `g`: The gravity vector.
    - `temperature_field`: The field name of the temperature field,
      defaults to `temperature`.
    - `TKE_field`: The field name of the turbulent kinetic energy (TKE) field,
      defaults to `TKE`.
    - `temperature_alphat_field`: The field name of the eddy diffusivity field
      for the temperature equation, defaults to `temperature_alphat`.
    - `TKE_alphat_field`: The field name of the eddy diffusivity field
      for the TKE equation, defaults to `TKE_alphat`.
    - `TKE_source_field`: The field name of the source terms in the TKE equation
      including shear production, buoyancy contribution and dissipation,
      defaults to `TKE_source`.
- `les_delta`: Selects the way to compute the LES filter length scale. Currently
  three alternatives are provided and the default one is `pointwise` if nothing
  is specified:
  - `pointwise`: Computes a local value based on the spacing of the GLL nodes.
  - `elementwise_average`: Computes a single value for the whole element based
    on the average spacing of the GLL nodes within the element.
  - `elementwise_max`: Computes a single value for the whole element based on
    the maximum spacing of the GLL nodes within the element. The `les_delta`
  field is added to the registry and written to the .fld files.
- `nut_field`: The name of the SGS eddy viscosity field added to the registry.
  Defaults to `nut`. This allows to have two different SGS models active, saved
  to different fields. For example, one for the scalar and one to the fluid.
- `extrapolation`: Whether or not extrapolate the velocity to
  compute the eddy viscosity.
  - `true`: extrapolate the velocity as the same order as
  the time scheme.
  - `false`: the default option, disable the extrapolation.
  In this case, the estimation of the eddy viscosity is of first order, while
  circumvent the risk of unstable extrapolation.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "les_model"
   "name": "les_model"
   "model": "smagorinsky",
   "delta_type": "pointwise",
   "output_control" : "never"
 }
 ~~~~~~~~~~~~~~~

 Please also note that for the dynamic Smagorinsky model, one needs to specify
 the test filter in the following way for the Boyd filter as the test filter
 (one could also use "nonBoyd" as the option):
 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "les_model"
   "model": "dynamic_smagorinsky",
   "test_filter": {
      "filter": {
        "type": "elementwise",
        "elementwise_filter_type": "Boyd"
      }
    }
 }
 ~~~~~~~~~~~~~~~
 And one could not change the default test filter's kernel through the case
 file. If one needs to do so, he/she needs to dig into the code in
 src/les/dynamic_smagorinksy.f90.

 ### User statistics {#user_stats}

 Allows to compute the time-average of an arbitrary collection of fields form
 the field registry. Just like the `fluid_stats` simcomp, it supports spatial
 averaging across homogeneous directions, both 1D and 2D.  The fields to average
 are prescribed via the `fields` keyword, and the optional averaging
 direction(s) via the `avg_direction`, which can be `x`, `y`, `z`, `xy`, `xz` or
  `yz`. Averaging across two directions will lead to the average being saved as
  a .csv, whereas a 2D .fld file will be produced when averaging across only one
  axis. The filename is controlled  by the `output_file` keyword and default to
  `user_stats`. We encourage reading the [statistics guide](@ref
  statistics-guide) for further details regarding how statistics are computed in
  Neko.

  Keep in mind that simcomps execute before `user%compute`, so if you update
  some custom averaged field in that routine, it will not affect the average
  until the next time step. You can consider using `user%preprocess` instead,
  which runs at the beginning of the time step.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "user_stats",
   "name": "user_stats",
   "fields": ["s"],
   "avg_direction": "xz",
   "output_file": "s_average"
 }
 ~~~~~~~~~~~~~~~

The statistics fields created by this simcomp are accessible from the
neko registry and retrievable under the following naming convention:
`name_in_registry = name_of_simcomp + "/mean_" + name_of_field`. Unless
specified, the name of the simcomp will default to `user_stats`.
For example, if `"fields": ["s", "my_field"]` and `"name": "my_stats"` then
the fields `"my_stats/mean_s"` and `"my_stats/mean_my_field"` will be added
to the registry.

### Spectral error indicator {#simcomp_speri}

Computes the spectral error indicator as developed by Mavriplis (1989)
(https://doi.org/10.1007/978-3-663-13975-1_34). This is an a posteriori error
measure, based on the local properties of the spectral solution. This method
formally only gives an indication of the error.

The spectral error indicator is computed for the 3 velocity fields, resulting
in 3 additional fields appended to the field files.

~~~~~~~~~~~~~~~{.json}
 {
   "type": "spectral_error"
   "name": "spectral_error"
 }
 ~~~~~~~~~~~~~~~

### Data streamer {#simcomp_data_streamer}

Enables data streaming of a set of given `fields` with the `ADIOS2` library. 
The simcomp is controlled by the following keywords:
- `"fields"`: A list of field names corresponding to the fields to stream 
  (must exist in the registry). The fields will be streamed in the order
  given in the list.
- `"stream_mesh"`: Whether or not to stream mesh coordinates, in the order
  `x`, `y`, `z`. The mesh coordinates will always be streamed first, in
  that exact order, before the fields in `"fields"`.

See the `cylinder` or `turb_pipe` examples for more details on how this 
simcomp cam be coupled to Python scripts for in-situ data processing.

@note This simcomp requires configuration of Neko with the ADIOS2 library
(`--with-adios2=DIR`).

~~~~~~~~~~~~~~~{.json}
 {
   "type": "data_streamer",
   "name": "spectral_error",
   "fields": ["u", "omega_z", "fluid_stats/mean_u"],
   "stream_mesh": true,
   "compute_control": "tsteps",
   "compute_value": 10
 }
 ~~~~~~~~~~~~~~~

### Field subsampler {#simcomp_field_subsampler}

Creates sub-sections of the domain from a `point_zone` and/or at a lower
`polynomial_order`. The fields are added to the registry under the name
`name_of_simcomp + "/" + name_of_base_field`. For example, 
`field_subsampler_u`.

The simcomp is controlled by the following keywords:
- `"source_fields"`: A list of names corresponding to the fields to subsample 
  (must exist in the registry).
- `point_zone` (optional): The name of the point zone to use to mask the fields.
- `polynomial_order` (optional): The new polynomial at which to interpolate 
  the fields. Must be different from the order used in the simulation.

The `field_subsampler` contains its own `field_writer`. Therefore, all the
keywords used by the latter can also be specified, with the exception of
`point_zone` and `fields` which will be handled by the `field_subsampler`.

~~~~~~~~~~~~~~~{.json}
 {
   "type": "field_subsampler",
   "name": "field_subsampler",
   "source_fields": ["u", "omega_z", "fluid_stats/mean_u"],
   "point_zone": "my_point_zone",
   "polynomial_order": 3,
   "compute_control": "tsteps",
   "compute_value": 10
 }
 ~~~~~~~~~~~~~~~
