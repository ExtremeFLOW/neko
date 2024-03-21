# Simulation components {#simcomps}
## What are simulation components?
Simulation components, or simcomps fo short,  incapsulate additional
functionality that may be useful for certain cases but not necessary to run the
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

## Controling execution and file output
Each simulation component is, by default, executed once per time step to perform
associated computations and output. However, this can be modified by using the
`compute_control` and `compute_value` parameters for the computation and the
`output_control and` and `output_value` for the output to disk. The parameters
for the `_control` values are the same as for the fluid and checkpointing.
Additionally, one can set `output_control` to `global` and `never`. The former
will sync the `output_` parameter to that of the fluid. Choosing `never` will
suppress output all together. If no parameters for the `output_` parameters are
 provided, they are set to be the same as for `compute_`.

For simcomps that compute 3D fields, the output can be either added to the main
`.fld` file, containing velocity and pressure, or saved to a separate file. For
the latter, the `output_filename` keyword should be provided. One can
additionally provide the `precision` keyword, which can be set to either
`single` or `double` to control the precision of the written data.

For example, in the `tgv` example case the `vorticity` component is executed
once per 50 time steps. The `output_` parameters are synced to that, and the
vorticity fields will be added to the main `.fld` file.
~~~~~~~~~~~~~~~{.json}
{
    "type": "vorticity",
    "compute_control": "tsteps",
    "compute_value": 50
}
~~~~~~~~~~~~~~~

 ## List of simulation components

 ### vorticity
 Computes the vorticity field an stores in the field registry as `omega_x`,
 `omega_y` and `omega_z`.
 Currently produces no output.

 ### lambda2
 Computes \f$ \lambda_2 \f$ for the velocity field and stores it in the normal output files as the first unused field.
 This means that \f$ \lambda_2 \f$ can be found in the temeprature field in then fld files if running without a scalar
 and s1 if neko is run with one scalar.

 ### probes
 Probes selected solution fields at the points given inside an input file. Example usage:
 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "probes",
   "compute_control": "simulationtime",
   "compute_value"    : 1,
   "points_file":  "probes.csv",
   "output_file":  "output.csv",
   "fields": ["w","s"]
 }
 ~~~~~~~~~~~~~~~
This probes the fields 'w', and 's' in the points described by points.csv and outputs into output.csv every 1 time units.

 ### field_writer
 Outputs registered 3D fields to an `.fld` file. Requires a list of field names
 in the `fields` keyword. Primarily to be used for outputting new fields defined
 in the user file. The fields are added to then `neko_field_registry` object and
 are expected to be updated in the user file, or, perhaps, by other simcomps.
 Since this simcomp does not compute anything `compute_` configuration is
 irrelevant.
 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "field_writer",
   "fields": ["my_field1", "my_field2"],
   "output_filename": "myfields",
   "precision": "double",
   "output_control" : "simulation_time",
   "output_value" : 1.0
 }
 ~~~~~~~~~~~~~~~
