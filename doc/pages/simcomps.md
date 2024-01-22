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
Each simulation component is, by default, executed once per time step to
perform associated computations and output.
However, this can be modified by using the `compute_control` and `compute_value`
parameters for the computation and the `output_control and` and 
`output_value` for the output to disk.
The parameters for the `_control` values are the same as for the fluid and 
checkpointing.

For example, in the `tgv` example case the `vorticity` component is executed 
once per 50 time steps. 
~~~~~~~~~~~~~~~{.json}
{
    "type": "vorticity",
    "compute_control": "tsteps",
    "compute_value": 50
}
~~~~~~~~~~~~~~~
If no parameters for the `output_` parameters are provided, they are set to be the
 same as for `compute_`.

 ## List of simulation components

 ### Vorticity
 Computes the vorticity field an stores in the field registry as `omega_x`,
 `omega_y` and `omega_z`.
 Currently produces no output.

 ### lambda2
 Computes \f$ \lambda_2 \f$ for the velocity field and stores it in the normal output files as the first unused field.
 This means that \f$ \lambda_2 \f$ can be found in the temeprature field in then fld files if running without a scalar
 and s1 if neko is run with one scalar.
 
 ### Probes
 Probes selected solution fields at given `x,y,z` coordinates. The input coordinates can be given in an input file or can be generated along a line spanning from `start_point` to `end_point`. The table below su

Name                | Description                                    |    Admissable values     | Default value
--------------------|------------------------------------------------|--------------------------|--------------
`input_type`        | Select how the probe locations should be read. | "file" or "line"         | "file"
`input_file`        | Path to the file containing the probe coordiates | String with file extension `.csv`         | -
`start_point`       | Coordinates of the starting point of the line. | Vector of 3 reals        | -
`end_point`         | Coordinates of the ending point of the line.   | Vector of 3 reals        | -
`n_samples`         | Number of points to create along the line.     | Integer                  | -
`compute_control`   | Defines the interpretation of `compute_value` to define the frequency of computing (and outputting) probes. | `nsamples`, `simulationtime`, `tsteps`, `never` | - 
`compute_value`     | The frequency of computing and outputting in terms of `compute_control`. | Positive real or integer | -

**Probes locations from an input file**

This example probes the fields 'w', and 's' in the points described by points.csv and outputs into output.csv every 1 time units.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "probes",
   "compute_control": "simulationtime",
   "compute_value"    : 1,
   "input_type": "file",
   "points_file":  "probes.csv",
   "output_file":  "output.csv",
   "fields": ["w","s"]
 }
 ~~~~~~~~~~~~~~~

With `probes.csv` containing, for example:

 ~~~~~~~~~~~~~~~{.csv}
 0.1,0.0,0.0
 0.1,0.5,0.0
 0.1,0.0,2.0
 0.1,0.5,2.0
 ~~~~~~~~~~~~~~~

**Probes along a line**


Example usage:
 ~~~~~~~~~~~~~~~{.json}
 {
   "type"            : "probes",
   "input_type"      : "line",
   "start_point"     : [0.0, 0.0, 0.0],
   "end_point"       : [1.0, 0.0, 0.0],
   "n_samples"       : 100,
   "output_file"     : "output.csv",
   "fields"          : ["w","s"],
   "compute_control" : "simulationtime",
   "compute_value"   : 1
 }
 ~~~~~~~~~~~~~~~
This probes the fields 'w', and 's' for 100 points on a line along the x-axis 
and outputs into output.csv every 1 time units.
