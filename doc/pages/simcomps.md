
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
 
Name              | Description                                                                         |    Admissable values               | Default value
------------------|---------------------------------------------------------|------------------------------------|--------------
`type`            | Selects the type of simulation component.               | `vorticity`, `lambda2` or `probes` | -
`compute_control` | Defines the interpretation of `compute_value` to define the frequency of computation of a simulation component. | `nsamples`, `simulationtime`, `tsteps` or `never` | `tsteps`
`compute_value`   | The frequency of computation in terms of `compute_control`. | Positive real or integer     | 1
`output_control` | Defines the interpretation of `output_value` to define the output frequency of a simulation component. | `nsamples`, `simulationtime`, `tsteps` or `never` | `compute_control`
`output_value`   | The output frequency in terms of `output_control`. | Positive real or integer     | `compute_value`

## List of simulation components

### Vorticity
Computes the vorticity field an stores in the field registry as `omega_x`,
`omega_y` and `omega_z`. Currently produces no output.

### lambda2
Computes \f$ \lambda_2 \f$ for the velocity field and stores it in the normal
output files as the first unused field. This means that \f$ \lambda_2 \f$ 
can be found in the temperature field in then fld files if running without
a scalar and s1 if neko is run with one scalar.
 
### Probes
Used for probing selected solution fields at given points. The coordinates of the probes 
can be specified in two fashions, from a `.csv` file or along a line with
specified `start` and `end` points.

**General parameters**

Name                | Description                                             |    Admissable values     | Default value
--------------------|---------------------------------------------------------|--------------------------|--------------
`input_type`        | Selects the type of input for the probe coordinates.    | `file` or `line`         | `file`
`output_file`       | File in which to output the probed fields.              | String ending with `.csv`| -
`fields`            | Labels of the fields to probe.                          | Vector of characters     | -

A different set of parameters will be looked up in the json object 
depending on the value set for `input_type`. They are as follows:

**Probes from a `.csv` file**

Name                | Description                                                                |    Admissable values     | Default value
--------------------|----------------------------------------------------------------------------|--------------------------|--------------
`points_file`       | File from which to read the probe coordinates.                             | String ending with `.csv`| -

Example usage:
 ~~~~~~~~~~~~~~~{.json}
 {
   "type"           : "probes",
   "points_file"    : "probes.csv",
   "output_file"    : "output.csv",
   "fields"         : ["w","s"],
   "compute_control": "simulationtime",
   "compute_value"  : 1,
 }
 ~~~~~~~~~~~~~~~
This probes the fields 'w', and 's' in the points described by points.csv and outputs into output.csv every 1 time units.

**Probes along a line**

Name                | Description                                    |    Admissable values     | Default value
--------------------|------------------------------------------------|--------------------------|--------------
`start_point`       | Coordinates of the starting point of the line. | Vector of 3 reals        | -
`end_point`         | Coordinates of the ending point of the line.   | Vector of 3 reals        | -
`n_samples`         | Number of points to create along the line.     | Integer                  | -

Example usage:
 ~~~~~~~~~~~~~~~{.json}
 {
   "type"            : "probes",
   "input_type"      : "line"
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
