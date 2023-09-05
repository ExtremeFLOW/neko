
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
The paramters for the `_control` values are the same as for the fluid and 
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
If no parameters for the `output_` parametersare provided, they areset to be the
 same as for `compute_`.

 ## List of simulation components

 ## Vorticity
 Computes the vorticity field an stores in the field registry as `omega_x`,
 `omega_y` and `omega_z`.
 Currently produces no output.





