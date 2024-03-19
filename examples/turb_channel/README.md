#Simple turbulent channel case at Re_tau=180
The channel starts with a slightly perturbed initial condition which becomes turbulent after around 10 time units.
If one wants to change the mesh and play around, please do so, by generating a new mesh with genmeshbox.
To compute thte spatial and temporal averages of the output fields, mean_fields and stats one can use the following contrib scripts:
- average_fields_in_time, averages a fld series in time
- average_field_in_space, for channelflow we can directly average in xz
- postprocess_fluid_stats, computes the mean gradients and the reynolds stresses. 
