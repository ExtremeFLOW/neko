# Statistics guide {#statistics-guide}

\tableofcontents

Statistics in the context of Neko, is the common name for fields that are averaged in time and possible also space.

The statistics module in Neko computes the temporal average of a wide range of fields.

In this page we use the following convention for a field
- \f$ u \f$, the instantaneous field.
- \f$ \langle u \rangle_t \f$, the temporal average of \f$ u \f$. 
- \f$ u = \langle u \rangle + u' \f$, the Reynolds decomposition of \f$ u \f$, where \f$ u' \f$ is the fluctuation of \f$ u \f$ around the mean field.

The temporal average of a field \f$u\f$ is the approximation of the integral

$$
\langle u \rangle_t = \int_{T_0}^{T_N} u dt
$$

In Neko, this is computed as

$$
\langle u \rangle_t = \sum_{i=0}^N u_i \Delta t_i
$$
where \f$ u_0 \f$ is the fields value at \f$ T_0 \f$ and \f$ N \f$ is the number of time steps needed to reach \f$ T_N \f$, \f$ T_N = T_0 + \sum_{i=0}^N \Delta t_i \f$.

# Fluid Statistics

In the fluid statistics in Neko, various averages of the different velocity components, derivatives and pressure are computed. In total, 44 "raw statistics" are computed that are required to compute the Reynolds stress budgets, mean fields, and the different terms in the turbulent kinetic energy equation.

## Using statistics
Statistics are enabled in the case file as a simcomp with the added argument `avg_direction`, `set_of_stats`, and `start_time`:

| Name                | Description                                                          | Admissible values | Default value |
| ------------------- | -------------------------------------------------------------------- | ----------------- | ------------- |
| `start_time`        | Time at which to start gathering statistics.                        | Positive real     | 0             |
| `avg_direction`        | Directions to compute spatial average.                         | x,y,z,xy,xz,yz  |  No spatial average           |
| `set_of_stats`        | What set of stats to compute.                         | basic, full  |  full         |
| `compute_value` | Interval, in timesteps or simulationtime, depending on compute\_control, for sampling the flow fields for statistics. | Positive real or int  | Not set (but recommended with every 50 timesteps or so  |
| `output_filename`        | User-specified filename to store output in.                       | filename  |  fluid_statsX*        |

\*The name of the written statistics file will by default be `fluid_statsX0.f0000X,..., fluid_statsX0.f0000Y` where X is the number of the first outputted statistic of the current run.

In addition, one can specify the usual controls for the output, which then outputs the averages computes from the last time the statistics were written to file. For example, if one wants to compute only the basic statistics and sample the fields every 4 time steps and compute and output batches every 20 time units and have an initial transient of 60 time units the following would work:

~~~~~~~~~~~~~~~{.json}
"simulation_components": 
  [
    {
      "type": "fluid_stats",
      "compute_control": "tsteps",
      "compute_value": 4,
      "output_control": "simulationtime",
      "output_value": 20,
      "start_time":60.0,
      "avg_direction":"xz",
      "set_of_stats":"basic"
    }
  ]
~~~~~~~~~~~~~~~

@attention For simulations requiring restarts, it is recommended to run each 
restart in a different output directory as a precaution to avoid potential overwritings of files.

Preferably set the initial transient to a multiple of output_value as otherwise the first output will be slightly shorter than the rest. The code related to fluid statistics are located in fluid_stats and fluid_stats_simcomp. 

The argument "avg_direction" is optional and if ignored we output 3d fields. The statistics are saved in a fld file according to the following in 2D and 3D. Observe that in 2D the mean Z-velocity is stored in a last scalar field. All other fields are kept the same. This is due to a limitation of the fld file format.

For 1D statistics a CSV file is outputted. The first column is the time at which the statistics are collected, the second column the spatial coordinate, and the rest of the data is stored in the order below. In this case all statistics are kept in the same order as in 3D. The name for these files are `fluid_statsX.csv,..., fluid_statsX.csv` where X is the number of the first outputted statistic of the current run.


## List of fields in output files

When only the basic set of stats is enabled, only stats 1-11 are computed. When 2D stats are enabled \f$ \langle w \rangle \f$ is stored in s7 for basic stats and in s40 for the full set of statistics.

| Number | Statistic | Stored in variable (for fld files) |
| ------ | --------- | ------------------ |
| 1 | \f$ \langle p \rangle \f$ | Pressure|
| 2 | \f$ \langle u \rangle \f$ | X-Velocity|
| 3 | \f$ \langle v \rangle \f$ | Y-Velocity|
| 4 | \f$ \langle w \rangle \f$ | Z-Velocity |
| 5 | \f$ \langle pp \rangle \f$ | Scalar 1 (s1)|
| 6 | \f$ \langle uu \rangle \f$ | Scalar 2 (s2)|
| 7 | \f$ \langle vv \rangle \f$ | s3|
| 8 | \f$ \langle ww \rangle \f$ | s4|
| 9 | \f$ \langle uv \rangle \f$ | s5|
|10 | \f$ \langle uw \rangle \f$ | s6|
| 11 | \f$ \langle vw \rangle \f$ | s7|
| 12| \f$ \langle uuu \rangle \f$ | s8|
| 13| \f$ \langle vvv \rangle \f$ | s9|
| 14 | \f$ \langle www \rangle \f$ | s10|
| 15 | \f$ \langle  uuv   \rangle \f$ | s11 |
| 16 | \f$ \langle  uuw   \rangle \f$ | s12 |
| 17 | \f$ \langle  uvv   \rangle \f$ | s13 |
| 18 | \f$ \langle  uvw   \rangle \f$ | s14 |
| 19 | \f$ \langle  vvw   \rangle \f$ | s15 |
| 20 | \f$ \langle  uww   \rangle \f$ | s16 |
| 21 | \f$ \langle  vww   \rangle \f$ | s17 |
| 22 | \f$ \langle  uuuu  \rangle \f$ | s18 |
| 23 | \f$ \langle  vvvv  \rangle \f$ | s19 |
| 24 | \f$ \langle wwww   \rangle \f$ | s20 |
| 25 | \f$ \langle  ppp   \rangle \f$ | s21 |
| 26 | \f$ \langle  pppp  \rangle \f$ | s22 |
| 27 | \f$ \langle  pu    \rangle \f$ | s23 |
| 28 | \f$ \langle  pv    \rangle \f$ | s24 |
| 29 | \f$ \langle  pw    \rangle \f$ | s25 |
| 30 | \f$ \langle  p \frac{\partial u} {\partial x} \rangle \f$ | s26 |
| 31 | \f$ \langle  p \frac{\partial u} {\partial y}\rangle \f$ | s27 |
| 32 | \f$ \langle  p \frac{\partial u} {\partial z}\rangle \f$ | s28 |
| 33 | \f$ \langle  p \frac{\partial v} {\partial x}\rangle \f$ | s29 |
| 34 | \f$ \langle  p \frac{\partial v} {\partial y}\rangle \f$ | s30 |
| 35 | \f$ \langle  p \frac{\partial v} {\partial z}\rangle \f$ | s31 |
| 36 | \f$ \langle  p \frac{\partial w} {\partial x}\rangle \f$ | s32 |
| 37 | \f$ \langle  p \frac{\partial w} {\partial y}\rangle \f$ | s33 |
| 38 | \f$ \langle  p \frac{\partial w} {\partial z}\rangle \f$ | s34 |
| 39 | \f$ \langle  e11   \rangle \f$ | s35 |
| 40 | \f$ \langle  e22   \rangle \f$ | s36 |
| 41 | \f$ \langle  e33   \rangle \f$ | s37 |
| 42 | \f$ \langle  e12   \rangle \f$ | s38 |
| 43 | \f$ \langle  e13   \rangle \f$ | s39 |
| 44 | \f$ \langle  e23   \rangle \f$ | s40 |

where \f$e11,e22...\f$ is computed as:
$$
\begin{aligned}
e11 &= \left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial u}{\partial y}\right)^2 + \left(\frac{\partial u}{\partial z}\right)^2 \\\\
e22 &= \left(\frac{\partial v}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2 + \left(\frac{\partial v}{\partial z}\right)^2 \\\\
e33 &= \left(\frac{\partial w}{\partial x}\right)^2 + \left(\frac{\partial w}{\partial y}\right)^2 + \left(\frac{\partial w}{\partial z}\right)^2 \\\\
e12 &= \frac{\partial u}{\partial x}  \frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}\frac{\partial v}{\partial y}+ \frac{\partial u}{\partial z}\frac{\partial v}{\partial z} \\\\
e13 &= \frac{\partial u}{\partial x}  \frac{\partial w}{\partial x} + \frac{\partial u}{\partial y}\frac{\partial w}{\partial y}+ \frac{\partial u}{\partial z}\frac{\partial w}{\partial z} \\\\
e23 &= \frac{\partial v}{\partial x}  \frac{\partial w}{\partial x} + \frac{\partial v}{\partial y}\frac{\partial w}{\partial y}+ \frac{\partial v}{\partial z}\frac{\partial w}{\partial z} \\\\
\end{aligned}
$$


# Postprocessing
These statistics are only the "raw statistics" in the sense that in general we are not interested in \f$ \langle uu\rangle \f$, but rather say \f$ \langle u'u'\rangle\f$. For this we need to postprocess the statistics. 

There is some rudimentary postprocessing to compute the spatial averages of fld files and also to combine the statistics collected from several runs (compute average in time) and also compute both the mean velocity gradient and the Reynolds stresses available in the contrib scripts. By running the contrib scripts without any arguments one gets a hint on their usage, (e.g. by running `./average_fields_in_time` will give the options). 

To further postprocess the statistics it is suggested to look into PyNekTools which introduces convenient functions for postprocessing, largely based on the PyMech library with the addition of an expanding set of tools for interpolation, computing derivatives and more advanced functionality.  PyNekTools is entirely parallelized in MPI and can also handle large data sets for postprocessing.



# Scalar Statistics

In the scalar statistics in Neko, various averages of a given scalar and the different velocity components, derivatives and pressure are computed, in a similar fashion to the fluid statistics. In total, 42 "raw statistics" are computed that are required to compute the mean scalar transport equation, skewness and kurtosis, as well as the scalar variance budget, and turbulent scalar flux budgets.

## Using statistics
Similar to fluid statistics, scalar statistics are enabled in the case file as a simcomp with an additional argument `field` for the name of the scalar field to be averaged:

| Name                | Description                                                          | Admissible values | Default value |
| ------------------- | -------------------------------------------------------------------- | ----------------- | ------------- |
| `field`        | Name of the scalar field to be averaged                       | String  |  `s`        |
| `start_time`        | Time at which to start gathering statistics.                        | Positive real     | 0             |
| `avg_direction`        | Directions to compute spatial average.                         | x,y,z,xy,xz,yz  |  No spatial average           |
| `set_of_stats`        | What set of stats to compute.                         | basic, full  |  full         |
| `compute_value` | Interval, in timesteps or simulationtime, depending on compute\_control, for sampling the flow fields for statistics. | Positive real or int  | Not set (but recommended with every 50 timesteps or so)  |
| `output_filename`        | User-specified filename to store output in.                       | filename  |  fluid_statsX*        |



\*The name of the written statistics file will by default be `scalar_statsX0.f0000X,..., scalar_statsX0.f0000Y` where X is the number of the first outputted statistic of the current run. Note that if you want to compute statistics for multiple scalars, you will need to specify an independent `output_filename` for each.

In addition, one can specify the usual controls for the output, in the same manner as for fluid statistics. For example, if one wants to compute only the basic statistics and sample the fields every 4 time steps and compute and output batches every 20 time units and have an initial transient of 60 time units the following would work:

~~~~~~~~~~~~~~~{.json}
"simulation_components": 
  [
    {
      "type": "scalar_stats",
      "field": "s",
      "compute_control": "tsteps",
      "compute_value": 4,
      "output_control": "simulationtime",
      "output_value": 20,
      "start_time": 60.0,
      "avg_direction":"xz",
      "set_of_stats":"basic"
    }
  ]
~~~~~~~~~~~~~~~



## List of fields in output files

When only the basic set of stats is enabled, only stats 1-5 are computed.

| Number | Statistic | Stored in variable (for fld files) |
| ------ | --------- | ------------------ |
| 1 | \f$ \langle s \rangle \f$ | Pressure|
| 2 | \f$ \langle us \rangle \f$ | X-Velocity|
| 3 | \f$ \langle vs \rangle \f$ | Y-Velocity|
| 4 | \f$ \langle ws \rangle \f$ | Z-Velocity |
| 5 | \f$ \langle ss \rangle \f$ | Scalar 1 (s1)|
| 6 | \f$ \langle sss \rangle \f$ | Scalar 2 (s2)|
| 7 | \f$ \langle ssss \rangle \f$ | s3|
| 8 | \f$ \langle uss \rangle \f$ | s4|
| 9 | \f$ \langle vss \rangle \f$ | s5|
|10 | \f$ \langle wss \rangle \f$ | s6|
| 11 | \f$ \langle uus \rangle \f$ | s7|
| 12| \f$ \langle vvs \rangle \f$ | s8|
| 13| \f$ \langle wws \rangle \f$ | s9|
| 14 | \f$ \langle uvs \rangle \f$ | s10|
| 15 | \f$ \langle  uws   \rangle \f$ | s11 |
| 16 | \f$ \langle  vws   \rangle \f$ | s12 |
| 17 | \f$ \langle  ps   \rangle \f$ | s13 |
| 18 | \f$ \langle  p \frac{\partial s}{\partial x}  \rangle \f$ | s14 |
| 19 | \f$ \langle  p \frac{\partial s}{\partial y}  \rangle \f$ | s15 |
| 20 | \f$ \langle  p \frac{\partial s}{\partial z}  \rangle \f$ | s16 |
| 21 | \f$ \langle  u \frac{\partial s}{\partial x}  \rangle \f$ | s17 |
| 22 | \f$ \langle  u \frac{\partial s}{\partial y}  \rangle \f$ | s18 |
| 23 | \f$ \langle  u \frac{\partial s}{\partial z}  \rangle \f$ | s19 |
| 24 | \f$ \langle  v \frac{\partial s}{\partial x}  \rangle \f$ | s20 |
| 25 | \f$ \langle  v \frac{\partial s}{\partial y}  \rangle \f$ | s21 |
| 26 | \f$ \langle  v \frac{\partial s}{\partial z}  \rangle \f$ | s22 |
| 27 | \f$ \langle  w \frac{\partial s}{\partial x}  \rangle \f$ | s23 |
| 28 | \f$ \langle  w \frac{\partial s}{\partial y}  \rangle \f$ | s24 |
| 29 | \f$ \langle  w \frac{\partial s}{\partial z}  \rangle \f$ | s25 |
| 30 | \f$ \langle  s \frac{\partial u}{\partial x} \rangle \f$ | s26 |
| 31 | \f$ \langle  s \frac{\partial u}{\partial y} \rangle \f$ | s27 |
| 32 | \f$ \langle  s \frac{\partial u}{\partial z} \rangle \f$ | s28 |
| 33 | \f$ \langle  s \frac{\partial v}{\partial x} \rangle \f$ | s29 |
| 34 | \f$ \langle  s \frac{\partial v}{\partial x} \rangle \f$ | s30 |
| 35 | \f$ \langle  s \frac{\partial v}{\partial x} \rangle \f$ | s31 |
| 36 | \f$ \langle  s \frac{\partial w}{\partial x} \rangle \f$ | s32 |
| 37 | \f$ \langle  s \frac{\partial w}{\partial x} \rangle \f$ | s33 |
| 38 | \f$ \langle  s \frac{\partial w}{\partial x} \rangle \f$ | s34 |
| 39 | \f$ \langle  e_{ss}   \rangle \f$ | s35 |
| 40 | \f$ \langle  e_{us}   \rangle \f$ | s36 |
| 41 | \f$ \langle  e_{vs}   \rangle \f$ | s37 |
| 42 | \f$ \langle  e_{ws}   \rangle \f$ | s38 |

where \f$e_{ss}, e_{us}, e_{vs}, e_{ws}\f$ are computed as:
$$
\begin{aligned}
e_{ss} &= \left(\frac{\partial s}{\partial x}\right)^2 + \left(\frac{\partial s}{\partial y}\right)^2 + \left(\frac{\partial s}{\partial z}\right)^2 \\\\
e_{us} &= \frac{\partial u}{\partial x}  \frac{\partial s}{\partial x} + \frac{\partial u}{\partial y}\frac{\partial s}{\partial y}+ \frac{\partial u}{\partial z}\frac{\partial s}{\partial z} \\\\
e_{vs} &= \frac{\partial v}{\partial x}  \frac{\partial s}{\partial x} + \frac{\partial v}{\partial y}\frac{\partial s}{\partial y}+ \frac{\partial v}{\partial z}\frac{\partial s}{\partial z} \\\\
e_{ws} &= \frac{\partial w}{\partial x}  \frac{\partial s}{\partial x} + \frac{\partial w}{\partial y}\frac{\partial s}{\partial y}+ \frac{\partial w}{\partial z}\frac{\partial s}{\partial z} \\\\
\end{aligned}
$$