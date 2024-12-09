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

In the statistics in Neko, various averages of the the different velocity components, derivatives and pressure are computed. In total, 44 "raw statistics" are computed that are required to compute the Reynolds stress budgets, mean fields, and the different terms in the turbulent kinetic energy equation.

## Using statistics
Statistics are enable in the the case file as a simcomp with the added argument `avg_direction`, `set_of_stats`, and `start_time`:

| Name                | Description                                                          | Admissible values | Default value |
| ------------------- | -------------------------------------------------------------------- | ----------------- | ------------- |
| `start_time`        | Time at which to start gathering statistics.                        | Positive real     | 0             |
| `avg_direction`        | Directions to compute spatial average.                         | x,y,z,xy,xz,yz  |  No spatial average           |
| `set_of_stats`        | What set of stats to compute.                         | basic, full  |  full         |
| `compute_value` | Interval, in timesteps or simulationtime, depending on compute\_control, for sampling the flow fields for statistics. | Positive real or int  | Not set (but recommended with every 50 timesteps or so  |

In addition to the usual controls for the output, which then outputs the averages computes from the last time the statistics were written to file.

For example, if one wants to compute only the basic statistics and sample the fields every 4 time steps and compute and output batches every 20 time units and have an initial transient of 60 time units the following would work:

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
    }1
  ]
~~~~~~~~~~~~~~~

Preferably set the initial transient to a multiple of output_value as otherwise the first output will be slightly shorter than the rest. The code related to fluid statistics are located in fluid_stats and fluid_stats_simcomp.

The argument "avg_direction" is optional and if ignored we output 3d fields. The statistics are saved in a fld file according to the following in 2D and 3D. Observe that in 2D the mean Z-velocity is stored in a last scalar field. All other fields are kept the same. This is due to a limitation of the fld file format.

For 1D statistics a CSV file is outputted. The first column is the time at which the statistics are collected, the second column the spatial coordinate, and the rest of the data is stored in the order below. In this case all statistics are kept in the same order as in 3D.

## List of fields in output files

When only the basic set of stats is enabled, only stats 1-11 are computed. When 2D stats are enabled \f$ \langle w \rangle \f$ is stored in s7 for basic stats and in s40 for the full set of statistics.

| Number | Statistic | Stored in variable (for fld files) |
| ------ | --------- | ------------------ |
| 1 | \f$ \langle p \rangle \f$ | Pressure|
| 2 | \f$ \langle u \rangle \f$ | X-Velocity|
| 3 | \f$ \langle v \rangle \f$ | Y-Velocity|
| 4 | \f$ \langle w \rangle \f$ | Z-Velocity |
| 5 | \f$ \langle pp \rangle \f$ | Temperature|
| 6 | \f$ \langle uu \rangle \f$ | Scalar 1 (s1)|
| 7 | \f$ \langle vv \rangle \f$ | Scalar 2 (s2)|
| 8 | \f$ \langle ww \rangle \f$ | s3|
| 9 | \f$ \langle uv \rangle \f$ | s4|
|10 | \f$ \langle uw \rangle \f$ | s5|
| 11 | \f$ \langle vw \rangle \f$ | s6|
| 12| \f$ \langle uuu \rangle \f$ | s7|
| 13| \f$ \langle vvv \rangle \f$ | s8|
| 14 | \f$ \langle www \rangle \f$ | s9|
| 15 | \f$ \langle  uuv   \rangle \f$ | s10 |
| 16 | \f$ \langle  uuw   \rangle \f$ | s11 |
| 17 | \f$ \langle  uvv   \rangle \f$ | s12 |
| 18 | \f$ \langle  uvw   \rangle \f$ | s13 |
| 19 | \f$ \langle  vvw   \rangle \f$ | s14 |
| 20 | \f$ \langle  uww   \rangle \f$ | s15 |
| 21 | \f$ \langle  vww   \rangle \f$ | s16 |
| 22 | \f$ \langle  uuuu  \rangle \f$ | s17 |
| 23 | \f$ \langle  vvvv  \rangle \f$ | s18 |
| 24 | \f$ \langle wwww   \rangle \f$ | s19 |
| 25 | \f$ \langle  ppp   \rangle \f$ | s20 |
| 26 | \f$ \langle  pppp  \rangle \f$ | s21 |
| 27 | \f$ \langle  pu    \rangle \f$ | s22 |
| 28 | \f$ \langle  pv    \rangle \f$ | s23 |
| 29 | \f$ \langle  pw    \rangle \f$ | s24 |
| 30 | \f$ \langle  p \frac{\partial u} {\partial x} \rangle \f$ | s25 |
| 31 | \f$ \langle  p \frac{\partial u} {\partial y}\rangle \f$ | s26 |
| 32 | \f$ \langle  p \frac{\partial u} {\partial z}\rangle \f$ | s27 |
| 33 | \f$ \langle  p \frac{\partial v} {\partial x}\rangle \f$ | s28 |
| 34 | \f$ \langle  p \frac{\partial v} {\partial y}\rangle \f$ | s29 |
| 35 | \f$ \langle  p \frac{\partial v} {\partial z}\rangle \f$ | s30 |
| 36 | \f$ \langle  p \frac{\partial w} {\partial x}\rangle \f$ | s31 |
| 37 | \f$ \langle  p \frac{\partial w} {\partial y}\rangle \f$ | s32 |
| 38 | \f$ \langle  p \frac{\partial w} {\partial z}\rangle \f$ | s33 |
| 39 | \f$ \langle  e11   \rangle \f$ | s34 |
| 40 | \f$ \langle  e22   \rangle \f$ | s35 |
| 41 | \f$ \langle  e33   \rangle \f$ | s36 |
| 42 | \f$ \langle  e12   \rangle \f$ | s37 |
| 43 | \f$ \langle  e13   \rangle \f$ | s38 |
| 44 | \f$ \langle  e23   \rangle \f$ | s39 |

where \f$e11,e22...\f$ is computed as:
$$
\begin{aligned}
e11 &= \left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial u}{\partial y}\right)^2 + \left(\frac{\partial u}{\partial z}\right)^2 \\\\
e22 &= \left(\frac{\partial v}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2 + \left(\frac{\partial v}{\partial z}\right)^2 \\\\
e33 &= \left(\frac{\partial w}{\partial x}\right)^2 + \left(\frac{\partial w}{\partial y}\right)^2 + \left(\frac{\partial w}{\partial z}\right)^2 \\\\
e12 &= \frac{\partial u}{\partial x}  \frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}\frac{\partial v}{\partial y}+ \frac{\partial u}{\partial z}\\frac{\partial v}{\partial z} \\\\
e13 &= \frac{\partial u}{\partial x}  \frac{\partial w}{\partial x} + \frac{\partial u}{\partial y}\frac{\partial w}{\partial y}+ \frac{\partial u}{\partial z}\\frac{\partial w}{\partial z} \\\\
e23 &= \frac{\partial v}{\partial x}  \frac{\partial w}{\partial x} + \frac{\partial v}{\partial y}\frac{\partial w}{\partial y}+ \frac{\partial v}{\partial z}\\frac{\partial w}{\partial z} \\\\
\end{aligned}
$$


# Postprocessing
These statistics are only the "raw statistics" in the sense that in general we are not interested in \f$ \langle uu\rangle \f$, but rather say \f$ \langle u'u'\rangle\f$. For this we need to postprocess the statistics. 

There is some rudimentary postprocessing to compute the spatial averages of fld files and also to combine the statistics collected from several runs (compute average in time) and also compute both the mean velocity gradient and the Reynolds stresses available in the contrib scripts. By running the contrib scripts without any arguments one gets a hint on their usage, (e.g. by running `./average_fields_in_time` will give the options). 

To further postprocess the statistics it is suggested to look into PyNekTools which introduces convenient functions for postprocessing, largely based on the PyMech library with the addition of an expanding set of tools for interpolation, computing derivatives and more advanced functionality.  PyNekTools is entirely parallelized in MPI and can also handle large data sets for postprocessing.


