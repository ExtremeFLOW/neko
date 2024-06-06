# Statistics guide {#statistics-guide}

\tableofcontents

Under development, updated incrementally

Statistics in the context of Neko, is the common name for fields that are averaged in time and possible also space.

The statistics module in Neko computes the temporal average of a wide range of fields.

In the following page we use the following convetion for a field
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
Statistics are enable in the the case file as the following:

| Name                | Description                                                          | Admissible values | Default value |
| ------------------- | -------------------------------------------------------------------- | ----------------- | ------------- |
| `enabled`           | Whether to enable the statistics computation.                        | `true` or `false` | `true`        |
| `start_time`        | Time at which to start gathering statistics.                         | Positive real     | 0             |
| `sampling_interval` | Interval, in timesteps, for sampling the flow fields for statistics. | Positive integer  | 10            |

In addition to the usual controls for the output, which then outputs the averages computes from the last time the statistics were written to file.

For example, if one wants to sample the fields every 4 time steps and compute the averages in time intervals of 20 and write the output every 20 time units, and start collecting statistics after an initial transient of 50 time units the following would work:

~~~~~~~~~~~~~~~{.json}
"statistics": {
        "enabled": true,
        "start_time": 50.0,
        "sampling_interval": 4,
        "output_control": "simulationtime",
        "output_value": 20,
  }
~~~~~~~~~~~~~~~
When the output is written one obtains two .fld files called mean_field and stats. 

## List of fields in output files 
In `mean_field` the following averages are stored. The stored in variable column is which field one finds the computed statistic if one opens the file in paraview or visit.

| Number | Statistic | Stored in variable |
| ------ | --------- | ------------------ |
| 1 | \f$ \langle p \rangle \f$ | Pressure|
| 2 | \f$ \langle u \rangle \f$ | X-Velocity|
| 3 | \f$ \langle v \rangle \f$ | Y-Velocity|
| 4 | \f$ \langle w \rangle \f$ | Z-Velocity|

In `stats` several other statistics are stored, and while not all might be interesting to your specific use case, with them most different budgets and quantities of interest can be computed. They are stored as the following:


| Number | Statistic | Stored in variable |
| ------ | --------- | ------------------ |
| 1 | \f$ \langle pp \rangle \f$ | Pressure|
| 2 | \f$ \langle uu \rangle \f$ | X-Velocity|
| 3 | \f$ \langle vv \rangle \f$ | Y-Velocity|
| 4 | \f$ \langle ww \rangle \f$ | Z-Velocity|
| 5 | \f$ \langle uv \rangle \f$ | Temperature|
| 6 | \f$ \langle uw \rangle \f$ | Scalar 1 (s1)|
| 7 | \f$ \langle vw \rangle \f$ | Scalar 2 (s2)|
| 8 | \f$ \langle uuu \rangle \f$ | s3|
| 9 | \f$ \langle vvv \rangle \f$ | s4|
| 10 | \f$ \langle www \rangle \f$ | s5|
| 11 | \f$ \langle  uuv   \rangle \f$ | s6 |
| 12 | \f$ \langle  uuw   \rangle \f$ | s7 |
| 13 | \f$ \langle  uvv   \rangle \f$ | s8 |
| 14 | \f$ \langle  uvw   \rangle \f$ | s9 |
| 15 | \f$ \langle  vvw   \rangle \f$ | s10 |
| 16 | \f$ \langle  uww   \rangle \f$ | s11 |
| 17 | \f$ \langle  vww   \rangle \f$ | s12 |
| 18 | \f$ \langle  uuuu  \rangle \f$ | s13 |
| 19 | \f$ \langle  vvvv  \rangle \f$ | s14 |
| 20 | \f$ \langle wwww   \rangle \f$ | s15 |
| 21 | \f$ \langle  ppp   \rangle \f$ | s16 |
| 22 | \f$ \langle  pppp  \rangle \f$ | s17 |
| 23 | \f$ \langle  pu    \rangle \f$ | s18 |
| 24 | \f$ \langle  pv    \rangle \f$ | s19 |
| 25 | \f$ \langle  pw    \rangle \f$ | s20 |
| 26 | \f$ \langle  p \frac{\partial u} {\partial x} \rangle \f$ | s21 |
| 27 | \f$ \langle  p \frac{\partial u} {\partial y}\rangle \f$ | s22 |
| 28 | \f$ \langle  p \frac{\partial u} {\partial z}\rangle \f$ | s23 |
| 29 | \f$ \langle  p \frac{\partial v} {\partial x}\rangle \f$ | s24 |
| 30 | \f$ \langle  p \frac{\partial v} {\partial y}\rangle \f$ | s25 |
| 31 | \f$ \langle  p \frac{\partial v} {\partial z}\rangle \f$ | s26 |
| 32 | \f$ \langle  p \frac{\partial w} {\partial x}\rangle \f$ | s27 |
| 33 | \f$ \langle  p \frac{\partial w} {\partial y}\rangle \f$ | s28 |
| 34 | \f$ \langle  p \frac{\partial w} {\partial z}\rangle \f$ | s29 |
| 35 | \f$ \langle  e11   \rangle \f$ | s30 |
| 36 | \f$ \langle  e22   \rangle \f$ | s31 |
| 37 | \f$ \langle  e33   \rangle \f$ | s32 |
| 38 | \f$ \langle  e12   \rangle \f$ | s33 |
| 39 | \f$ \langle  e13   \rangle \f$ | s34 |
| 40 | \f$ \langle  e23   \rangle \f$ | s35 |

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
Of course, these statistics are only the "raw statistics" in the sense that in general we are not interested in \f$ \langle uu\rangle \f$, but rather say the rms of the velocity fluctuation. FOr this we need to postprocess the statistics. 

There is some rudimentary postprocessing to compute the spatial averages of fld filesa and also to combine the statistics collected from several runs (compute average in time) and also compute both the rms values and the Reynolds stresses available among the contrib scripts. By running the contrib scripts without any arguments one gets a hint on their usage, and also the text below gives a guide on how to postprocess the raw statistics. The postprocessing part of Neko is expanding and changing quite a lot at the moment, where we currently envision primarily using python for the postprocessing of the final statistics.



## Notes on the statistics calculation in Neko
***Daniele Massaro, Martin Karp (KTH)***

1) Run your simulations and collect mean_field* and stats* files by having the statistics object added to the case file,  and specifying the write interval to something suitable.

2) For each RUN_i, you get a set of mean_field* and stats* files. You can average them for each single RUN_i, or average all of them only once (after re-ordering them properly). If you follow the second approach, go to step 4. 
Here, for each RUN_i, we compute the averaged means with "average_fields_in_time":
--mean
`srun --unbuffered /your/location/neko/bin/average_fields_in_time meanXX.fld T0 mean_p.fld`
where T0 is the initial time. To get some hints on the input for the script one can simply run `./average_fields_in_time` without any arguments. For RUN_1 the time T0 can be taken from the log of the first simulation, or from the header of the first mean_field* file; in this way you discard that file. For RUN_i, with i>1, it can be taken from header of the last file mean_field* of the previous simulation RUN_{i-1}. 
In the command line, for the name "meanXX.fld", XX indicates the number of the nek5000 file. In mean_fieldXX.nek5000 you set the number of the first mean0* file to read and the number of steps corresponding to the number of files. In this way, the code generates a mean_p0.f00000 and mean_post0.nek5000. It is suggested to rename mean_p0.f00000 as mean_p0.f0000i and move it to a separate folder where you take the average with all the others. 
--stats
`srun --unbuffered /your/location/neko/bin/average_fields_in_time statXX.fld T0 stat_p.fld`
T0 is the same as before. In stat0.nek5000 you set the number of the first stat0* file to read and the number of steps corresponds to the number of files. It is suggested to rename stat_p0.f00000 as stat_p0.f0000i and move it to a separate folder where you take the average with all the others. 
Repeat this for each RUN_i folder. Eventually, given n RUN_i folders, you will get n mean_p* and stat_p* files.

3) Take the average of the averaged runs. Now, the time average over all the n simulations is taken. The procedure is similar, but changing the output name is recommended to  avoid over-writing.
-- mean
`srun --unbuffered /your/location/neko/bin/average_fields mean_p0.fld T0 mean_post.fld`
where T0 is the initial time which has been used to compute mean_p* for RUN_1.
-- stats
`srun --unbuffered /your/location/neko/bin/average_fields stat_p0.fld T0 stat_post.fld`
where T0 is the initial time which has been used to compute mean_p* for RUN_1.





4) Compute Reynolds stress tensors and other statistical moments (see the list).
`srun --unbuffered /your/location/neko/bin/postprocess_fluid_stats mesh_file.nmsh mean_post0.fld stat_post0.fld`

5) We also provide a tool to average the resulting field in a homogenous direction in `bin/average_field_in_space`. The required arguments are shown if one runs the program without any input. Currently it requires the number of elements in the homogenous direction as an input argument, e.g. 
`./average_field_in_space mesh.nmsh field.fld x 18 outfield.fld`
if we want to average a field in the x direction on a mesh with 18 elements in x and output the averaged field in outfield0.nek5000.







