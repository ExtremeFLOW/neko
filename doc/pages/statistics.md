
# Statistics guide {#statistics-guide}

Under development, updated incrementally
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







