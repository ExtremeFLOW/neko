#!/bin/bash

genmeshbox 0 0.5 0 0.1 0 0.1 10 1 1 .false. .true. .true.

# using csv as input and output
# python3 particles_init.py
# neko rebounding_particles_csv.case
# python3 trajectory_plot_csv.py

# using json as input and hdf5 as output
neko rebounding_particles_h5.case
python3 trajectory_plot_h5.py
