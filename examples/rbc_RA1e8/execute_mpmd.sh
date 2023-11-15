#!/bin/bash
mpirun -n 1 ./neko rayleigh.case : -n 4 python3 insitu/run_spMoDePy.py
