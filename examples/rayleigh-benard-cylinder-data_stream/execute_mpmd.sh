#!/bin/bash

mpirun -n 1 ./neko rayleigh.case : -n 2 python3 run_spMoDePy.py
