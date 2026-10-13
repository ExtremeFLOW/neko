#!/bin/bash
set -e

# Set the paths to your executables here

GMSH_PATH="/Gmsh-Path/bin/gmsh"
GMSH2NMSH_PATH="/Neko-Path/bin/gmsh2nmsh"

"$GMSH_PATH" double_cylinder.geo -

# Physical surfaces 9 and 8 are made periodic
"$GMSH2NMSH_PATH" 3D_ext_cyl.msh double_oscillating_cylinders.nmsh --periodic=9:8
cp double_oscillating_cylinders.nmsh ../

echo "Process completed successfully."
