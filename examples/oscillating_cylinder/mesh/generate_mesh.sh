#!/bin/bash
set -e

# Set the paths to your executables here

GMSH_PATH="/Gmsh-Path/bin/gmsh"
GMSH2NMSH_PATH="/Neko-Path/bin/gmsh2nmsh"

"$GMSH_PATH" cylinder.geo -

# Physical surfaces 3 and 4 are made periodic
"$GMSH2NMSH_PATH" 3D_ext_cyl.msh oscillating_cylinder.nmsh --periodic=3:4
cp oscillating_cylinder.nmsh ../

echo "Process completed successfully."
