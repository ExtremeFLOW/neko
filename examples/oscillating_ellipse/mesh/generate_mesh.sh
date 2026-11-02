#!/bin/bash
set -e

# Set the paths to your executables here

GMSH_PATH="/Gmsh-Path/bin/gmsh"
GMSH2NMSH_PATH="/Neko-Path/bin/gmsh2nmsh"

"$GMSH_PATH" ellipse.geo -

# Physical surfaces 6 and 7 are made periodic
"$GMSH2NMSH_PATH" inclined_ellipse_3D.msh oscillating_ellipse.nmsh --periodic=6:7
cp oscillating_ellipse.nmsh ../

echo "Process completed successfully."
