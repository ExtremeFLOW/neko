#!/bin/bash

# Set the paths to your executables here

GMSH_PATH="/Gmsh-Path/bin/gmsh"
GMSH2NEK_PATH="/Neko-Path/contrib/gmsh2nek/gmsh2nek"
REA2NBIN_PATH="/Neko-Path/bin/rea2nbin"

"$GMSH_PATH" ellipse.geo -

"$GMSH2NEK_PATH" <<EOF
3
inclined_ellipse_3D
1
6 7
0.0 0.0 0.50
inclined_ellipse_3D
EOF

"$REA2NBIN_PATH" inclined_ellipse_3D.re2 oscillating_ellipse.nmsh
cp oscillating_ellipse.nmsh ../

echo "Process completed successfully."