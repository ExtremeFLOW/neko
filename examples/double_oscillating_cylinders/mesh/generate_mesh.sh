#!/bin/bash

# Set the paths to your executables here

GMSH_PATH="/Gmsh-Path/"
GMSH2NEK_PATH="/Neko-Path/contrib/gmsh2nek/"
REA2NBIN_PATH="/Neko-Path/bin/rea2nbin"

"$GMSH_PATH" double_cylinder.geo -

"$GMSH2NEK_PATH" <<EOF
3
3D_ext_cyl
1
9 8
0.0 0.0 -0.50
3D_ext_cyl
EOF

"$REA2NBIN_PATH" 3D_ext_cyl.re2 double_oscillating_cylinders.nmsh
cp double_oscillating_cylinders.nmsh ../

echo "Process completed successfully."

