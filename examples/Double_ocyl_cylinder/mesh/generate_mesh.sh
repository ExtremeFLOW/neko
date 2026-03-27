#!/bin/bash

gmsh double_cylinder.geo -

gmsh2nek <<EOF
3
3D_ext_cyl
0
1
9 8
cylinder
EOF

rea2nbin cylinder.re2 cylinder.nmsh
cp cylinder.nmsh ../
