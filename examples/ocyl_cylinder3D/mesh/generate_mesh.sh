#!/bin/bash

gmsh cylinder.geo -

gmsh2nek <<EOF
3
3D_ext_cyl
0
1
3 4
cylinder
EOF

rea2nbin cylinder.re2 cylinder.nmsh
cp cylinder.nmsh ../
