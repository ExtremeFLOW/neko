#!/bin/bash

gmsh ellipse.geo -

gmsh2nek <<EOF
3
inclined_ellipse_3D
0
1
6 7
ellipse
EOF

rea2nbin ellipse.re2 ellipse.nmsh
cp ellipse.nmsh ../
