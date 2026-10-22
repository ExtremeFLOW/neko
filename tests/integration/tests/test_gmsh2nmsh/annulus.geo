// Annular pipe section with curved walls, 16 second order hexahedra.
// gmsh annulus.geo -3 -order 2 -format msh41 -bin -o annulus_v41_bin.msh
R1 = 0.5; R2 = 1.0; L = 1.0;
Point(1) = {0, 0, 0};
Point(2) = {R1, 0, 0}; Point(3) = {0, R1, 0}; Point(4) = {-R1, 0, 0}; Point(5) = {0, -R1, 0};
Point(6) = {R2, 0, 0}; Point(7) = {0, R2, 0}; Point(8) = {-R2, 0, 0}; Point(9) = {0, -R2, 0};
Circle(1) = {2, 1, 3}; Circle(2) = {3, 1, 4}; Circle(3) = {4, 1, 5}; Circle(4) = {5, 1, 2};
Circle(5) = {6, 1, 7}; Circle(6) = {7, 1, 8}; Circle(7) = {8, 1, 9}; Circle(8) = {9, 1, 6};
Line(9) = {2, 6}; Line(10) = {3, 7}; Line(11) = {4, 8}; Line(12) = {5, 9};
Curve Loop(1) = {9, 5, -10, -1}; Plane Surface(1) = {1};
Curve Loop(2) = {10, 6, -11, -2}; Plane Surface(2) = {2};
Curve Loop(3) = {11, 7, -12, -3}; Plane Surface(3) = {3};
Curve Loop(4) = {12, 8, -9, -4}; Plane Surface(4) = {4};
Transfinite Curve{1:8} = 3; Transfinite Curve{9:12} = 3;
Transfinite Surface{1:4}; Recombine Surface{1:4};
out[] = Extrude {0, 0, L} { Surface{1:4}; Layers{1}; Recombine; };
Physical Surface("inner", 1) = {out[5], out[11], out[17], out[23]};
Physical Surface("outer", 2) = {out[3], out[9], out[15], out[21]};
Physical Surface("bottom", 3) = {1, 2, 3, 4};
Physical Surface("top", 4) = {out[0], out[6], out[12], out[18]};
Physical Volume("fluid", 10) = {out[1], out[7], out[13], out[19]};
