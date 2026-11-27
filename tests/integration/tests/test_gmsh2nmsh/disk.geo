// 2D annulus in the x-y plane, 16 second order quadrilaterals. One sector
// is meshed clockwise and the physical tags are above 20.
// gmsh disk.geo -2 -order 2 -format msh22 -o disk_v22.msh
R1 = 0.5; R2 = 1.0;
Point(1) = {0, 0, 0};
Point(2) = {R1, 0, 0}; Point(3) = {0, R1, 0}; Point(4) = {-R1, 0, 0}; Point(5) = {0, -R1, 0};
Point(6) = {R2, 0, 0}; Point(7) = {0, R2, 0}; Point(8) = {-R2, 0, 0}; Point(9) = {0, -R2, 0};
Circle(1) = {2, 1, 3}; Circle(2) = {3, 1, 4}; Circle(3) = {4, 1, 5}; Circle(4) = {5, 1, 2};
Circle(5) = {6, 1, 7}; Circle(6) = {7, 1, 8}; Circle(7) = {8, 1, 9}; Circle(8) = {9, 1, 6};
Line(9) = {2, 6}; Line(10) = {3, 7}; Line(11) = {4, 8}; Line(12) = {5, 9};
Curve Loop(1) = {9, 5, -10, -1}; Plane Surface(1) = {1};
// Reversed loop: this sector gets clockwise (left-handed) quads
Curve Loop(2) = {2, 11, -6, -10}; Plane Surface(2) = {2};
Curve Loop(3) = {11, 7, -12, -3}; Plane Surface(3) = {3};
Curve Loop(4) = {12, 8, -9, -4}; Plane Surface(4) = {4};
Transfinite Curve{1:8} = 3; Transfinite Curve{9:12} = 3;
Transfinite Surface{1:4}; Recombine Surface{1:4};
Physical Curve("inner", 101) = {1, 2, 3, 4};
Physical Curve("outer", 102) = {5, 6, 7, 8};
Physical Surface("fluid", 103) = {1, 2, 3, 4};
