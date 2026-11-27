// 4 x 3 x 2 box with one physical surface per side.
// gmsh box.geo -3 -format msh41 -o box_v41.msh
// gmsh box.geo -3 -format msh22 -bin -o box_v22_bin.msh
Point(1) = {0, 0, 0, 1}; Point(2) = {2, 0, 0, 1};
Point(3) = {2, 1.5, 0, 1}; Point(4) = {0, 1.5, 0, 1};
Point(5) = {0, 0, 1, 1}; Point(6) = {2, 0, 1, 1};
Point(7) = {2, 1.5, 1, 1}; Point(8) = {0, 1.5, 1, 1};
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
Line(5) = {5, 6}; Line(6) = {6, 7}; Line(7) = {7, 8}; Line(8) = {8, 5};
Line(9) = {1, 5}; Line(10) = {2, 6}; Line(11) = {3, 7}; Line(12) = {4, 8};
Curve Loop(1) = {4, 9, -8, -12}; Plane Surface(1) = {1};
Curve Loop(2) = {2, 11, -6, -10}; Plane Surface(2) = {2};
Curve Loop(3) = {1, 10, -5, -9}; Plane Surface(3) = {3};
Curve Loop(4) = {3, 12, -7, -11}; Plane Surface(4) = {4};
Curve Loop(5) = {1, 2, 3, 4}; Plane Surface(5) = {5};
Curve Loop(6) = {5, 6, 7, 8}; Plane Surface(6) = {6};
Surface Loop(1) = {1, 2, 3, 4, 5, 6}; Volume(1) = {1};
Transfinite Curve{1, 3, 5, 7} = 5;
Transfinite Curve{2, 4, 6, 8} = 4;
Transfinite Curve{9, 10, 11, 12} = 3;
Transfinite Surface{:}; Recombine Surface{:}; Transfinite Volume{1};
Physical Surface("xmin", 1) = {1}; Physical Surface("xmax", 2) = {2};
Physical Surface("ymin", 3) = {3}; Physical Surface("ymax", 4) = {4};
Physical Surface("zmin", 5) = {5}; Physical Surface("zmax", 6) = {6};
Physical Volume("fluid", 7) = {1};
