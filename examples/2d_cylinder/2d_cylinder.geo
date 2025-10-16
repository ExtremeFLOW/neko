//
// Block structured mesh for 2D flow around cylinder
//
//////////////////////////////////////////////////////////////////////
// Parameter section
//////////////////////////////////////////////////////////////////////
// I use open cascade
SetFactory("OpenCASCADE");

// mesh parameters
cyl_rad = DefineNumber[ 0.5 , Min 0.25, Max 1.0, Step 0.05, Name "Cylinder radius" ];
box_square = DefineNumber[ 3.0 , Min 2.0, Max 10.0, Step 1.0, Name "Size of a square region around cylinder" ];
box_min_x = DefineNumber[ -15.0 , Min -30.0, Max -10.0, Step 1.0, Name "Domain minimum x" ];
box_max_x = DefineNumber[  35.0 , Min  10.0, Max  50.0, Step 1.0, Name "Domain maximum x" ];
box_width = DefineNumber[ 30.0 , Min 20.0, Max 50.0, Step 1.0, Name "Domain with" ];
// mesh parameters
// number of splits at the cylinder surface
srf_nsplit = DefineNumber[ 4 , Min 2, Max 20, Step 1, Name "Numer of splits; cylinder surface" ];
// number of splits at the cylinder surface from the wake side
wsrf_nsplit = DefineNumber[ 6 , Min 2, Max 20, Step 1, Name "Numer of splits; cylinder surface from the wake side" ];
// number of splits normal to the cylinder surface
nsrf_nsplit = DefineNumber[ 6 , Min 4, Max 20, Step 1, Name "Numer of splits; normal to the cylinder surface" ];
// number of splits in front of cylinder
front_nsplit = DefineNumber[ 6 , Min 2, Max 20, Step 1,   Name "Numer of splits; front of cylinder" ];
// number of splits in wake
wake_nsplit = DefineNumber[ 15 , Min 2, Max 50, Step 1,   Name "Numer of splits; wake" ];
// spanwise number of splits
spnw_nsplit = DefineNumber[ 10 , Min 2, Max 50, Step 1,   Name "Spanwise numer of splits" ];

// progression normal to the cylinder surface
wnprog = DefineNumber[ 1.4 , Min 1, Max 2, Step 0.1,   Name "Progression normal to the cylinder surface" ];
// progression in the wake region
wwprog = DefineNumber[ 1.05 , Min 1, Max 2, Step 0.01,   Name "Progression in the wake region" ];
// progression in the front region
wfprog = DefineNumber[ 1.08 , Min 1, Max 2, Step 0.01,   Name "Progression in the front region" ];

// Element scale at cylinder surface
cs_el_sc = 0.1;

//////////////////////////////////////////////////////////////////////
// Square region around cylinder
//////////////////////////////////////////////////////////////////////
// Points
pts_centre = newp;
Point(newp) = {0.0,0.0,0.0,cs_el_sc};

x_pos = cyl_rad*Sin(Pi/4.0);
pts_arc_1 = newp;
Point(newp) = {-x_pos,-x_pos,0.0,cs_el_sc};
pts_arc_2 = newp;
Point(newp) = {-x_pos,x_pos,0.0,cs_el_sc};

x_pos = box_square;
pts_sqr_1 = newp;
Point(newp) = {-x_pos,-x_pos,0.0,cs_el_sc};
pts_sqr_2 = newp;
Point(newp) = {-x_pos,x_pos,0.0,cs_el_sc};

// Curves
list() = {newl};
Circle(newl) = {pts_arc_1,pts_centre,pts_arc_2};
list() += {newl};
Line(newl) = {pts_arc_2,pts_sqr_2};
list() += {newl};
Line(newl) = {pts_sqr_2,pts_sqr_1};
list() += {newl};
Line(newl) = {pts_sqr_1,pts_arc_1};

// Surfaces
surface_list() = {};
crvl = newll;
Curve Loop (crvl) = {list()};
surface_list() += {news};
Surface(news) = {crvl};


// Replicate rotate
langle = 0.5*Pi;
surface_list() += Rotate {{0.0,0.0,1.0},{0.0,0.0,0.0},langle} { Duplicata{ Surface{surface_list(0)}; }};
surface_list() += Rotate {{0.0,0.0,1.0},{0.0,0.0,0.0},-langle} { Duplicata{ Surface{surface_list(0)}; }};
surface_list() += Rotate {{0.0,0.0,1.0},{0.0,0.0,0.0},2*langle} { Duplicata{ Surface{surface_list(0)}; }};

// Remove duplicates
Coherence;

// Extract edges
list() = Unique(Abs(Boundary { Surface{surface_list()}; }));
//For il In {0:# list() - 1}
//   Printf("list %g %g", il,list(il));
//EndFor

// Extrude
// -x
ltmp[] = Extrude{box_square+box_min_x,0,0} { Curve{list(1)}; };
surface_list() += {ltmp[1]};
// +x
ltmp[] = Extrude{box_max_x-box_square,0,0} { Curve{list(10)}; };
surface_list() += {ltmp[1]};
// -y
ltmp[] = Extrude{0,box_square-box_width*0.5,0,0} { Curve{list(5)}; };
surface_list() += {ltmp[1]};
// +y
ltmp[] = Extrude{0,box_width*0.5-box_square,0,0} { Curve{list(7)}; };
surface_list() += {ltmp[1]};

// Remove duplicates
Coherence;

// Extract edges
list() = Unique(Abs(Boundary { Surface{surface_list()}; }));
//For il In {0:# list() - 1}
//   Printf("list %g %g", il,list(il));
//EndFor

// Extrude
// -x -y
ltmp[] = Extrude{0,box_square-box_width*0.5,0,0} { Curve{list(13)}; };
surface_list() += {ltmp[1]};
// -x +y
ltmp[] = Extrude{0,box_width*0.5-box_square,0,0} { Curve{list(12)}; };
surface_list() += {ltmp[1]};
// +x -y
ltmp[] = Extrude{0,box_square-box_width*0.5,0,0} { Curve{list(15)}; };
surface_list() += {ltmp[1]};
// +x +y
ltmp[] = Extrude{0,box_width*0.5-box_square,0,0} { Curve{list(16)}; };
surface_list() += {ltmp[1]};

// Remove duplicates
Coherence;

// Extract edges
list() = Unique(Abs(Boundary { Surface{surface_list()}; }));
//For il In {0:# list() - 1}
//   Printf("list %g %g", il,list(il));
//EndFor

//////////////////////////////////////////////////////////////////////
// Physical properties section
//////////////////////////////////////////////////////////////////////
// Volume
Physical Surface("Fluid",1) = {surface_list()};

// Wall
Physical Curve("Wal",2) = {list(3),list(6),list(9),list(11)};

// Inflow
Physical Curve("Inf",3) = {list(14),list(24),list(26)};

// Outflow
Physical Curve("Out",4) = {list(17),list(28),list(30)};

// Top periodic
Physical Curve("Pet",5) = {list(23),list(27),list(31)};

// Bottom periodic
Physical Curve("Peb",6) = {list(20),list(25),list(29)};

//////////////////////////////////////////////////////////////////////
// Transfinite division section
//////////////////////////////////////////////////////////////////////
// Edge division
//
Transfinite Curve{list(3),list(6),list(9),list(1),list(5),list(7),list(14),list(20),list(23)} = (srf_nsplit+1) Using Progression 1;

Transfinite Curve{list(11),list(10),list(17)} = (wsrf_nsplit+1) Using Progression 1;

Transfinite Curve{list(18),list(19),list(21),list(22),list(24),list(26),list(28),list(30)} = (spnw_nsplit+1) Using Progression 1;

Transfinite Curve{list(12),list(13),list(25),list(27)} = (front_nsplit+1) Using Progression wfprog;

Transfinite Curve{list(15),list(16),list(29),list(31)} = (wake_nsplit+1) Using Progression wwprog;

Transfinite Curve{-list(0),list(2),-list(4),list(8)} = (nsrf_nsplit+1) Using Progression wnprog;

// Surface division
Transfinite Surface{surface_list()};
Recombine Surface{surface_list()};


//////////////////////////////////////////////////////////////////////
// Meshing section
//////////////////////////////////////////////////////////////////////
Mesh 1;
Mesh 2;
Mesh 3;

SetOrder 2;

RenumberMeshElements;

//////////////////////////////////////////////////////////////////////
// Mesh saving section
//////////////////////////////////////////////////////////////////////
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 0;
Mesh.Binary = 0;

Save "ext_cyl.msh";

