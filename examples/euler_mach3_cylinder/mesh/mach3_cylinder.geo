// Structured quadrilateral mesh for the Ryujin Mach-3 cylinder benchmark.
//
// Reference geometry:
//   channel [0, 4] x [-1, 1]
//   cylinder centre (0.6, 0), diameter 0.5
//
// The transfinite counts below produce 25,020 quadrilateral elements.

SetFactory("OpenCASCADE");

cyl_x = 0.6;
cyl_y = 0.0;
cyl_rad = 0.25;
box_square = 0.42;
box_min_x = 0.0;
box_max_x = 4.0;
box_min_y = -1.0;
box_max_y = 1.0;

// Transfinite interval counts.
srf_nsplit = 40;
wsrf_nsplit = 40;
nsrf_nsplit = 27;
front_nsplit = 25;
wake_nsplit = 158;
spnw_nsplit = 30;

// Mild grading away from the cylinder and through the long wake.
wnprog = 1.02;
wwprog = 1.015;
wfprog = 1.0;

cs_el_sc = 0.01;

//////////////////////////////////////////////////////////////////////
// O-grid block around the cylinder
//////////////////////////////////////////////////////////////////////
pts_centre = newp;
Point(newp) = {cyl_x, cyl_y, 0.0, cs_el_sc};

x_pos = cyl_rad * Sin(Pi / 4.0);
pts_arc_1 = newp;
Point(newp) = {cyl_x - x_pos, cyl_y - x_pos, 0.0, cs_el_sc};
pts_arc_2 = newp;
Point(newp) = {cyl_x - x_pos, cyl_y + x_pos, 0.0, cs_el_sc};

x_pos = box_square;
pts_sqr_1 = newp;
Point(newp) = {cyl_x - x_pos, cyl_y - x_pos, 0.0, cs_el_sc};
pts_sqr_2 = newp;
Point(newp) = {cyl_x - x_pos, cyl_y + x_pos, 0.0, cs_el_sc};

list() = {newl};
Circle(newl) = {pts_arc_1, pts_centre, pts_arc_2};
list() += {newl};
Line(newl) = {pts_arc_2, pts_sqr_2};
list() += {newl};
Line(newl) = {pts_sqr_2, pts_sqr_1};
list() += {newl};
Line(newl) = {pts_sqr_1, pts_arc_1};

surface_list() = {};
crvl = newll;
Curve Loop(crvl) = {list()};
surface_list() += {news};
Surface(news) = {crvl};

langle = 0.5 * Pi;
surface_list() += Rotate {{0.0, 0.0, 1.0}, {cyl_x, cyl_y, 0.0}, langle} {
  Duplicata { Surface{surface_list(0)}; }
};
surface_list() += Rotate {{0.0, 0.0, 1.0}, {cyl_x, cyl_y, 0.0}, -langle} {
  Duplicata { Surface{surface_list(0)}; }
};
surface_list() += Rotate {{0.0, 0.0, 1.0}, {cyl_x, cyl_y, 0.0}, 2 * langle} {
  Duplicata { Surface{surface_list(0)}; }
};

Coherence;
list() = Unique(Abs(Boundary { Surface{surface_list()}; }));

//////////////////////////////////////////////////////////////////////
// Four side blocks
//////////////////////////////////////////////////////////////////////
ltmp[] = Extrude{box_min_x - (cyl_x - box_square), 0, 0} {
  Curve{list(1)};
};
surface_list() += {ltmp[1]};

ltmp[] = Extrude{box_max_x - (cyl_x + box_square), 0, 0} {
  Curve{list(10)};
};
surface_list() += {ltmp[1]};

ltmp[] = Extrude{0, box_min_y - (cyl_y - box_square), 0} {
  Curve{list(5)};
};
surface_list() += {ltmp[1]};

ltmp[] = Extrude{0, box_max_y - (cyl_y + box_square), 0} {
  Curve{list(7)};
};
surface_list() += {ltmp[1]};

Coherence;
list() = Unique(Abs(Boundary { Surface{surface_list()}; }));

//////////////////////////////////////////////////////////////////////
// Four corner blocks
//////////////////////////////////////////////////////////////////////
ltmp[] = Extrude{0, box_min_y - (cyl_y - box_square), 0} {
  Curve{list(13)};
};
surface_list() += {ltmp[1]};

ltmp[] = Extrude{0, box_max_y - (cyl_y + box_square), 0} {
  Curve{list(12)};
};
surface_list() += {ltmp[1]};

ltmp[] = Extrude{0, box_min_y - (cyl_y - box_square), 0} {
  Curve{list(15)};
};
surface_list() += {ltmp[1]};

ltmp[] = Extrude{0, box_max_y - (cyl_y + box_square), 0} {
  Curve{list(16)};
};
surface_list() += {ltmp[1]};

Coherence;
list() = Unique(Abs(Boundary { Surface{surface_list()}; }));

//////////////////////////////////////////////////////////////////////
// Physical zones used by mach3_cylinder.json
//////////////////////////////////////////////////////////////////////
Physical Surface("Fluid", 1) = {surface_list()};
Physical Curve("Inlet", 1) = {list(14), list(24), list(26)};
Physical Curve("Outlet", 2) = {list(17), list(28), list(30)};
Physical Curve("Top", 3) = {list(23), list(27), list(31)};
Physical Curve("Bottom", 4) = {list(20), list(25), list(29)};
Physical Curve("Cylinder", 7) = {list(3), list(6), list(9), list(11)};

//////////////////////////////////////////////////////////////////////
// Transfinite divisions
//////////////////////////////////////////////////////////////////////
Transfinite Curve{
  list(3), list(6), list(9), list(1), list(5), list(7),
  list(14), list(20), list(23)
} = (srf_nsplit + 1) Using Progression 1;

Transfinite Curve{list(11), list(10), list(17)} =
  (wsrf_nsplit + 1) Using Progression 1;

Transfinite Curve{
  list(18), list(19), list(21), list(22), list(24), list(26),
  list(28), list(30)
} = (spnw_nsplit + 1) Using Progression 1;

Transfinite Curve{list(12), list(13), list(25), list(27)} =
  (front_nsplit + 1) Using Progression wfprog;

Transfinite Curve{list(15), list(16), list(29), list(31)} =
  (wake_nsplit + 1) Using Progression wwprog;

Transfinite Curve{-list(0), list(2), -list(4), list(8)} =
  (nsrf_nsplit + 1) Using Progression wnprog;

Transfinite Surface{surface_list()};
Recombine Surface{surface_list()};

Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 0;
Mesh.Binary = 0;
Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 0;
