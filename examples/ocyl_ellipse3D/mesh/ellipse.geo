/*
   -----------  gmsh script for generating 2D/3D inclined ellipse --------------
   -----------    Best mesh quality for theta \approx 45 degree   --------------
   -----------      For circle: set theta=45, and a=b=radious     --------------
   By Mohammad Moniripiri, momp@kth.se
   FLOW, KTH Mechanics, Stockholm, Sweden

   Based on the code by Daniele Massaro (cylinder example, Neko)

*/

meshDim = 3;

// =================================================
// MESH CONTROL
// =================================================
Nxc  = 15;   Rxc  = 1.00;
Nxw  = 10;   Rxw  = 1.22;
Nxlr = 7;    Rxlr = 1.0;
Nxll = 6;    Rxll = 1.0;

Nyc  = Nxc;  Ryc  = 1.0;
Nytb = 5;    Rytb = 1.0;
Nytp = 5;    Rytp = 1.0;

Nb   = 15;   Rb   = 0.9;
Nc   = Nxc;  Rc   = 1.0;

// =================================================
// GEOMETRY PARAMETERS
// =================================================

// Internal box
xc = 2;
yc = 2;

// Ellipse
a = 0.5;                  // semi-major axis
b = 0.25;                  // semi-minor axis
theta = 45*Pi/180;        // inclination

// External box
xe = 3.5;
ye = 4;

// Wake
xw = 10;

// 3D
If (meshDim == 3)
  Lz = 0.5;
  Nz = 1;
EndIf

// =================================================
// POINTS
// =================================================

// INTERNAL BOX
Point(1) = {-xc, -yc, 0};
Point(2) = { xc, -yc, 0};
Point(3) = {-xc,  yc, 0};
Point(4) = { xc,  yc, 0};

// ELLIPSE AXIS POINTS (before rotation)
Point(5) = {-a,  0, 0};
Point(6) = { 0, -b, 0};
Point(7) = { a,  0, 0};
Point(8) = { 0,  b, 0};
Point(9) = { 0,  0, 0}; // center

// Rotate ellipse points
Rotate {{0,0,1},{0,0,0},theta} { Point{5:8}; }

// EXTERNAL BOX
Point(10) = {-xe, -ye, 0};
Point(11) = { xe, -ye, 0};
Point(12) = {-xe,  ye, 0};
Point(13) = { xe,  ye, 0};

Point(14) = {-xe, -yc, 0};
Point(15) = { xe, -yc, 0};
Point(16) = {-xe,  yc, 0};
Point(17) = { xe,  yc, 0};

Point(22) = { xc, -ye, 0};
Point(23) = {-xc, -ye, 0};
Point(24) = { xc,  ye, 0};
Point(25) = {-xc,  ye, 0};

// WAKE
Point(18) = {xw,  ye, 0};
Point(19) = {xw, -ye, 0};
Point(20) = {xw,  yc, 0};
Point(21) = {xw, -yc, 0};

// =================================================
// LINES
// =================================================

// INTERNAL BOX
Line(1) = {1,2}; Transfinite Line {1} = Nxc Using Progression Rxc;
Line(2) = {3,4}; Transfinite Line {2} = Nxc Using Progression Rxc;
Line(3) = {1,3}; Transfinite Line {3} = Nyc Using Progression Ryc;
Line(4) = {2,4}; Transfinite Line {4} = Nyc Using Progression Ryc;

// ELLIPSE (native)
Ellipse(5) = {5,9,6}; Transfinite Line {5} = Nc Using Progression Rc;
Ellipse(6) = {6,9,7}; Transfinite Line {6} = Nc Using Progression Rc;
Ellipse(7) = {7,9,8}; Transfinite Line {7} = Nc Using Progression Rc;
Ellipse(8) = {8,9,5}; Transfinite Line {8} = Nc Using Progression Rc;

// DIAGONALS
Line(9)  = {1,5}; Transfinite Line {9}  = Nb Using Progression Rb;
Line(10) = {2,6}; Transfinite Line {10} = Nb Using Progression Rb;
Line(11) = {4,7}; Transfinite Line {11} = Nb Using Progression Rb;
Line(12) = {3,8}; Transfinite Line {12} = Nb Using Progression Rb;

// RIGHT
Line(17) = {4,17}; Transfinite Line {17} = Nxlr Using Progression Rxlr;
Line(18) = {2,15}; Transfinite Line {18} = Nxlr Using Progression Rxlr;
Line(19) = {15,17};Transfinite Line {19} = Nyc Using Progression Ryc;

// LEFT
Line(20) = {1,14}; Transfinite Line {20} = Nxll Using Progression Rxll;
Line(21) = {3,16}; Transfinite Line {21} = Nxll Using Progression Rxll;
Line(22) = {14,16};Transfinite Line {22} = Nyc Using Progression Ryc;

// BOTTOM
Line(23) = {15,11}; Transfinite Line {23} = Nytb Using Progression Rytb;
Line(24) = {15,14}; Transfinite Line {24} = Nxw Using Progression Rxw;
Line(25) = {14,10}; Transfinite Line {25} = Nytb Using Progression Rytb;

Line(40) = {1,23}; Transfinite Line {40} = Nytb Using Progression Rytb;
Line(41) = {2,22}; Transfinite Line {41} = Nytb Using Progression Rytb;
Line(42) = {23,22};Transfinite Line {42} = Nxc Using Progression Rxc;
Line(43) = {23,10};Transfinite Line {43} = Nxll Using Progression Rxll;
Line(44) = {22,11};Transfinite Line {44} = Nxlr Using Progression Rxlr;

// TOP
Line(26) = {16,17}; Transfinite Line {26} = Nxw Using Progression Rxw;
Line(27) = {13,17}; Transfinite Line {27} = Nytp Using Progression Rytp;
Line(28) = {12,16}; Transfinite Line {28} = Nytp Using Progression Rytp;

Line(45) = {25,3}; Transfinite Line {45} = Nytp Using Progression Rytp;
Line(46) = {25,24};Transfinite Line {46} = Nxc Using Progression Rxc;
Line(47) = {24,4}; Transfinite Line {47} = Nytp Using Progression Rytp;
Line(48) = {25,12};Transfinite Line {48} = Nxll Using Progression Rxll;
Line(49) = {24,13};Transfinite Line {49} = Nxlr Using Progression Rxlr;

// WAKE
Line(32) = {11,19}; Transfinite Line {32} = Nxw Using Progression Rxw;
Line(34) = {13,18}; Transfinite Line {34} = Nxw Using Progression Rxw;
Line(35) = {17,20}; Transfinite Line {35} = Nxw Using Progression Rxw;
Line(36) = {15,21}; Transfinite Line {36} = Nxw Using Progression Rxw;
Line(37) = {21,20}; Transfinite Line {37} = Nyc Using Progression Ryc;
Line(38) = {18,20}; Transfinite Line {38} = Nytp Using Progression Rytp;
Line(39) = {21,19}; Transfinite Line {39} = Nytb Using Progression Rytb;

// =================================================
// SURFACES
// =================================================

Curve Loop(1)  = {10,6,-11,-4};   Plane Surface(1)  = {1};
Curve Loop(2)  = {11,7,-12,2};    Plane Surface(2)  = {2};
Curve Loop(3)  = {12,8,-9,3};     Plane Surface(3)  = {3};
Curve Loop(4)  = {5,-10,-1,9};    Plane Surface(4)  = {4};
Curve Loop(5)  = {4,17,-19,-18};  Plane Surface(5)  = {5};
Curve Loop(6)  = {-3,20,22,-21};  Plane Surface(6)  = {6};
Curve Loop(7)  = {40,43,-25,-20}; Plane Surface(7)  = {7};
Curve Loop(8)  = {41,-42,-40,1};  Plane Surface(8)  = {8};
Curve Loop(9)  = {-23,-18,41,44}; Plane Surface(9)  = {9};
Curve Loop(10) = {-27,-49,47,17}; Plane Surface(10) = {10};
Curve Loop(11) = {47,-2,-45,46};  Plane Surface(11) = {11};
Curve Loop(12) = {-48,45,21,-28}; Plane Surface(12) = {12};
Curve Loop(13) = {-27,34,38,-35}; Plane Surface(13) = {13};
Curve Loop(14) = {19,35,-37,-36}; Plane Surface(14) = {14};
Curve Loop(15) = {-23,36,39,-32}; Plane Surface(15) = {15};

Transfinite Surface "*";
Recombine Surface "*";

// =================================================
// MESH GENERATION (2D)
// =================================================
If (meshDim == 2)
   Physical Curve("wall") = {5,6,7,8};
   Physical Curve("inlet") = {22,28,25};
   Physical Curve("outlet") = {37,38,39};
   Physical Curve("top") = {48,46,49,34};
   Physical Curve("bottom") = {43,42,44,32};
   Physical Surface("fluid") = {1:15};

   Recombine Surface "*";
   Transfinite Surface "*";
   Coherence;

   Mesh 1;
   Mesh 2;

   SetOrder 2;

   Mesh.Format = 1;
   Mesh.MshFileVersion = 2.2;
   Mesh.SaveAll = 0;
   Mesh.Binary = 0;

   Save "inclined_ellipse_2D.msh";
EndIf

// =================================================
// 3D EXTRUSION AND MESH SAVE
// =================================================
If (meshDim == 3)


      Extrude {0,0,Lz} {
         Surface{1:15};
         Layers{Nz};
         Recombine;}


   Color Black{Surface{:};}

   Color Black{Volume{:};}



   // TRANSFINITE

   Transfinite Surface "*";
   Recombine Surface "*";
   Transfinite Volume "*";


   Physical Volume("fluid")      = {1:15};
   Physical Surface("wall", 1)   =  {84, 106, 62, 124};
   Physical Surface("inlet", 2)  = {312, 176, 198};
   Physical Surface("outlet", 3) = {374, 352, 330};
   Physical Surface("top", 4)    =  {326, 260, 300, 290};
   Physical Surface("bottom", 5) =  {378, 246, 216, 194};

   Physical Surface("back", 6)  = {1, 4, 3, 2, 12, 11, 10, 5, 9, 8, 7, 6, 15, 14, 13};
   Physical Surface("front", 7)   = {379, 357, 335, 269, 159, 247, 225, 203, 181, 313, 291, 71, 93, 115, 137};


   Recombine Volume "*";
   Coherence;

   Mesh 1;
   Mesh 2;
   Mesh 3;

   SetOrder 2;

   Mesh.Format = 1;
   Mesh.MshFileVersion = 2.2;
   Mesh.SaveAll = 0;
   Mesh.Binary = 0;

   Save "inclined_ellipse_3D.msh";
EndIf
