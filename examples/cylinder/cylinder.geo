/*
   -----------  gmsh script for generating 3D cylinder  --------------

   -------------------------------------------------------------------
   By Daniele Massaro, dmassaro@kth.se
      Linne FLOW Centre, KTH Mechanics, Stockholm, Sweden
   Updated ever so slightly by Martin Karp, makarp@kth.se, 2024
*/




meshDim=3;  

//  ----Nxll---|----Nxc---|----Nxlr---|-----------------Nxw--------------|
//  |          |          |           |                                  |
//  Nytp      Nytp       Nytp        Nytp                                |
//  |          |          |           |                                  |
//  ----Nxll---|----Nxc---|----Nxlr---|-----------------Nxw--------------|
//  |          |          |           |                                  |
//  Nyc       Nyc  Nb,Nc Nyc         Nyc                                Nyc
//  |          |          |           |                                  |
//  ----Nxll---|----Nxc---|----Nxlr---|-----------------Nxw--------------|
//  |          |          |           |
//  Nytb      Nytb       Nytb        Nytb                                |
//  |          |          |           |                                  |
//  ----Nxll---|----Nxc---|----Nxlr---|-----------------Nxw--------------|


// NOTE: Number of elemenets = number of points-1
// 
Nxc  = 5 ;      Rxc  = 1.00;   // Progression ratio towards right (elements clustered on right side of Nxc box if Rxc<0)
Nxw  = 10;      Rxw  = 1.05;   // Progression ratio towards right  (elements clustered on right side of Nxw box if Rxw<0)
Nxlr =  5;      Rxlr = 1.0;    // Progression ratio towards right (elements clustered on right side of Nxlr box if Rxlr<0)
Nxll =  4;      Rxll = 1.35;    // Progression ratio towards left (elements clustered on left side of Nxll box if Rxlr<0)

Nyc  = Nxc;     Ryc  = 1.0;    // Progression ratio towards top (elements clustered on top side of Nyc box if Ryc<0)
Nytb =  5;       Rytb = 1.0;   // Progression ratio towards bottom (elements clustered on bottom side of Nytb box if Rytb<0)
Nytp =  5;       Rytp = 1.0;   // Progression ratio towards bottom (elements clustered on bottom side of Nytp box if Rytp<0)

Nb   =  6;      Rb   = 0.9;    // Progression ratio towards center (elements clustered on the central part side of Nb box if Rb<0)
Nc  = Nxc;      Rc   = 1.00;   // Progression on the cirle (keep 1)




// ---- Coordinates ----
// Internal box (core)
xc = 4;
yc = 4;
R  = 0.5; D = 1;
// Cylinder xr=yr= R*sqrt(2)/2
xr = 0.35355339; yr = 0.35355339;  
// External box
xe = 10;
ye = 15;
// Wake (yw=ye)
xw = 30;
// -----------------------------

// ---- 3D cylinder ----
If (meshDim==3)
   Lz = 6;
   Nz = 6;
EndIf
// -----------------------------


// ---- Points ----
// INTERNAL BOX
Point(1) = {-xc, -yc, 0, 1.0};
Point(2) = {xc, -yc, 0, 1.0};
Point(3) = {-xc, yc, 0, 1.0};
Point(4) = {xc, yc, 0, 1.0};

// CYLINDER POINTS
Point(5) = {-xr, -yr, 0, 1.0};
Point(6) = {xr, -yr, 0, 1.0};
Point(7) = {xr, yr, 0, 1.0};
Point(8) = {-xr, yr, 0, 1.0};
Point(9) = {0, 0, 0, 1.0};

// EXTERNAL BOX
Point(10) = {-xe, -ye, 0, 1.0};
Point(11) = {xe, -ye, 0, 1.0};
Point(12) = {-xe, ye, 0, 1.0};
Point(13) = {xe, ye, 0, 1.0};

Point(14) = {-xe, -yc, 0, 1.0};
Point(15) = {xe, -yc, 0, 1.0};
Point(16) = {-xe, yc, 0, 1.0};
Point(17) = {xe, yc, 0, 1.0};

Point(22) = {xc, -ye, 0, 1.0};
Point(23) = {-xc, -ye, 0, 1.0};
Point(24) = {xc, ye, 0, 1.0};
Point(25) = {-xc, ye, 0, 1.0};

// WAKE BOX
Point(18) = {xw, ye, 0, 1.0};
Point(19) = {xw, -ye, 0, 1.0};

Point(20) = {xw, yc, 0, 1.0};
Point(21) = {xw, -yc, 0, 1.0};
// -----------------------------




// ---- Lines ----
// INTERNAL BOX LINES
Line(1) = {1, 2}; Transfinite Line {1} = Nxc Using Progression Rxc;
Line(2) = {3, 4}; Transfinite Line {2} = Nxc Using Progression Rxc;
Line(3) = {1, 3}; Transfinite Line {3} = Nyc Using Progression Ryc;
Line(4) = {2, 4}; Transfinite Line {4} = Nyc Using Progression Ryc;

// BIG CYLINDER LINES
Circle(5) = {5, 9, 6}; Transfinite Line {5} = Nc Using Progression Rc;
Circle(6) = {6, 9, 7}; Transfinite Line {6} = Nc Using Progression Rc;
Circle(7) = {7, 9, 8}; Transfinite Line {7} = Nc Using Progression Rc;
Circle(8) = {8, 9, 5}; Transfinite Line {8} = Nc Using Progression Rc;

// DIAGONAL LINES
Line(9) = {1, 5}; Transfinite Line {9} = Nb Using Progression Rb;
Line(10) = {2, 6}; Transfinite Line {10} = Nb Using Progression Rb;
Line(11) = {4, 7}; Transfinite Line {11} = Nb Using Progression Rb;
Line(12) = {3, 8}; Transfinite Line {12} = Nb Using Progression Rb;


// RIGHT BOX
Line(17) = {4, 17}; Transfinite Line {17} = Nxlr Using Progression Rxlr;
Line(18) = {2, 15}; Transfinite Line {18} = Nxlr Using Progression Rxlr;
Line(19) = {15, 17}; Transfinite Line {19} = Nyc Using Progression Ryc;

// LEFT BOX
Line(20) = {1, 14};  Transfinite Line {20} = Nxll Using Progression Rxll;
Line(21) = {3, 16};  Transfinite Line {21} = Nxll Using Progression Rxll;
Line(22) = {14, 16}; Transfinite Line {22} = Nyc Using Progression Ryc;

// BOTTOM BOX
Line(23) = {15, 11}; Transfinite Line {23} = Nytb Using Progression Rytb;
Line(24) = {15, 14}; Transfinite Line {24} = Nxw Using Progression Rxw;
Line(25) = {14, 10}; Transfinite Line {25} = Nytb Using Progression Rytb;

Line(40) = {1, 23};  Transfinite Line {40} = Nytb Using Progression Rytb;
Line(41) = {2, 22}; Transfinite Line {41} = Nytb Using Progression Rytb;
Line(42) = {23, 22}; Transfinite Line {42} = Nxc Using Progression Rxc;
Line(43) = {23, 10}; Transfinite Line {43} = Nxll Using Progression Rxll;
Line(44) = {22, 11}; Transfinite Line {44} = Nxlr Using Progression Rxlr;

// TOP BOX
Line(26) = {16, 17}; Transfinite Line {26} = Nxw Using Progression Rxw;
Line(27) = {13, 17}; Transfinite Line {27} = Nytb Using Progression Rytp;
Line(28) = {12, 16}; Transfinite Line {28} = Nytb Using Progression Rytp;

Line(45) = {25, 3};  Transfinite Line {45} = Nytb Using Progression Rytp;
Line(46) = {25, 24}; Transfinite Line {46} = Nxc Using Progression Rxc;
Line(47) = {24, 4}; Transfinite Line {47} = Nytb Using Progression Rytp;
Line(48) = {25, 12}; Transfinite Line {48} = Nxll Using Progression Rxll;
Line(49) = {24, 13}; Transfinite Line {49} = Nxlr Using Progression Rxlr;

// WAKE BOX
Line(32) = {11, 19}; Transfinite Line {32} = Nxw Using Progression Rxw;

Line(34) = {13, 18}; Transfinite Line {34} = Nxw Using Progression Rxw;

Line(35) = {17, 20}; Transfinite Line {35} = Nxw Using Progression Rxw;
Line(36) = {15, 21}; Transfinite Line {36} = Nxw Using Progression Rxw;
Line(37) = {21, 20}; Transfinite Line {37} = Nyc Using Progression Ryc;
Line(38) = {18, 20}; Transfinite Line {38} = Nytb Using Progression Rytp;
Line(39) = {21, 19}; Transfinite Line {39} = Nytb Using Progression Rytb;
// -----------------------------



// ---- Surfaces ----
//+
Curve Loop(1) = {10, 6, -11, -4};
Plane Surface(1) = {1};
//+
Curve Loop(2) = {11, 7, -12, 2};
Plane Surface(2) = {2};
//+
Curve Loop(3) = {12, 8, -9, 3};
Plane Surface(3) = {3};
//+
Curve Loop(4) = {5, -10, -1, 9};
Plane Surface(4) = {4};
//+
Curve Loop(5) = {4, 17, -19, -18};
Plane Surface(5) = {5};
//+
Curve Loop(6) = {-3, 20, 22, -21};
Plane Surface(6) = {6};
//+
Curve Loop(7) = {40, 43, -25, -20};
Plane Surface(7) = {7};
//+
Curve Loop(8) = {41, -42, -40, 1};
Plane Surface(8) = {8};
//+
Curve Loop(9) = {-23, -18, 41, 44};
Plane Surface(9) = {9};
//+
Curve Loop(10) = {-27, -49, 47, 17};
Plane Surface(10) = {10};
//+
Curve Loop(11) = {47, -2, -45, 46};
Plane Surface(11) = {11};
//+
Curve Loop(12) = {-48, 45, 21, -28};
Plane Surface(12) = {12};
//+
Curve Loop(13) = {-27, 34, 38, -35};
Plane Surface(13) = {13};
//+
Curve Loop(14) = {19, 35, -37, -36};
Plane Surface(14) = {14};
//+
Curve Loop(15) = {-23, 36, 39, -32};
Plane Surface(15) = {15};
Color Black{Surface{1:15};}

// TRANSFINITE
Transfinite Surface "*";
// RECOMBINE
Recombine Surface "*";


If (meshDim==2)

   //+ Physical surface [need to be defined in some CFD codes, as Nek5000]
   Physical Curve("wall") = {7, 6, 5, 8};
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

   Save "ext_cyl.msh";

EndIf


If (meshDim==3)


   Extrude {0, 0, Lz} {
     Surface{1:15}; 
     Layers{Nz};
     Recombine; 
   }
   Color Black{Surface{:};}
   Color Black{Volume{:};}

   // TRANSFINITE
   Transfinite Surface "*";
   Recombine Surface "*";
   // TRANSFINITE
   Transfinite Volume "*";


   //+ Physical surface [need to be defined in some CFD codes, as Nek5000]
   Physical Surface("inlet") = {312, 176, 198};
   Physical Surface("outlet") = {374, 352, 330};
   Physical Surface("top") = {379, 357, 335, 269, 159, 247, 225, 203, 181, 313, 291, 71, 93, 115, 137};
   Physical Surface("bottom") = {1, 4, 3, 2, 12, 11, 10, 5, 9, 8, 7, 6, 15, 14, 13};
   Physical Surface("front") = {378, 246, 216, 194};
   Physical Surface("back") = {326, 260, 300, 290};
   Physical Surface("wall") = {84, 106, 62, 124};
   Physical Volume("Fluid") =  {1:15};




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

   Save "3D_ext_cyl.msh";




EndIf











