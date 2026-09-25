meshDim=3;  

// ====================  PARAMETERS  =======================

// Grid
N_inlet = 12;   R_inlet = 1.1;   
Nx_cyl = 10;   Rx_cyl = 1.0;   
Nx_bridge = 10;   Rx_bridge = 1.0;   
N_wake = 25;   R_wake = 1.05;  

Ny_core   = 10;   Ry_core   = 1.0;   // Middle Row
Ny_outer  = 12;   Ry_outer  = 1.2;   // Top/Bottom Rows

Nb = 8;    Rb = 1.0;   // Butterfly Radial
Nc = 10;   Rc = 1.0;   // Butterfly Circumferential

// Dimensions 
R = 0.5;                
xr = R*Sqrt(2)/2;        
yc = 1.5;                
xc = 1.5;               

Sep = 7.0;                
x_c1 = -Sep/2;            
x_c2 =  Sep/2;            

// Boundaries
x_inlet = -10.0;
x_outlet =  15.0;
y_top =  6.0;
y_bot = -6.0;

If (meshDim==3)
   Lz = 0.5;
   Nz = 1; 
EndIf

// ======================  POINTS  =========================
x0 = x_inlet;
x1 = x_c1 - xc; x2 = x_c1 + xc; 
x3 = x_c2 - xc; x4 = x_c2 + xc; 
x5 = x_outlet;
y0 = y_bot; y1 = -yc; y2 = yc; y3 = y_top;

Point(1) = {x0, y0, 0, 1.0}; Point(2) = {x0, y1, 0, 1.0}; Point(3) = {x0, y2, 0, 1.0}; Point(4) = {x0, y3, 0, 1.0};
Point(5) = {x1, y0, 0, 1.0}; Point(6) = {x1, y1, 0, 1.0}; Point(7) = {x1, y2, 0, 1.0}; Point(8) = {x1, y3, 0, 1.0};
Point(9) = {x2, y0, 0, 1.0}; Point(10)= {x2, y1, 0, 1.0}; Point(11)= {x2, y2, 0, 1.0}; Point(12)= {x2, y3, 0, 1.0};
Point(13)= {x3, y0, 0, 1.0}; Point(14)= {x3, y1, 0, 1.0}; Point(15)= {x3, y2, 0, 1.0}; Point(16)= {x3, y3, 0, 1.0};
Point(17)= {x4, y0, 0, 1.0}; Point(18)= {x4, y1, 0, 1.0}; Point(19)= {x4, y2, 0, 1.0}; Point(20)= {x4, y3, 0, 1.0};
Point(21)= {x5, y0, 0, 1.0}; Point(22)= {x5, y1, 0, 1.0}; Point(23)= {x5, y2, 0, 1.0}; Point(24)= {x5, y3, 0, 1.0};
Point(100) = {x_c1, 0, 0, 1.0}; 
Point(101) = {x_c1-xr, -xr, 0, 1.0}; Point(102) = {x_c1+xr, -xr, 0, 1.0};
Point(103) = {x_c1+xr,  xr, 0, 1.0}; Point(104) = {x_c1-xr,  xr, 0, 1.0};
Point(200) = {x_c2, 0, 0, 1.0}; 
Point(201) = {x_c2-xr, -xr, 0, 1.0}; Point(202) = {x_c2+xr, -xr, 0, 1.0};
Point(203) = {x_c2+xr,  xr, 0, 1.0}; Point(204) = {x_c2-xr,  xr, 0, 1.0};

// ======================  LINES  ==========================
// Horizontal Lines (Bottom)
Line(1) = {1, 5}; Transfinite Line {1} = N_inlet Using Progression R_inlet;
Line(5) = {5, 9}; Transfinite Line {5} = Nx_cyl Using Progression Rx_cyl;
Line(9) = {9, 13}; Transfinite Line {9} = Nx_bridge Using Progression Rx_bridge;
Line(13)= {13,17}; Transfinite Line {13}= Nx_cyl Using Progression Rx_cyl;
Line(17)= {17,21}; Transfinite Line {17}= N_wake Using Progression R_wake;
// Horizontal Lines (Core Bottom) 
Line(2) = {2, 6}; Transfinite Line {2} = N_inlet Using Progression R_inlet;
Line(6) = {6, 10}; Transfinite Line {6} = Nx_cyl Using Progression Rx_cyl;
Line(10)= {10,14}; Transfinite Line {10}= Nx_bridge Using Progression Rx_bridge;
Line(14)= {14,18}; Transfinite Line {14}= Nx_cyl Using Progression Rx_cyl;
Line(18)= {18,22}; Transfinite Line {18}= N_wake Using Progression R_wake;
// Horizontal Lines (Core Top) 
Line(3) = {3, 7}; Transfinite Line {3} = N_inlet Using Progression R_inlet;
Line(7) = {7, 11}; Transfinite Line {7} = Nx_cyl Using Progression Rx_cyl;
Line(11)= {11,15}; Transfinite Line {11}= Nx_bridge Using Progression Rx_bridge;
Line(15)= {15,19}; Transfinite Line {15}= Nx_cyl Using Progression Rx_cyl;
Line(19)= {19,23}; Transfinite Line {19}= N_wake Using Progression R_wake;
// Horizontal Lines (Top) 
Line(4) = {4, 8}; Transfinite Line {4} = N_inlet Using Progression R_inlet;
Line(8) = {8, 12}; Transfinite Line {8} = Nx_cyl Using Progression Rx_cyl;
Line(12)= {12,16}; Transfinite Line {12}= Nx_bridge Using Progression Rx_bridge;
Line(16)= {16,20}; Transfinite Line {16}= Nx_cyl Using Progression Rx_cyl;
Line(20)= {20,24}; Transfinite Line {20}= N_wake Using Progression R_wake;
// Vertical Lines (Lower) 
Line(21) = {1, 2}; Transfinite Line {21} = Ny_outer Using Progression 1/Ry_outer;
Line(24) = {5, 6}; Transfinite Line {24} = Ny_outer Using Progression 1/Ry_outer;
Line(27) = {9, 10}; Transfinite Line {27} = Ny_outer Using Progression 1/Ry_outer;
Line(30) = {13,14}; Transfinite Line {30} = Ny_outer Using Progression 1/Ry_outer;
Line(33) = {17,18}; Transfinite Line {33} = Ny_outer Using Progression 1/Ry_outer;
Line(36) = {21,22}; Transfinite Line {36} = Ny_outer Using Progression 1/Ry_outer;
// Vertical Lines (Core) 
Line(23) = {2, 3}; Transfinite Line {23} = Ny_core Using Progression Ry_core;
Line(26) = {6, 7}; Transfinite Line {26} = Ny_core Using Progression Ry_core;
Line(29) = {10,11}; Transfinite Line {29} = Ny_core Using Progression Ry_core;
Line(32) = {14,15}; Transfinite Line {32} = Ny_core Using Progression Ry_core;
Line(35) = {18,19}; Transfinite Line {35} = Ny_core Using Progression Ry_core;
Line(38) = {22,23}; Transfinite Line {38} = Ny_core Using Progression Ry_core;
// Vertical Lines (Upper)
Line(22) = {3, 4}; Transfinite Line {22} = Ny_outer Using Progression Ry_outer;
Line(25) = {7, 8}; Transfinite Line {25} = Ny_outer Using Progression Ry_outer;
Line(28) = {11,12}; Transfinite Line {28} = Ny_outer Using Progression Ry_outer;
Line(31) = {15,16}; Transfinite Line {31} = Ny_outer Using Progression Ry_outer;
Line(34) = {19,20}; Transfinite Line {34} = Ny_outer Using Progression Ry_outer;
Line(37) = {23,24}; Transfinite Line {37} = Ny_outer Using Progression Ry_outer;

// Cylinders
// Cyl 1
Circle(40) = {101, 100, 102}; Transfinite Line {40} = Nc Using Progression Rc;
Circle(41) = {102, 100, 103}; Transfinite Line {41} = Nc Using Progression Rc;
Circle(42) = {103, 100, 104}; Transfinite Line {42} = Nc Using Progression Rc;
Circle(43) = {104, 100, 101}; Transfinite Line {43} = Nc Using Progression Rc;

Line(50) = {6, 101};  Transfinite Line {50} = Nb Using Progression Rb;
Line(51) = {10, 102}; Transfinite Line {51} = Nb Using Progression Rb;
Line(52) = {11, 103}; Transfinite Line {52} = Nb Using Progression Rb;
Line(53) = {7, 104};  Transfinite Line {53} = Nb Using Progression Rb;

// Cyl 2
Circle(44) = {201, 200, 202}; Transfinite Line {44} = Nc Using Progression Rc;
Circle(45) = {202, 200, 203}; Transfinite Line {45} = Nc Using Progression Rc;
Circle(46) = {203, 200, 204}; Transfinite Line {46} = Nc Using Progression Rc;
Circle(47) = {204, 200, 201}; Transfinite Line {47} = Nc Using Progression Rc;

Line(54) = {14, 201}; Transfinite Line {54} = Nb Using Progression Rb;
Line(55) = {18, 202}; Transfinite Line {55} = Nb Using Progression Rb;
Line(56) = {19, 203}; Transfinite Line {56} = Nb Using Progression Rb;
Line(57) = {15, 204}; Transfinite Line {57} = Nb Using Progression Rb;

// ====================  SURFACES  =========================
// BOTTOM ROW 
// Inlet Bot
Curve Loop(101) = {1, 24, -2, -21}; Plane Surface(101)={101}; Transfinite Surface{101}; Recombine Surface{101};
// Cyl 1 Bot
Curve Loop(201) = {5, 27, -6, -24}; Plane Surface(201)={201}; Transfinite Surface{201}; Recombine Surface{201};
// Bridge Bot
Curve Loop(301) = {9, 30, -10, -27}; Plane Surface(301)={301}; Transfinite Surface{301}; Recombine Surface{301};
// Cyl 2 Bot
Curve Loop(401) = {13, 33, -14, -30}; Plane Surface(401)={401}; Transfinite Surface{401}; Recombine Surface{401};
// Outlet Bot
Curve Loop(501) = {17, 36, -18, -33}; Plane Surface(501)={501}; Transfinite Surface{501}; Recombine Surface{501};

// TOP ROW 
// Inlet Top
Curve Loop(103) = {3, 25, -4, -22}; Plane Surface(103)={103}; Transfinite Surface{103}; Recombine Surface{103};
// Cyl 1 Top
Curve Loop(203) = {7, 28, -8, -25}; Plane Surface(203)={203}; Transfinite Surface{203}; Recombine Surface{203};
// Bridge Top
Curve Loop(303) = {11, 31, -12, -28}; Plane Surface(303)={303}; Transfinite Surface{303}; Recombine Surface{303};
// Cyl 2 Top
Curve Loop(403) = {15, 34, -16, -31}; Plane Surface(403)={403}; Transfinite Surface{403}; Recombine Surface{403};
// Outlet Top
Curve Loop(503) = {19, 37, -20, -34}; Plane Surface(503)={503}; Transfinite Surface{503}; Recombine Surface{503};

// Inlet Mid
Curve Loop(102) = {2, 26, -3, -23}; Plane Surface(102)={102}; Transfinite Surface{102}; Recombine Surface{102};
// Bridge Mid
Curve Loop(302) = {10, 32, -11, -29}; Plane Surface(302)={302}; Transfinite Surface{302}; Recombine Surface{302};
// Outlet Mid
Curve Loop(502) = {18, 38, -19, -35}; Plane Surface(502)={502}; Transfinite Surface{502}; Recombine Surface{502};

// Cyl 1
Curve Loop(210) = {6, 51, -40, -50}; Plane Surface(210)={210}; // Bot
Curve Loop(211) = {29, 52, -41, -51}; Plane Surface(211)={211}; // Right
Curve Loop(212) = {-7, 53, -42, -52}; Plane Surface(212)={212}; // Top
Curve Loop(213) = {-26, 50, -43, -53}; Plane Surface(213)={213}; // Left
Transfinite Surface{210, 211, 212, 213}; Recombine Surface{210, 211, 212, 213};

// Cyl 2
Curve Loop(410) = {14, 55, -44, -54}; Plane Surface(410)={410}; // Bot
Curve Loop(411) = {35, 56, -45, -55}; Plane Surface(411)={411}; // Right
Curve Loop(412) = {-15, 57, -46, -56}; Plane Surface(412)={412}; // Top
Curve Loop(413) = {-32, 54, -47, -57}; Plane Surface(413)={413}; // Left
Transfinite Surface{410, 411, 412, 413}; Recombine Surface{410, 411, 412, 413};

// ====================  EXTRUSION  ========================
If (meshDim==3)
   Extrude {0, 0, Lz} {
      Surface{101:103, 201, 203, 210:213, 301:303, 401, 403, 410:413, 501:503};
      Layers{Nz};
      Recombine;
   }
   
   Physical Volume("fluid") = {1:21}; 
   
   Physical Surface("inlet") = Surface In BoundingBox{x_inlet-0.1, -100, -100, x_inlet+0.1, 100, 100};
   Physical Surface("outlet") = Surface In BoundingBox{x_outlet-0.1, -100, -100, x_outlet+0.1, 100, 100};
   Physical Surface("top") = Surface In BoundingBox{-100, y_top-0.1, -100, 100, y_top+0.1, 100};
   Physical Surface("bottom") = Surface In BoundingBox{-100, y_bot-0.1, -100, 100, y_bot+0.1, 100};
   Physical Surface("cyl1_wall") = {696, 674, 652, 630}; 
   Physical Surface("cyl2_wall") = {872, 850, 828, 894};
   Physical Surface("front") = {101:103, 201, 203, 210:213, 301:303, 401, 403, 410:413, 501:503};
   Physical Surface("back") = {525, 547, 569, 613, 679, 657, 701, 635, 591, 723, 745, 767, 811, 877, 899, 855, 833, 789, 965, 943, 921};

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
