// =====================================================================
//  Flow over a wall-mounted hump  (low-Mach, heated downstream wall)
//
//  Large 2-D channel with a smooth cosine hump on the lower wall,
//  meshed structured (transfinite -> quads) and extruded in z to hexes.
//  Streamwise is split into 3 zones so elements PACK over the hump
//  (capture the curvature) and cluster just downstream (separation),
//  while the wall-normal spacing is fine at the lower wall.
//
//    gmsh -3 hump.geo -o hump.msh        # all-hex mesh
//    gmsh2nek ; rea2nbin hump.re2 hump.nmsh
// =====================================================================
SetFactory("Built-in");

// ---- domain (fairly large) -------------------------------------------
Lx = 25.0;  Ly = 3.0;  Lz = 4.0;
xc = 6.0;   wh = 2.0;   Hh = 1.0;       // hump centre / half-width / height
x1 = xc - wh;  x2 = xc + wh;            // hump foot points = 4, 8

// ---- element counts per streamwise zone ------------------------------
NxA = 16;   // upstream  [0,  x1]  (graded, dense toward the hump)
NxB = 48;   // hump      [x1, x2]  (FINE, captures the curvature)
NxC = 80;   // downstream[x2, Lx]  (graded, dense just after the hump)
Ny  = 32;   // wall-normal
Nz  = 30;   // spanwise
rGrow = 1.09;   // wall-normal growth -> fine cells at the lower wall
pA    = 0.90;   // upstream  streamwise progression (<1: dense toward x1)
pC    = 1.028;  // downstream streamwise progression (>1: dense after x2)

// ---- corner / foot points --------------------------------------------
pA0=newp; Point(pA0)={0 ,0 ,0,1};  pB0=newp; Point(pB0)={x1,0 ,0,1};
pC0=newp; Point(pC0)={x2,0 ,0,1};  pD0=newp; Point(pD0)={Lx,0 ,0,1};
pAt=newp; Point(pAt)={0 ,Ly,0,1};  pBt=newp; Point(pBt)={x1,Ly,0,1};
pCt=newp; Point(pCt)={x2,Ly,0,1};  pDt=newp; Point(pDt)={Lx,Ly,0,1};

// ---- hump spline interior points -------------------------------------
nb=120; ph[]={};
For i In {1:nb-1}
  xx = x1 + (x2-x1)*i/nb;
  yy = 0.5*Hh*(1.0 + Cos(Pi*(xx-xc)/wh));
  pp=newp; Point(pp)={xx,yy,0,1}; ph[]+=pp;
EndFor

// ---- edges -----------------------------------------------------------
cbA=newc; Line(cbA)  ={pA0,pB0};                 // bottom upstream
cbB=newc; Spline(cbB)={pB0, ph[], pC0};          // bottom hump
cbC=newc; Line(cbC)  ={pC0,pD0};                 // bottom downstream
ctA=newc; Line(ctA)  ={pAt,pBt};                 // top upstream
ctB=newc; Line(ctB)  ={pBt,pCt};                 // top hump
ctC=newc; Line(ctC)  ={pCt,pDt};                 // top downstream
vL =newc; Line(vL)   ={pA0,pAt};                 // inlet  (x=0)
vR =newc; Line(vR)   ={pD0,pDt};                 // outlet (x=Lx)

ll=newll; Curve Loop(ll)={cbA,cbB,cbC, vR, -ctC,-ctB,-ctA, -vL};
ss=news;  Plane Surface(ss)={ll};

// ---- structured mesh: matching bottom/top streamwise distributions ---
Transfinite Curve {cbA, ctA} = NxA+1 Using Progression pA;   // dense toward hump
Transfinite Curve {cbB, ctB} = NxB+1;                        // fine, uniform
Transfinite Curve {cbC, ctC} = NxC+1 Using Progression pC;   // dense after hump
Transfinite Curve {vL,  vR } = Ny+1  Using Progression rGrow;// fine at lower wall
Transfinite Surface {ss} = {pA0, pD0, pDt, pAt};
Recombine Surface {ss};

// ---- extrude in z -> hexes -------------------------------------------
// out[0]=back(z=Lz), out[1]=volume, out[2..]=laterals in loop-curve order
out[] = Extrude {0,0,Lz} { Surface{ss}; Layers{Nz}; Recombine; };

// loop order {cbA,cbB,cbC, vR, -ctC,-ctB,-ctA, -vL} ->
//  out[2]=cbA out[3]=cbB out[4]=cbC out[5]=vR out[6]=ctC out[7]=ctB out[8]=ctA out[9]=vL
Physical Volume ("fluid")      = {out[1]};
Physical Surface("lower_wall") = {out[2], out[3], out[4]};   // heated part via x in driver
Physical Surface("outlet")     = {out[5]};
Physical Surface("upper_wall") = {out[6], out[7], out[8]};
Physical Surface("inlet")      = {out[9]};
Physical Surface("front")      = {ss};        // z=0  (periodic pair)
Physical Surface("back")       = {out[0]};    // z=Lz (periodic pair)

Mesh.RecombineAll = 1;
Mesh.Recombine3DAll = 1;
