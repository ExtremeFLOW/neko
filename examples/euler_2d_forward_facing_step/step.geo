// Forward-facing step 3D geometry for GMSH

// Step parameters
step_height = 0.2;
step_position = 0.6;
tunnel_height = 1.0;
tunnel_length = 3.0;
depth = 0.1;

// Mesh size (adjust for coarser/finer mesh)
mesh_size = 0.05;

// 2D Geometry (XY plane)
Point(1) = {0, 0, 0, mesh_size};          // Inlet bottom left
Point(2) = {0, tunnel_height, 0, mesh_size}; // Inlet top left
Point(3) = {step_position, tunnel_height, 0, mesh_size}; // Top wall at step
Point(4) = {step_position, 0, 0, mesh_size};  // Bottom at step base
Point(5) = {step_position, step_height, 0, mesh_size}; // Step top
Point(6) = {tunnel_length, step_height, 0, mesh_size};  // Outlet bottom
Point(7) = {tunnel_length, tunnel_height, 0, mesh_size}; // Outlet top

// Connect points with lines
Line(1) = {1, 4};  // Bottom before step
Line(2) = {4, 5};  // Step vertical face
Line(3) = {5, 6};  // Bottom after step
Line(4) = {6, 7};  // Outlet vertical
Line(5) = {7, 3};  // Top after step
Line(6) = {3, 2};  // Top before step
Line(7) = {2, 1};  // Inlet vertical

// Create surface
Line Loop(1) = {7, 1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// Extrude to 3D (thin in Z-direction)
extruded[] = Extrude {0, 0, depth} {
  Surface{1};
  Layers{1};
  Recombine;
};

// Boundary labels
Physical Surface("inlet") = {19};          // Extruded from line 7 (inlet)
Physical Surface("outlet") = {35};         // Extruded from line 4 (outlet)
Physical Surface("top_wall") = {39, 43}; // Extruded from lines 5-6 (top walls)
Physical Surface("bottom_wall") = {23, 27, 31}; // Extruded from lines 1-3 (bottom walls, including step)
Physical Surface("front") = {1};                   // Original surface (front)
Physical Surface("back") = {extruded[0]};          // Extruded surface (back)

// Physical groups
Physical Volume("fluid") = {extruded[1]};

// Mesh control for second-order elements
Mesh.ElementOrder = 2; // Use second-order elements (quad8/9, hex20/27)
Mesh.SecondOrderLinear = 0; // Use curved elements (set to 1 for linear elements)
Mesh.RecombinationAlgorithm = 1; // Blossom recombination
Mesh.RecombineAll = 1; // Recombine all surfaces
Mesh.Algorithm = 6; // Frontal-Delaunay for quads
Mesh.Optimize = 1; // Optimize mesh quality
Mesh.Smoothing = 5; // Smooth the mesh