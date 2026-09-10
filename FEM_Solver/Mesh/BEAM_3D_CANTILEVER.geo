Mesh.MshFileVersion = 2;

// ------------------------- PARAMETERS ------------------------- //
H = 0.25;    // Height of beam (y-axis)
L  = 20 * H;   // Length of beam (x-axis)
W = H/2;    // Width / depth of beam (z-axis) from -w to w is the length
cl = 0.05;   // Mesh characteristic length

// ------------------------- POINTS ------------------------- //
Point(1) = {0, 0, -W, cl}; // Bottom-left-front
Point(2) = {L, 0, -W, cl}; // Bottom-right-front
Point(3) = {L, H, -W, cl}; // Bottom-right-back
Point(4) = {0, H, -W, cl}; // Bottom-left-back
Point(5) = {0, 0, W, cl}; // Top-left-front
Point(6) = {L, 0, W, cl}; // Top-right-front
Point(7) = {L, H, W, cl}; // Top-right-back
Point(8) = {0, H, W, cl}; // Top-left-back

// ------------------------- LINES ------------------------- //
Line(1)  = {1, 2}; Line(2)  = {2, 3}; Line(3)  = {3, 4}; Line(4)  = {4, 1};
Line(5)  = {5, 6}; Line(6)  = {6, 7}; Line(7)  = {7, 8}; Line(8)  = {8, 5};
Line(9)  = {1, 5}; Line(10) = {2, 6}; Line(11) = {3, 7}; Line(12) = {4, 8};

// ------------------------- SURFACES ------------------------- //
// Define 6 faces for the cube
Line Loop(13) = {-1, -2, -3, -4};   Plane Surface(14) = {13};  // Bottom
Line Loop(15) = {5, 6, 7, 8};   Plane Surface(16) = {15};  // Top
Line Loop(17) = {1, 10, -5, -9}; Plane Surface(18) = {17}; // Front
Line Loop(19) = {2, 11, -6, -10}; Plane Surface(20) = {19}; // Right
Line Loop(21) = {3, 12, -7, -11}; Plane Surface(22) = {21}; // Back
Line Loop(23) = {4, 9, -8, -12}; Plane Surface(24) = {23}; // Left

// ------------------------- VOLUME ------------------------- //
Surface Loop(25) = {14, 16, 18, 20, 22, 24};
Volume(26) = {25};

// -------------------- TRANSFINITE MESH -------------------- //
// Transfinite lines (set divisions = 2 for simplicity)
/*
Transfinite Line {1:12} = 6;

// Transfinite surfaces (defined by corner points)
Transfinite Surface {14} = {1,2,3,4};
Transfinite Surface {16} = {5,6,7,8};
Transfinite Surface {18} = {1,2,6,5};
Transfinite Surface {20} = {2,3,7,6};
Transfinite Surface {22} = {3,4,8,7};
Transfinite Surface {24} = {4,1,5,8};

// Transfinite volume (defined by 8 corner points)
Transfinite Volume {26} = {1,2,3,4,5,6,7,8};
*/
// -------------------- PHYSICAL GROUPS -------------------- //
Physical Volume(100) = {26}; 	// Cube // 100 - STEEL, 130 - SILICONE

//-> DIRICHLET + NEUMANN:
Physical Surface(1) = {24};			//DIR_CANTILEVER
Physical Surface(10) = {20}; 			//NEUMANN_LOAD
Physical Surface(11) = {14, 16, 18, 22}; 	//NEUMANN_FREE

// ------------------------- MESH ALGORITHM ------------------------- //
//Mesh.Algorithm = 1; // Structured mesh generation