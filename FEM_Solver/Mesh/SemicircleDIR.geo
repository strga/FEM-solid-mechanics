Mesh.MshFileVersion=2;

Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {2,0,0};
Point(4) = {3,0,0};
Point(5) = {4,0,0};
Point(6) = {4,2,0};
Point(7) = {0,2,0};
Point(8) = {0,2,0};

Line(1) = {1,2};
Line(2) = {4,5};
Line(3) = {5,6};
Line(4) = {6,7};
Line(5) = {7,1};

Circle(6) = {4,3,2};

Curve Loop(1) = {4, 5, 1, -6, 2, 3};

Plane Surface(1) = {1};

Physical Curve(1) = {5, 1, 6, 2, 3, 4};

Physical Surface(5) = {1};