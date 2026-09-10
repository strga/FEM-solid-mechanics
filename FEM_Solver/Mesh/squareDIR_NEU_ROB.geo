Mesh.MshFileVersion=2;

Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {1,1,0};
Point(4) = {0,1,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Line{1,2,3,4} = 17;
Transfinite Surface(1);

Physical Curve(1) = {3,4};
Physical Curve(10) = {1};
Physical Curve(20) = {2};

Physical Surface(5) = {1};