// Gmsh project created on Fri Aug  4 10:27:14 2023
SetFactory("OpenCASCADE");
a=.5;
Point(1) = {0, 0, 0, a};
//+
Point(2) = {2, 0, 0, a};
//+
Point(3) = {2, 1, 0, a};
//+
Point(4) = {0, 1, 0, a};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
Circle(9) = {1, .5, 0., 0.05, 0, 2*Pi};
//+
Circle(10) = {1, .5, 0., 0.1, 0, 2*Pi};
//+
Circle(11) = {1, .5, 0., 0.12, 0, 2*Pi};

//+
Curve Loop(2) = {4, 1, 2, 3};
//+
Curve Loop(3) = {11};
//+
Plane Surface(2) = {2, 3}; //Recombine Surface{2};
//+
Curve Loop(4) = {11};
//+
Curve Loop(5) = {10};
//+
Plane Surface(3) = {4,5};// Recombine Surface{3};
//+
Curve Loop(6) = {10};
//+
Curve Loop(7) = {9};
//+
Plane Surface(4) = {6, 7}; //Recombine Surface{4};
//+
Physical Curve("Left") = {4};
//+
Physical Curve("Right") = {2};
//Physical Curve("Circle") = {9};

Physical Curve("Neumann") = {1,3};
//+
Physical Surface("Global") = {2};
Physical Surface("Gluing") = {3};
Physical Surface("Local") = {4};
