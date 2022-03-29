// Gmsh project created on Mon Mar 14 17:08:21 2022
//+
a=1;
Point(1) = {0, 0, 0, a};
//+
Point(2) = {1, 0, 0, a};
//+
Point(3) = {2, 0, 0, a};
//+
Point(4) = {3, 0, 0, a};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Physical Point("Point0") = {1};
//+
//Physical Point("Point1") = {2};
//+
//Physical Point("Point2") = {3};
//+
Physical Point("Point3") = {4};
//+
Physical Curve("Global") = {1};
//+
Physical Curve("Gluing") = {2};
//+
Physical Curve("Local") = {3};
