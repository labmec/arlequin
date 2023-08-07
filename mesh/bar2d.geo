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

Transfinite Curve {1, 2, 3} = 2 Using Progression 1;
//+
Extrude {0, 1, 0} {
  Curve{1,2,3}; Layers{1}; Recombine;
}

//+
Physical Surface("Global") = {7};
//+
Physical Surface("Gluing") = {11};
//+
Physical Surface("Local") = {15};

Physical Curve("Left") = {5};
Physical Curve("Right") = {14};

Physical Curve("Neumann") = {1,2,3,4,8,12};
