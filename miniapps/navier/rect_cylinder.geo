// -----------------------------------------------------------------------------
//
//  Gmsh GEO tutorial 4
//
//  Built-in functions, holes in surfaces, annotations, entity colors
//
// -----------------------------------------------------------------------------

// As usual, we start by defining some variables:

Lc1 = 0.125;

Point(1) = {0, 0, 0, Lc1}; Point(2) = {1, 0, 0, Lc1};
Point(3) = {0, 2, 0, Lc1}; Point(4) = {1, 2, 0, Lc1};

//Point(1) = {0, 0, 0, Lc1}; Point(2) = {2, 0, 0, Lc1};
//Point(3) = {0, 1, 0, Lc1}; Point(4) = {2, 1, 0, Lc1};

Point(5) = {0.5, 0.5, 0, Lc1};
Point(6) = {0.45, 0.5, 0, Lc1};
Point(7) = {0.5, 0.55, 0, Lc1};
Point(8) = {0.55, 0.5, 0, Lc1};
Point(9) = {0.5, 0.45, 0, Lc1};

//Point(5) = {0.51, 0.5, 0, Lc1};
//Point(6) = {0.46, 0.5, 0, Lc1};
//Point(7) = {0.5, 0.55, 0, Lc1};
//Point(8) = {0.56, 0.5, 0, Lc1};
//Point(9) = {0.5, 0.45, 0, Lc1};


Line(1)  = {1 , 2};
Line(2)  = {2 , 4};
Line(3)  = {4 , 3};
Line(4)  = {3 , 1};

Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

// The third elementary entity is the surface. In order to define a simple
// rectangular surface from the four curves defined above, a curve loop has
// first to be defined. A curve loop is also identified by a tag (unique amongst
// curve loops) and defined by an ordered list of connected curves, a sign being
// associated with each curve (depending on the orientation of the curve to form
// a loop):

Curve Loop(1) = {4, 1, 2, 3};
Curve Loop(2) = {5, 6, 7, 8};

// We can then define the surface as a list of curve loops (only one here,
// representing the external contour, since there are no holes--see `t4.geo' for
// an example of a surface with a hole):

Plane Surface(1) = {1, 2};

Physical Curve(1) = {1};
Physical Curve(2) = {2};
Physical Curve(3) = {3};
Physical Curve(4) = {4};
Physical Curve(5) = {5, 6, 7, 8};
Physical Surface(1) = {1};

// Point(10) = {0.2, 0, 0, Lc1};
// Point(11) = {0.4, 0, 0, Lc1};
// Point(12) = {0.6, 0, 0, Lc1};
// Point(13) = {0.8, 0, 0, Lc1};

// Point{10} In Surface {1};
// Point{11} In Surface {1};
// Point{12} In Surface {1};
// Point{13} In Surface {1};

Mesh.Algorithm = 11; // Delaunay for quads
Mesh.AnisoMax = 0;
Mesh.ElementOrder = 3;
Mesh.MshFileVersion = 2.2;
// Mesh.LcIntegrationPrecision = 1.0e-15;
Mesh 2;
Save "rect-cylinder.msh";
//Exit;
