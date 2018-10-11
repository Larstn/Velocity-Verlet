Point(1) = {-0.2, 0, 0.0, 0.0105};
Point(2) = {-0.0095, 0, 0.0, 0.00063};
Point(3) = {0, 0, 0.0, 0.00063};
Point(4) = {0.0095, 0, 0.0, 0.00063};
Point(5) = {0.08, 0, 0.0, 0.0105};
Point(6) = {0.08, 0.045, 0.0, 0.0168};
Point(7) = {0, 0.045, 0.0, 0.0168};
Point(8) = {-0.2, 0.045, 0.0, 0.0168};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line Loop(100000) = {1, 2, 3, 4, 5, 6, 7, 8};
Point(9) = {-0.0095, 0.003, 0.0, 0.00063};
Point(10) = {-0.0075, 0.003, 0.0, 0.00063};
Point(11) = {-0.0075, 0.011, 0.0, 0.00126};
Point(12) = {-0.0095, 0.011, 0.0, 0.00126};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};
Line Loop(100001) = {9, 10, 11, 12};
Point(13) = {-0.0025, 0.003, 0.0, 0.00063};
Point(14) = {0.0025, 0.003, 0.0, 0.00063};
Point(15) = {0.0025, 0.011, 0.0, 0.00126};
Point(16) = {-0.0025, 0.011, 0.0, 0.00126};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};
Line Loop(100002) = {13, 14, 15, 16};
Point(17) = {0.0075, 0.003, 0.0, 0.00063};
Point(18) = {0.0095, 0.003, 0.0, 0.00063};
Point(19) = {0.0095, 0.011, 0.0, 0.00126};
Point(20) = {0.0075, 0.011, 0.0, 0.00126};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 17};
Line Loop(100003) = {17, 18, 19, 20};
Plane Surface(1) = {100000, 100001, 100002, 100003};
Physical Surface("Everwhere") = {1};
Physical Line(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Line(2) = {9, 10, 11, 12};
Physical Line(3) = {13, 14, 15, 16};
Physical Line(4) = {17, 18, 19, 20};