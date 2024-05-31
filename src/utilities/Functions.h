#include "../cell_types/Point.h"
#include <vector>
using namespace std;

#ifndef Functions_h
#define Functions_h

double advecting_circle(Point seed, double t, Point v = Point(0.5, 0.3), Point p0 = Point(0.5, 0.5), double r = 0.1);
double L1_error(vector<double> a, vector<double> b);
double advecting1D_stepfunc(Point seed, double t, double v, double a0, double b0);


#endif