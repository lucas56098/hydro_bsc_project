#include "../cell_types/Point.h"
#include <vector>
using namespace std;

#ifndef Functions_h
#define Functions_h

// advection
double advecting_circle(Point seed, double t, Point v = Point(0.5, 0.3), Point p0 = Point(0.5, 0.5), double r = 0.1);
double advecting1D_stepfunc(Point seed, double t, double v, double a0, double b0);

// shallow water
double swe1D_dam_break(Point seed, double t, double hl = 2, double hr = 1, double hm = 1.45384089237457292398403296829201281070709228515625, double g = 1, double x0 = 0.5);

// L1 error
double L1_error(vector<double> a, vector<double> b);

// other
string format_time(double seconds);

#endif