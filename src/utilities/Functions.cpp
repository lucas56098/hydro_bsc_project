#include "Functions.h"
#include <cmath>
#include <vector>
#include <iostream>

// Returns analytic solution of advecting circle for a given seedpoint
double advecting_circle(Point seed, double t, Point v, Point p0, double r) {

    Point p = Point(p0.x + t * v.x, p0.y + t * v.y);
    double dist = sqrt((p.x - seed.x)*(p.x - seed.x) + (p.y - seed.y)*(p.y - seed.y));

    if (dist <= r) {
        return 1;
    }
    return 0;

}
// Returns analytic solution of advecting stepfunc in 1D for a given seedpoint
double advecting1D_stepfunc(Point seed, double t, double v, double a0, double b0) {

    if (seed.x >= a0 + v*t && seed.x <= b0 + v*t) {
        return 1;
    }
    return 0;

}

// Calculates L1_error given two vectors
double L1_error(vector<double> a, vector<double> b) {

    // make sure that a and b have the same size
    if (a.size() != b.size()) {
        cerr << "Unable to calculate L1_error since a.size() != b.size()" << endl;
        exit(EXIT_FAILURE);
    }

    // sum up total error
    double total_error = 0;
    for (int i = 0; i < a.size(); i++) {
        total_error += abs(a[i] - b[i]);
    }

    // divide by N
    double L1 = 1.0/(a.size()) * total_error;

    return L1;
}