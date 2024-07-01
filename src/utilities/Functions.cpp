#include "Functions.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>


// ANALYTIC SOLUTIONS --------------------------------------------------------------------
// ADVECTION - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Returns analytic solution of advecting circle for a given seedpoint
double advecting_circle(Point seed, double t, Point v, Point p0, double r) {

    // calc midpoint of circle
    Point p = Point(p0.x + t * v.x, p0.y + t * v.y);

    // if dist to point < radius return 1
    double dist = sqrt((p.x - seed.x)*(p.x - seed.x) + (p.y - seed.y)*(p.y - seed.y));
    if (dist <= r) {
        return 1;
    }
    return 0;

}


// Returns analytic solution of advecting stepfunc in 1D for a given seedpoint
double advecting1D_stepfunc(Point seed, double t, double v, double a0, double b0) {

    // if x value between a(t), b(t) set Q to 1
    if (seed.x >= a0 + v*t && seed.x <= b0 + v*t) {
        return 1;
    }
    return 0;

}


// SHALLOW WATER EQUATIONS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Returns analytic solution of SWE dam break 1D
// note hm needs to be manually calculated before e.g. using method in vis_tk solving for that equations root.
// of course thats bad practice in general but i dont see any benefit for my project in writing a sophisticated root finder.
// hm default value is correct for hl = 2, hm = 1
double swe1D_dam_break(Point seed, double t, double hl, double hr, double hm, double g, double x0) {

    // calc wave speeds
    double cl = sqrt(g*hl);
    double cr = sqrt(g*hr);    
    double cm = sqrt(g*hm);

    // calc positions of shocks/rarefications
    double x = seed.x;
    double xa = x0 - cl * t;                // leftmost part of rarefication
    double xb = x0 + t * (2*cl - 3*cm);     // rightmost part of rarefication
    double xc = x0 + t * ((2 * cm * cm * (cl - cm))/(cm*cm - cr*cr)); // right shock

    // return height according to characteristics
    if (x<xa) {
        return hl;
    } else if (xa < x && x < xb) {
        return (4)/(9*g) * (cl - ((x - x0)/(2*t)))*(cl - ((x - x0)/(2*t)));
    } else if (xb < x && x < xc) {
        return hm;
    } else if (xc < x) {
        return hr;
    }
    return INFINITY;
}


// OTHER ------------------------------------------------------------------------------
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

// formats time from seconds into mm:ss
string format_time(double seconds) {

    // get minutes and seconds
    int minutes = static_cast<int>(seconds) / 60;
    int sec = static_cast<int>(seconds) % 60;
    
    // put them into string
    stringstream ss;
    ss << setfill('0') << setw(2) << minutes << ":" << setfill('0') << setw(2) << sec;
    return ss.str();
}