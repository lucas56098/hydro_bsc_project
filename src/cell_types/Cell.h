#include "Point.h"
#include <vector>
#include "../Eigen/Dense"
using namespace std;
using namespace Eigen;

#ifndef Cell_h
#define Cell_h

class Cell;

// structure to store a edge out of the edges of the cell
struct edge
    {
        // start and endpoint of the edge
        Point a;
        Point b;

        // other
        double length;
        bool is_boundary;
        Cell* neighbour;
    };

class Cell {

public:
    Cell();
    Cell(Point seedin, vector<edge> edgesin);
    ~Cell();

    Point seed;
    Point centroid;
    vector<edge> edges;
    double volume;
    int index;

    virtual double getQ() const {
        throw runtime_error("Q is not defined for base Cell");
    }

    virtual double get_h() const {
        throw runtime_error("h is not defined for base Cell");
    }

    virtual double get_u() const {
        throw runtime_error("u is not defined for base Cell");
    }

    virtual double get_v() const {
        throw runtime_error("v is not defined for base Cell");
    }

    virtual double get_rho() const {
        throw runtime_error("rho is not defined for base Cell");
    }

    virtual double get_E() const {
        throw runtime_error("E is not defined for base Cell");
    }

    virtual double get_gamma() const {
        throw runtime_error("adiabatic index is not defined for base Cell");
    }

    virtual double get_P() const {
        throw runtime_error("P is not defined for base Cell");
    }

    virtual VectorXd get_coeff() const {
        throw runtime_error("coeff is not defined for base Cell");
    }


};

#endif