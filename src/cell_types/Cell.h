#include "Point.h"
#include <vector>
using namespace std;

#ifndef Cell_h
#define Cell_h

class Cell;

// structure to store a face out of the edges of the cell
struct face
    {
        // start and endpoint of the face
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
    Cell(Point seedin, vector<face> edgesin);
    ~Cell();

    Point seed;
    vector<face> edges;
    double volume;

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

};

#endif