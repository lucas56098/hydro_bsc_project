#include "Point.h"
#include <vector>
using namespace std;


#ifndef Cell_h
#define Cell_h

class Cell;

struct face
    {
        Point a;
        Point b;
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
        throw std::runtime_error("Q is not defined for base Cell");
    }

};

#endif