#ifndef SWE_Cell_h
#define SWE_Cell_h
#include "Cell.h"
#include <iostream>

// cell to solve the shallow water equation
// assumptions: rho = const, flat bottom of the sea
class SWE_Cell : public Cell {
public:
    SWE_Cell();
    SWE_Cell(Point seedin, std::vector<edge> edgesin);
    ~SWE_Cell();

    double h; // column height
    double u; // v.x
    double v; // v.y

    double get_h() const {
        return h;
    }

    double get_u() const {
        return u;
    }

    double get_v() const {
        return v;
    }

};

#endif
