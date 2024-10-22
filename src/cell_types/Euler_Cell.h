#ifndef Euler_Cell_h
#define Euler_Cell_h
#include "Cell.h"
#include <iostream>

// cell to solve the euler equations
class Euler_Cell : public Cell {
public:
    Euler_Cell();
    Euler_Cell(Point seedin, std::vector<edge> edgesin);
    ~Euler_Cell();

    double rho;     // density
    double u;       // v.x
    double v;       // v.y
    double E;       // energy density
    double gamma;   // adiabatic index (is assumed to be constant an the same for every cell)

    double get_rho() const {
        return rho;
    }

    double get_u() const {
        return u;
    }

    double get_v() const {
        return v;
    }

    double get_E() const {
        return E;
    }

    double get_gamma() const {
        return gamma;
    }

    double get_P() const {
        return (gamma - 1) * (E - (0.5 * rho * (u*u + v*v)));
    }

};

#endif
