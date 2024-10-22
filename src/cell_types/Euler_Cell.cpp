#include "Euler_Cell.h"

Euler_Cell::Euler_Cell() {}

Euler_Cell::Euler_Cell(Point seedin, std::vector<edge> edgesin) {

    seed = seedin;
    edges = edgesin;

    rho = 1;
    u = 0;
    v = 0;
    E = 1;
    gamma = 5.0/3.0;

}

Euler_Cell::~Euler_Cell() {}