#include "SWE_Cell.h"

SWE_Cell::SWE_Cell() {}

SWE_Cell::SWE_Cell(Point seedin, std::vector<edge> edgesin) {

    seed = seedin;
    edges = edgesin;
    h = 1;
    u = 0;
    v = 0;

}

SWE_Cell::~SWE_Cell() {}