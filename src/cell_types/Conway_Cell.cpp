#include "Conway_Cell.h"

Conway_Cell::Conway_Cell() {}

Conway_Cell::Conway_Cell(Point seedin, std::vector<edge> edgesin) {

    seed = seedin;
    edges = edgesin;
    Q = 0;

}

Conway_Cell::~Conway_Cell() {}