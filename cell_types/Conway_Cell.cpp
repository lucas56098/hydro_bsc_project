#include "Conway_Cell.h"

Conway_Cell::Conway_Cell() {}

Conway_Cell::Conway_Cell(Point seedin, std::vector<face> edgesin, int Qin) {

    seed = seedin;
    edges = edgesin;
    Q = Qin;

}

Conway_Cell::~Conway_Cell() {}