#include "Q_Cell.h"

Q_Cell::Q_Cell() {}

Q_Cell::Q_Cell(Point seedin, std::vector<face> edgesin) {

    seed = seedin;
    edges = edgesin;
    Q = 0;

}

Q_Cell::~Q_Cell() {}