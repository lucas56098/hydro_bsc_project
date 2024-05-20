#include "Q_Cell.h"

Q_Cell::Q_Cell() {}

Q_Cell::Q_Cell(Point seedin, std::vector<face> edgesin) {

    seed = seedin;
    edges = edgesin;
    Q = 0;
    cell_type = 1;

}

Q_Cell::~Q_Cell() {}