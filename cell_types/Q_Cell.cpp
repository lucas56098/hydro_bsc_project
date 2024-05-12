#include "Q_Cell.h"

Q_Cell::Q_Cell() {}

Q_Cell::Q_Cell(Point seedin, std::vector<face> edgesin, double Qin) {

    seed = seedin;
    edges = edgesin;
    Q = Qin;

}

Q_Cell::~Q_Cell() {}