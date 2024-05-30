#include "Cell.h"
#include "Point.h"
#include <iostream>

Cell::Cell() {}

Cell::Cell(Point seedin, vector<face> edgesin) {

    seed = seedin;
    edges = edgesin;

}

Cell::~Cell() {}