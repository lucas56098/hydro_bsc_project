#include "Cell.h"
#include "Point.h"
#include <iostream>

Cell::Cell() {}

Cell::Cell(Point seedin, vector<face> edgesin) {

    seed = seedin;
    edges = edgesin;

}

// function to calculate flux of Q based on current grid
void Cell::calc_deltaQs() {

    for (int i = 0; i<edges.size(); i++) {
        if (edges[i].is_boundary == false) {
            edges[i].deltaQ = 0.5 * (edges[i].neighbour->Q - Q) * edges[i].length;
        } else {
            edges[i].deltaQ = 0;        // solid boundaries e.g. no flux through them
        }
    }

}

// apply calculated fluxes to grid
void Cell::apply_deltaQs() {

    for (int i = 0; i<edges.size(); i++) {
        Q += edges[i].deltaQ;
    }

}

Cell::~Cell() {}