#ifndef Conway_Cell_h
#define Conway_Cell_h
#include "Cell.h"
#include <iostream>


class Conway_Cell : public Cell {
public:
    Conway_Cell();
    Conway_Cell(Point seedin, std::vector<edge> edgesin);
    ~Conway_Cell();

    int Q;

    double getQ() const {
        return Q;
    }

};

#endif
