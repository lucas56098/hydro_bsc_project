#ifndef Q_Cell_h
#define Q_Cell_h
#include "Cell.h"
#include <iostream>


class Q_Cell : public Cell {
public:
    Q_Cell();
    Q_Cell(Point seedin, std::vector<edge> edgesin);
    ~Q_Cell();

    double Q;

    double getQ() const {
        return Q;
    }

};

#endif
