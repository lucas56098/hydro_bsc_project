#include "Point.h"
#include <vector>
using namespace std;


#ifndef Cell_h
#define Cell_h

class Cell;

struct face
    {
        Point a;
        Point b;
        double length;
        double deltaQ;
        bool is_boundary;
        Cell* neighbour;
    };

class Cell {

public:
    Cell();
    Cell(Point seedin, vector<face> edgesin);
    ~Cell();

    Point seed;
    vector<face> edges;
    double Q;

    void calc_deltaQs();
    void apply_deltaQs();


};

#endif