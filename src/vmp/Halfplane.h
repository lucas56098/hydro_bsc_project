#include "Point.h"
#include <vector>
using namespace std;

#ifndef Halfplane_h
#define Halfplane_h

class Halfplane {

public:
    Halfplane();
    Halfplane(Point inseed1, Point inseed2, int index_1, int index_2);
    Halfplane(Point inseed1, Point inseed2, int index_1, int index_2, bool is_boundary);
    ~Halfplane();
    Point midpoint;
    Point hp_vec;
    bool boundary;
    int index1;
    int index2;
    void calc_half_plane_vec(Point &seed1, Point &seed2);
    void calc_midpoint(Point &seed1, Point &seed2);
    


private:

};

#endif