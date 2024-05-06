#include <vector>
#include <cmath>
#include "Halfplane.h"


Halfplane::Halfplane() {}

Halfplane::Halfplane(Point inseed1, Point inseed2, int index_1, int index_2) {
    
    boundary = false;

    index1 = index_1;
    index2 = index_2;

    calc_half_plane_vec(inseed1, inseed2);
    calc_midpoint(inseed1, inseed2);
}

Halfplane::Halfplane(Point inseed1, Point inseed2, int index_1, int index_2, bool is_boundary) {
    
    boundary = is_boundary;

    index1 = index_1;
    index2 = index_2;

    calc_half_plane_vec(inseed1, inseed2);
    calc_midpoint(inseed1, inseed2);
}

Halfplane::~Halfplane() {}

// returns normalized vector orthogonal to vector seed2-seed1
void Halfplane::calc_half_plane_vec(Point &seed1, Point &seed2) {

    double hp_vec_x = seed2.y - seed1.y;
    double hp_vec_y = - (seed2.x - seed1.x);

    double norm = sqrt(hp_vec_x * hp_vec_x + hp_vec_y * hp_vec_y);

    hp_vec = Point(hp_vec_x/norm, hp_vec_y/norm);
}

// returns midpoint between the two seeds
void Halfplane::calc_midpoint(Point &seed1, Point &seed2) {
    
    double x_mid = 0.5 * (seed1.x + seed2.x);
    double y_mid = 0.5 * (seed1.y + seed2.y);

    midpoint =  Point(x_mid, y_mid);
}
