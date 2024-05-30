#include <vector>
#include "Point.h"
#include "Halfplane.h"

#ifndef VoronoiCell_h
#define VoronoiCell_h

struct intersection
    {
        Point intersect_pt;
        Halfplane* intersecting_with;
        double dist_to_midpoint; //distance signed relative to half_plane_vec
    };

class VoronoiCell {

public:
    VoronoiCell();
    VoronoiCell(Point in_seed, int index);
    ~VoronoiCell();
    int index;
    Point seed;
    vector<Halfplane> halfplanes;
    vector<Halfplane> edges;
    vector<Point> verticies;
    void intersect_two_halfplanes(Halfplane &hp1, Halfplane &hp2, vector<intersection> &intersections);
    void construct_cell(vector<Point> pts, vector<int> indices);
    bool check_equidistance_condition(vector<Point> seeds);
    double get_area();
    void generate_halfplane_vector(vector<Point> pts, vector<int> indices);
    double get_signed_angle(Point u, Point v);
    long long calculate_cell_memory(bool use_capacity);
private:
    void search_hp_closest_to_seed(Halfplane &first_hp);


};

#endif