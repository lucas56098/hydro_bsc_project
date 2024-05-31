#include <iostream>
#include <cmath>
#include <vector>
#include "VoronoiCell.h"
#include "Point.h"
#include "Halfplane.h"

VoronoiCell::VoronoiCell() {}

VoronoiCell::VoronoiCell(Point in_seed, int in_index) {
    
    index = in_index;
    seed = in_seed;
    halfplanes.reserve(10);
    edges.reserve(10);
    verticies.reserve(10);

}

VoronoiCell::~VoronoiCell() {}

// intersect two halfplanes, add that intersection to the first halfplane
void VoronoiCell::intersect_two_halfplanes(Halfplane &hp1, Halfplane &hp2, vector<intersection> &intersections) {

    // solve linear system of equations for the two lines
    //calculate Determinant D
    double D = hp1.hp_vec.x * hp2.hp_vec.y - hp1.hp_vec.y * hp2.hp_vec.x;
    //calculate Dx
    double Dx = (hp2.midpoint.x - hp1.midpoint.x) * hp2.hp_vec.y - 
                (hp2.midpoint.y -  hp1.midpoint.y) * hp2.hp_vec.x;
    //calculate Dy
    double Dy = hp1.hp_vec.x * (hp2.midpoint.y -  hp1.midpoint.y) -
                hp1.hp_vec.y * (hp2.midpoint.x - hp1.midpoint.x);

    double x;
    double y;


    // if D !=0 -> exact solution x = Dx/D, y = Dy/D
    if (D != 0) {
        x = Dx/D;
        y = Dy/D;
        
        // calculate intersection
        double intersect_pt_x = hp1.midpoint.x + x*hp1.hp_vec.x;
        double intersect_pt_y = hp1.midpoint.y + x*hp1.hp_vec.y;

        // add intersection to hp1
        intersection intersect;
        intersect.intersect_pt = Point(intersect_pt_x, intersect_pt_y);
        intersect.dist_to_midpoint = x;
        intersect.intersecting_with = &hp2;

        intersections.push_back(intersect);


    }
    // if D = 0 -> check if Dx and Dy are 0 -> if: infinite sol, not: no sol
    else if (D == 0)
    {
        if (Dx == 0 && Dy == 0) {
            cout << "infinite solutions: that shouldnt happen" << endl;
            
        } else if (!(hp1.boundary || hp2.boundary)) {
            cout << "no solution while trying to intersect two halfplanes" << endl;
                    
        }
    }
    
    
}

// generate all halfplanes + boundary halfplanes
void VoronoiCell::generate_halfplane_vector(vector<Point> pts,  vector<int> indices) {
    
    // generate boundary halfplanes
    halfplanes.push_back(Halfplane(Point(0.5,0.5), Point(0.5,1.5), -1, -2, true));
    halfplanes.push_back(Halfplane(Point(0.5,0.5), Point(1.5, 0.5), -1, -3,true));
    halfplanes.push_back(Halfplane(Point(0.5,0.5), Point(0.5,-0.5), -1, -4, true));
    halfplanes.push_back(Halfplane(Point(0.5,0.5), Point(-0.5,0.5), -1, -5, true));

    // generate usual halfplanes
    for (int i = 0; i<pts.size(); i++) {
        if (!(index == indices[i])) {
            halfplanes.push_back(Halfplane(seed, pts[i], index, indices[i]));
        }
    }

}

// search for closest point
void VoronoiCell::search_hp_closest_to_seed(Halfplane &first_hp) {

    double dist_min = 42;

    // go through halfplanes and find closest one
    for (int i = 0; i<halfplanes.size(); i++) {
        double dist = sqrt((seed.x - halfplanes[i].midpoint.x)*(seed.x - halfplanes[i].midpoint.x)+
                            (seed.y - halfplanes[i].midpoint.y) * (seed.y - halfplanes[i].midpoint.y));
        if (dist_min > dist) {
            dist_min = dist;
            first_hp = halfplanes[i];
        }
    }
}

// get the signed angle between two vectors (signed means angle from -180deg to +180deg continously, + == right)
double VoronoiCell::get_signed_angle(Point u, Point v) {
    double crossp = u.x*v.y - u.y*v.x;
    double dotp = u.x * v.x + u.y * v.y;

    double angle = -atan2(crossp, dotp);

    return angle;
}

// algorithm to construct the cell
void VoronoiCell::construct_cell(vector<Point> pts, vector<int> indices) {

    // generate all halfplanes
    generate_halfplane_vector(pts, indices);

    // set first and current_hp to hp closest to seed
    Halfplane current_hp;
    search_hp_closest_to_seed(current_hp);
    Halfplane first_hp = current_hp;
    //Point last_vertex_midpoint = Point(42,42);
    int last_vertex_index2 = -42;

    // intitalize last_vertex  for the first time
    Point last_vertex = current_hp.midpoint;

    int counter = 0;
    // step by step generate voronoi cell
    do {

    // determine signed distance between last vertex and current hp_midpoint
    double last_vertex_dist_to_midpoint = (last_vertex.x - current_hp.midpoint.x)*current_hp.hp_vec.x
                                        + (last_vertex.y - current_hp.midpoint.y)*current_hp.hp_vec.y;

    double smallest_pos_dist = 42;
    Halfplane next_hp;
    Point vertex;
    bool need_to_check_for_degeneracy = false;

    vector<intersection> intersections;

    //intersect halfplanes for current_hp
    for (int j = 0; j<halfplanes.size(); j++) {
            if (!(current_hp.index2 == halfplanes[j].index2)) {
                intersect_two_halfplanes(current_hp, halfplanes[j], intersections);
            }
        }

    // find intersection with smallest positive signed distance to last_vertex
    for (int i = 0; i<intersections.size(); i++) {

        // calculate signed relative distance
        double rel_dist =  intersections[i].dist_to_midpoint - last_vertex_dist_to_midpoint;

        bool same_vertex = ((*intersections[i].intersecting_with).index2 == last_vertex_index2);

        // if distance <= smallest_pos_distance and positive replace intersection with closer one
        if (rel_dist <= smallest_pos_dist && rel_dist > 0.0000001 && !same_vertex) {
            // check if there could be a degenerate case
            if (rel_dist - 0.000001 < smallest_pos_dist && rel_dist + 0.000001 > smallest_pos_dist) {
                need_to_check_for_degeneracy = true;
            }
            // update nearest intersection candidate
            smallest_pos_dist = rel_dist;
            next_hp = *intersections[i].intersecting_with;
            vertex = intersections[i].intersect_pt;

        
        }
    }
     
    // if degeneracy possible check for it and choose the correct halfplane
    if (need_to_check_for_degeneracy) {
        vector<Halfplane> deg_hp_list;
        
        for (int i = 0; i<intersections.size(); i++) {
            
            // calculate signed relative distance
            double rel_dist =  intersections[i].dist_to_midpoint - last_vertex_dist_to_midpoint;
            
            // put degenerate cases into vector
            if (rel_dist == smallest_pos_dist) {
                deg_hp_list.push_back(*intersections[i].intersecting_with);
            }


        }

        // calculate the signed angle between current_hp vec and deg_hp_list vec
        double max_angle = 0;

        for (int i = 0; i<deg_hp_list.size(); i++) {

            double angle = get_signed_angle(current_hp.hp_vec, deg_hp_list[i].hp_vec);

            // choose the hp with highest signed angle as the next hp
            if (angle > max_angle) {
                max_angle = angle;
                next_hp = deg_hp_list[i];
            }
        }
    }
    
    
    // when done save current halfplane as edge and next closest intersection as vertex
    edges.push_back(current_hp);
    verticies.push_back(vertex);


    last_vertex_index2 = current_hp.index2;

    // update last vertex and current hp
    last_vertex = vertex;
    current_hp = next_hp;
    counter += 1;

    // continue with the steps above until the half plane is the same as the first one

    } while (!(first_hp.midpoint.x == current_hp.midpoint.x && first_hp.midpoint.y == current_hp.midpoint.y) && (10000 > counter));

    if (counter >= 10000) {
        cout << "failed to generate cell" << endl;
    }

    // some memory management
    halfplanes.clear();

}

// check all vertecies of the cell for equidistance conditions
bool VoronoiCell::check_equidistance_condition(vector<Point> seeds) {

    bool correct_cell = true;

    // check if every vertex has at least three equidistant points
    for (int i = 0; i < verticies.size(); i++) {

        // only check vertex if it is not on boundary within error margin
        bool condition = verticies[i].x > 0.000000000000001 && verticies[i].y > 0.000000000000001 && verticies[i].x < 0.999999999999999 && verticies[i].y < 0.999999999999999;

        if (condition) {

            int nr_equidist_pts = 0;

            double min_dist = 42;

            // determine minimum distance
            for (int j = 0; j < seeds.size(); j++) {

                // calc distance between seed and vertex
                double dist = sqrt((seeds[j].x - verticies[i].x)*(seeds[j].x - verticies[i].x) +
                                    (seeds[j].y - verticies[i].y)*(seeds[j].y - verticies[i].y));

                if (dist < min_dist) {
                    min_dist = dist;
                }
            }

            // search for seeds with minimum distance
            for (int j = 0; j < seeds.size(); j++) {

                // calc distance between seed and vertex
                double dist = sqrt((seeds[j].x - verticies[i].x)*(seeds[j].x - verticies[i].x) +
                                    (seeds[j].y - verticies[i].y)*(seeds[j].y - verticies[i].y));


                // if seed has minimum distance add to equidistant points nr
                if (min_dist - 0.000000000000001 < dist && dist < min_dist + 0.000000000000001) {
                    nr_equidist_pts += 1;
                }
            }

            // check equidistance condition (= all vertecies not on boundary must have at least three equidistant seeds)
            if (nr_equidist_pts < 3) {
                correct_cell = false;
            }

        }

    }

    return correct_cell;
}

// function to calculate the area of the voronoi cell 
double VoronoiCell::get_area() {

    double area = 0;

    for (int i = 0; i < verticies.size(); i++) {
        area += -0.5 * ((verticies[i].y + verticies[(i+1)%(verticies.size())].y) * 
                       (verticies[i].x - verticies[(i+1)%(verticies.size())].x));
    }
    return area;
}

long long VoronoiCell::calculate_cell_memory(bool use_capacity) {

    long long total_size;
    if (use_capacity) {
        total_size = sizeof(VoronoiCell) + sizeof(Halfplane)*(halfplanes.capacity() + edges.capacity()) + sizeof(Point)*verticies.capacity();
    } else {
        total_size = sizeof(VoronoiCell) + sizeof(Halfplane)*(halfplanes.size() + edges.size()) + sizeof(Point)*verticies.size();
    }

    return total_size;
}


Point VoronoiCell::get_centroid() {

    double A = get_area();

    double sum_x = 0;
    double sum_y = 0;

    for (int i = 0; i<verticies.size(); i++) {
        sum_x += (verticies[i].x + verticies[(i+1)%verticies.size()].x)*(verticies[i].x * verticies[(i+1)%verticies.size()].y - verticies[(i+1)%verticies.size()].x * verticies[i].y);
        sum_y += (verticies[i].y + verticies[(i+1)%verticies.size()].y)*(verticies[i].x * verticies[(i+1)%verticies.size()].y - verticies[(i+1)%verticies.size()].x * verticies[i].y);
    }

    double C_x = -sum_x/(6*A);
    double C_y = -sum_y/(6*A);

    return Point(C_x, C_y);
}