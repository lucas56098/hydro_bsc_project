#include <iostream>
#include <vector>
#include "Point.h"
#include "Cell.h"
#include "Mesh.h"
#include "vmp/VoronoiMesh.h"
#include <random>
#include <algorithm>
using namespace std;


// Function out of vmp main.cpp
// RANDOM POINTS: function to get a sort index
int get_sort_index(Point pt, int sort_grid_size, int sort_scheme) {

    double nr = static_cast<double>(sort_grid_size);

    int index;

    // sort by x-y modulo grid
    if (sort_scheme == 1) {

        if (static_cast<int>(nr * nr * pt.y)%2==0) {

            index = (sort_grid_size - static_cast<int>(pt.x * nr)) + static_cast<int>(nr * nr * pt.y);

        } else {

            index = static_cast<int>(pt.x * nr) + static_cast<int>(nr * nr * pt.y);

        }

    // sort radially outward
    } else if (sort_scheme == 2) {

        index = static_cast<int>((nr*nr)*sqrt((pt.x-0.5)*(pt.x-0.5) + (pt.y-0.5)*(pt.y-0.5)));

    // sort radially inward
    } else if (sort_scheme == 3) {

        index = static_cast<int>((nr*nr)*(sqrt(0.5) - sqrt((pt.x-0.5)*(pt.x-0.5) + (pt.y-0.5)*(pt.y-0.5))));

    // all other numbers -> do not sort
    } else {

        index = 1;

    }

    return index;
}

// Function out of vmp main.cpp
// RANDOM POINTS: generates seed points to use for mesh generation
vector<Point> generate_seed_points(int N, bool fixed_random_seed, double min, int max, int rd_seed, bool sort_pts, int sort_precision, int sort_scheme) {
    vector<Point> points;

    unsigned int random_seed;
    default_random_engine eng;

    // set either fixed or changing random seed
    if (fixed_random_seed) {
        //cout << "specify random_seed: ";
        //cin >> random_seed;
        random_seed = rd_seed;
    } else {
        random_device rd;
        random_seed = rd();

    }

    // define uniform random distribution
    eng = default_random_engine(random_seed);
    uniform_real_distribution<double> distr(min, max);

    // generate random coordinates for Points
    for (int i = 0; i < N; ++i) {
        double x = distr(eng);
        double y = distr(eng);
        points.push_back(Point(x, y));
    }

    // if this is true the points will be sorted
    if (sort_pts) {
        vector<int> indices;
        vector<int> sort_indices;

        // get sort indices
        for (int i = 0; i < points.size(); i++) {
            indices.push_back(get_sort_index(points[i], sort_precision, sort_scheme));
            sort_indices.push_back(i);
        }

        // combine data into pairs
        vector<pair<int, int> > combined;
        for (int i = 0; i < indices.size(); ++i) {
            combined.push_back(make_pair(indices[i], sort_indices[i]));
        }    

        // sort combined data by sort indices
        sort(combined.begin(), combined.end());

        // get sorted_pts
        vector<Point> sorted_pts;
        for (int i = 0; i < combined.size(); i++) {
            sorted_pts.push_back(points[combined[i].second]);
        }

        return sorted_pts;
    }

    return points;
}


// MAIN :  -------------------------------------------------------------------------------------------------------
int main () {
    
    Mesh grid;
    bool cartesian = true;
    int sim_steps = 2000;
    int save_iter = 40;

    // GRID -----------------------------
    if (cartesian) {
        
        // generate simple uniform grid
        grid.generate_uniform_grid(Point(0,0), 45, 45, 1.0/45.0);
    
    } else {

        // genereate seedpoints and vmesh
        vector<Point> pts = generate_seed_points(2000, true, 0, 1, 42, true, 100, 2);
        VoronoiMesh vmesh(pts);
        vmesh.do_point_insertion();

        // generate mesh from vmesh
        grid.generate_from_vmesh(vmesh);

    }
    cout << "grid generated" << endl;

    // INITIAL CONDITIONS ---------------

    grid.cells[1011].Q = 1;

    //for (int i =0; i<1000; i++) {
    //    grid.cells[i].Q = 1;
    //}

    //for (int i = 0; i < 100*100; i++) {
    //    if (i%100 > 30 && i%100 < 70 && i/100 > 30 && i/100 < 70) {
    //        grid.cells[i].Q = 1;
    //    }
    //}

    // SIMULATE change of Q on mesh -----
    for (int i = 0; i<sim_steps; i++) {
        
        // save meshfiles
        if (i%save_iter == 0) {
            grid.save_mesh(i);
            cout << i << endl;
        }

        grid.updateQ();

        // optional: continuous input
        for (int i =0; i<1; i++) {
            grid.cells[1011].Q += 1;
        }
    }

    cout << "done" << endl;

    return 0;    

}


