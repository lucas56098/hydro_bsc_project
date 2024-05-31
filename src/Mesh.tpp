#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "Mesh.h"
#include "cell_types/Cell.h"
#include "cell_types/Point.h"
#include "cell_types/Q_Cell.h"
#include "cell_types/Conway_Cell.h"
#include "vmp/VoronoiMesh.h"
#include "utilities/Functions.h"
#include <algorithm>
#include <type_traits>
#include <random>

template <typename CellType>
Mesh<CellType>::Mesh() {}


// PRIVATE Helper functions for Point Generation --------------------------------------------------
// Function out of vmp main.cpp
// RANDOM POINTS: function to get a sort index
template <typename CellType>
int Mesh<CellType>::get_sort_index(Point pt, int sort_grid_size, int sort_scheme) {

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
template <typename CellType>
vector<Point> Mesh<CellType>::generate_seed_points(int N, bool fixed_random_seed, double min, int max, int rd_seed, bool sort_pts, int sort_precision, int sort_scheme) {
    vector<Point> points;

    unsigned int random_seed;
    default_random_engine eng;

    // set either fixed or changing random seed
    if (fixed_random_seed) {
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


// GRID GENERATION: -------------------------------------------------------------------------------
// calls the generate Mesh functions depending on specified options (cartesian, 1D/2D, N_row)
template <typename CellType>
void Mesh<CellType>::generate_grid(bool cartesian, bool is_1D, int N_row, int lloyd_iterations) {

    if (cartesian) {
        // generate cartesian mesh
        if (is_1D) {
            // do it in 1D
            this->generate_uniform_grid1D(Point(0, 0), N_row, 1.0/static_cast<double>(N_row));
            is_cartesian = true;
        } else {
            // do it in 2D
            this->generate_uniform_grid2D(Point(0, 0), N_row, N_row, 1.0/static_cast<double>(N_row), 1.0/static_cast<double>(N_row));
            is_cartesian = true;
        }
    } else {
        // generate voronoi mesh
        if (is_1D) {
            // do it in 1D
            vector<Point> pts = generate_seed_points(N_row, true, 0, 1, 42, true, 100, 1);
            this->generate_vmesh1D(pts);
            is_cartesian = false;
        } else {
            // do it in 2D
            vector<Point> pts = generate_seed_points(N_row * N_row, true, 0, 1, 42, true, 100, 1);
            this->generate_vmesh2D(pts, lloyd_iterations);
            is_cartesian = false;
        }

    }
    cout << "grid generated" << endl;


}


// generates a uniform grid with all the neighbour relations and so on
template <typename CellType>
void Mesh<CellType>::generate_uniform_grid2D(Point start, int n_hor, int n_vert, double distx, double disty) {

    for (int b = 0; b<n_vert; b++) {
        for (int a = 0; a<n_hor; a++) {

            // set the correct seed
            Point seedin(start.x + 0.5*distx + a*distx, start.y + 0.5*disty + b*disty);
            
            // start to define the faces
            vector<face> edgesin;
            face f0; face f1; face f2; face f3;

            // set a and b for the faces
            f0.a = Point(start.x + a*distx, start.y + b*disty);
            f0.b = Point(start.x + a*distx, start.y + b*disty + disty);
            f1.a = Point(start.x + a*distx, start.y + b*disty + disty);
            f1.b = Point(start.x + a*distx + distx, start.y + b*disty + disty);
            f2.a = Point(start.x + a*distx + distx, start.y + b*disty + disty);
            f2.b = Point(start.x + a*distx + distx, start.y + b*disty);
            f3.a = Point(start.x + a*distx + distx, start.y + b*disty);
            f3.b = Point(start.x + a*distx, start.y + b*disty);

            // set correct boundary flags
            f0.is_boundary = false; f1.is_boundary = false; f2.is_boundary = false; f3.is_boundary = false;
            if (a == 0) {f0.is_boundary = true;}
            if (a == n_hor -1) {f2.is_boundary = true;}
            if (b == 0) {f3.is_boundary = true;}
            if (b == n_vert -1) {f1.is_boundary = true;}

            // push faces in edges vector
            edgesin.push_back(f0);
            edgesin.push_back(f1);
            edgesin.push_back(f2);
            edgesin.push_back(f3);

            // set face length to dist

            edgesin[0].length = disty;
            edgesin[1].length = distx;
            edgesin[2].length = disty;
            edgesin[3].length = distx;

            // push back new cell in cells vector
            cells.emplace_back(seedin, edgesin);

        }
    }

    // now that all cells exist we define the neighbour relations
    for (int i = 0; i<cells.size(); i++) {

        cells[i].volume = distx * disty;

        for (int j = 0; j<cells[i].edges.size(); j++) {
            if (cells[i].edges[j].is_boundary == false) {
                if (j == 0) {
                    cells[i].edges[j].neighbour = &cells[i-1];
                } else if (j == 1) {
                    cells[i].edges[j].neighbour = &cells[i+n_hor];
                } else if (j == 2) {
                    cells[i].edges[j].neighbour = &cells[i+1];
                } else if (j == 3) {
                    cells[i].edges[j].neighbour = &cells[i-n_hor];
                }
                
            }
        }
    }
}


// generates vmesh using vmp and converts it into data usable for this mesh type
template <typename CellType>
void Mesh<CellType>::generate_vmesh2D(vector<Point> pts, int lloyd_iterations) {

    VoronoiMesh vmesh(pts);
    vmesh.do_point_insertion();

    // do iterations of lloyds algorithm as preprocessing
    for (int i = 0; i<lloyd_iterations; i++) {

            vector<Point> centroids;
        for (int i = 0; i<vmesh.vcells.size(); i++) {
            centroids.push_back(vmesh.vcells[i].get_centroid());
        }

        vmesh = VoronoiMesh(centroids);
        vmesh.do_point_insertion();

    }

    // loop through all cells to set everything but neighbour relations
    for (int i = 0; i<vmesh.vcells.size(); i++) {

        // set seed and define edge vector
        Point seedin(vmesh.vcells[i].seed.x, vmesh.vcells[i].seed.y);
        vector<face> edgesin;

        for (int j = 0; j<vmesh.vcells[i].edges.size(); j++) {

            face f;

            // set face quantities
            f.a = vmesh.vcells[i].verticies[((vmesh.vcells[i].edges.size()-1) + j)%vmesh.vcells[i].edges.size()];
            f.b = vmesh.vcells[i].verticies[j];
            f.length = sqrt((f.a.x - f.b.x)*(f.a.x - f.b.x) + (f.a.y - f.b.y)*(f.a.y - f.b.y));
            f.is_boundary = (vmesh.vcells[i].edges[j].index1 < 0);

            edgesin.push_back(f);

        }

        // push back new cell in cells
        cells.emplace_back(seedin, edgesin);

    }

    // loop through everything and set neighbour relations
    for (int i = 0; i<vmesh.vcells.size(); i++) {

        cells[i].volume = vmesh.vcells[i].get_area();

        for (int j = 0; j<vmesh.vcells[i].edges.size(); j++) {

            // neighbour index as used in vmesh
            int neigbour_index;
            neigbour_index = vmesh.vcells[i].edges[j].index2;

            // exclude boundaries
            if (neigbour_index >= 0) {

                // neighbour as defined before over index now with adress
                cells[i].edges[j].neighbour = &cells[neigbour_index];

            }

        }

    }


}


// generates a 1D uniform grid with all the neighbour relations and so on
template <typename CellType>
void Mesh<CellType>::generate_uniform_grid1D(Point start, int n, double dist) {
    generate_uniform_grid2D(start, n, 1, dist, 1);
}


// generates a 1D voronoi mesh, only works if points are between 0 and 1
template <typename CellType>
void Mesh<CellType>::generate_vmesh1D(vector<Point> pts) {

    for (int i = 0; i < pts.size(); i++) {
        pts[i].y = 0.5;
    }

    vector<Point> sorted_pts = pts;
    sort(sorted_pts.begin(), sorted_pts.end(), [](const Point& a, const Point& b) {
        return a.x < b.x;
    });

    for (int i = 0; i < sorted_pts.size(); i++) {

        // set the correct seed
        Point seedin = sorted_pts[i];

        // start to define the faces
        vector<face> edgesin;
        face f0; face f1; face f2; face f3;

        // set a and b for the faces and boundary flags
        f0.is_boundary = false;
        f1.is_boundary = true;
        f2.is_boundary = false;
        f3.is_boundary = true;

        double distl;
        double distr;

        // manually set left part of faces
        if (i == 0) {
            // we are in the leftmost cell
            f0.is_boundary = true;
            f0.a = Point(0, 0);
            f0.b = Point(0, 1);
            f1.a = Point(0, 1);
            f3.b = Point(0,0);
            distl = seedin.x;
        } else {
            // just a random cell
            f0.a = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 0);
            f0.b = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 1);
            f1.a = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 1);
            f3.b = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 0);
            distl = sorted_pts[i].x - sorted_pts[i-1].x;
        }
        
        // manually set rigth part of faces
        if (i == sorted_pts.size()-1) {
            // we are in the rightmost cell
            f2.is_boundary = true;
            f1.b = Point(1, 1);
            f2.a = Point(1, 1);
            f2.b = Point(1,0);
            f3.a = Point(1,0);
            distr = 1 - seedin.x;
        } else {
            // just a random cell
            f1.b = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 1);
            f2.a = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 1);
            f2.b = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 0);
            f3.a = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 0);
            distr = sorted_pts[i+1].x - sorted_pts[i].x;
        }

        // push back faces
        edgesin.push_back(f0);
        edgesin.push_back(f1);
        edgesin.push_back(f2);
        edgesin.push_back(f3);

        // set correct lengths
        edgesin[0].length = 1;
        edgesin[2].length = 1;
        edgesin[1].length = distl + distr;
        edgesin[3].length = distl + distr;

        // push back new cell in cells vector
        cells.emplace_back(seedin, edgesin);

    }

    // now that all cells exist we define the neighbour relations
    for (int i = 0; i < cells.size(); i++) {

        cells[i].volume = cells[i].edges[0].length * cells[i].edges[1].length;

        if (i != 0) {
            cells[i].edges[0].neighbour = &cells[i-1];
        }
        if (i != cells.size()-1) {
            cells[i].edges[2].neighbour = &cells[i+1];
        }

    }

}


// SET INITIAL CONDITIONS -------------------------------------------------------------------------
// sets the inital value for the first N cells of the cell vector to value
template <typename CellType>
void Mesh<CellType>::initialize_Q_cells(int a, int b, double value, int step) {

    if(cells.size()< b) {
        cerr << "ERROR: initalize_cells(a, b, value), tried to initalize more cells then there are! b = " << b  << " > cells.size() =" << cells.size() << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = a; i < b; i+=step) {
        cells[i].Q = value;
    }

}

// sets the initial condition according to given function
template <typename CellType>
void Mesh<CellType>::initalize_Q_circle(Point p0, double r) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Q_Cell> == false) {
        cerr << "initalize_Q_circle called with wrong cell type, you must use Q_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < cells.size(); i++) {
        cells[i].Q = advecting_circle(cells[i].seed, 0, Point(0, 0), p0, r);
    }

}

// SAVE DATA TO FILE ------------------------------------------------------------------------------
// saves mesh into csv file readable for python script
template <typename CellType>
void Mesh<CellType>::save_mesh(int file_nr, string name, double dt) {

    string filename;
    filename = "../src/files/" + name + to_string(file_nr) + ".csv"; 

    ofstream output_file(filename);
    
    // get maximum edge number for column correction later on
    int max_edge_nr = 0;
    for (int i = 0; i<cells.size(); i++) {
        if (cells[i].edges.size()>=max_edge_nr) {
            max_edge_nr = cells[i].edges.size();
        }
    }

    // save the mesh in the following format: ax1, ay1, bx1, by1 ; ax2, ... ; ..; | Q
    output_file << "seed.x, seed.y | a.x, a.y, b.x, b.y ; a.x ... ; | Quantities" << endl;

    for (int i = 0; i<cells.size(); i++) {

        output_file << cells[i].seed.x << "," << cells[i].seed.y << "|";

        for (int j = 0; j<cells[i].edges.size(); j++) {

            output_file << cells[i].edges[j].a.x << ','
                        << cells[i].edges[j].a.y
                        << ";";
        }

        // correct for empty columns such that further values are always at the same column
        for (int a = 0; a<(max_edge_nr-cells[i].edges.size()); a++) {
            output_file << ";";
        }

        output_file << "|";

        output_file << file_nr * dt << ",";

        // if cell type is Q_cell or Conway Cell save Q
        if constexpr (is_same_v<CellType, Q_Cell> || is_same_v<CellType, Conway_Cell>) {
            output_file << cells[i].Q;
        }

        output_file << endl;

    }

    output_file.close();

}


// function to save the change in the summed up Q value over the whole grid
template <typename CellType>
void Mesh<CellType>::save_Q_diff(double t, bool reset_file, bool is_density) {

    double total_Q = 0;

    // calculate total_Q summed up over mesh 
    //(depending on wether it is a density multiply with volume) 
    if (is_density) {
        for (int i = 0; i < cells.size(); i++) {
            total_Q += cells[i].getQ()*cells[i].volume;
        }
    } else {
        for (int i = 0; i < cells.size(); i++) {
            total_Q += cells[i].getQ();
        }
    }
    
    // open new file or in append mode
    ofstream output_file;
    if (reset_file) {
        output_file = ofstream("../src/files/total_Q_diff.csv");
        total_Q_initial = total_Q;
    } else {
        output_file = ofstream("../src/files/total_Q_diff.csv", ios::app);
    }

    // write caluclated difference in file
    output_file << t << "," << total_Q-total_Q_initial << endl;

    output_file.close();
    

}

// function to calculate and store the L1 error of an advecting circle, requires Q_cells
template <typename CellType>
void Mesh<CellType>::save_L1_adv_circle(double t, bool reset_file, Point v, Point p0, double r) {

    vector<double> Qs_num;
    vector<double> Qs_ana;
    Qs_num.reserve(cells.size());
    Qs_ana.reserve(cells.size());

    for (int i = 0; i<cells.size(); i++) {
        Qs_num.push_back(cells[i].Q);
        Qs_ana.push_back(advecting_circle(cells[i].seed, t, v, p0, r));
    }

    // open new file or in append mode
    ofstream output_file;
    if (reset_file) {
        output_file = ofstream("../src/files/L1_error.csv");
    } else {
        output_file = ofstream("../src/files/L1_error.csv", ios::app);
    }

    // write caluclated difference in file
    output_file << t << "," << L1_error(Qs_num, Qs_ana) << endl;
    output_file.close();
}


// function to calculate and store the L1 error of an advecting circle, requires Q_cells
template <typename CellType>
void Mesh<CellType>::save_L1_adv_1Dstepfunc(double t, bool reset_file, double v, double a0, double b0) {

    vector<double> Qs_num;
    vector<double> Qs_ana;
    Qs_num.reserve(cells.size());
    Qs_ana.reserve(cells.size());

    for (int i = 0; i<cells.size(); i++) {
        Qs_num.push_back(cells[i].Q);
        Qs_ana.push_back(advecting1D_stepfunc(cells[i].seed, t, v, a0, b0));
    }

    // open new file or in append mode
    ofstream output_file;
    if (reset_file) {
        output_file = ofstream("../src/files/L1_error.csv");
    } else {
        output_file = ofstream("../src/files/L1_error.csv", ios::app);
    }

    // write caluclated difference in file
    output_file << t << "," << L1_error(Qs_num, Qs_ana) << endl;
    output_file.close();
}



template <typename CellType>
Mesh<CellType>::~Mesh() {}