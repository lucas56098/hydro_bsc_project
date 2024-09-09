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
#include "cell_types/SWE_Cell.h"
#include "cell_types/Euler_Cell.h"
#include "vmp/VoronoiMesh.h"
#include "utilities/Functions.h"
#include <algorithm>
#include <type_traits>
#include <random>

template <typename CellType>
Mesh<CellType>::Mesh(int N_bfunc) {
    N_basisfunc = N_bfunc;
}


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
    
    /*
    // optional point density for KH-instability
    uniform_real_distribution<int> distra(min, max);
    normal_distribution<double> distr1(0.7, 0.03);
    normal_distribution<double> distr2(0.3, 0.03);


    // generate random coordinates for Points
    for (int i = 0; i < N; ++i) {
        double x = distr(eng);

        double a;
        if (distra(eng) > 0) {
            a = distr1(eng);
        } else {
            a = distr2(eng);
        }

        double y = std::max(std::min(a, 0.9999), 0.0001) + 0.0009999999 * distr(eng);

        points.push_back(Point(x, y));
    }*/
    


    // optional point density for quad shock
    uniform_real_distribution<int> distra(1, 2);
    for (int i = 0; i < N; i++) {
        double x;
        double y;
        double scalea = 0.1*distr(eng) + 0.65;
        double scaleb = 0.2*distr(eng) + 0.65;
        if (distra(eng) % 3 != 0) {
            x = scalea*distr(eng) + 0.20;
            y = scaleb*distr(eng) + 0.05;
        } else {
            x = distr(eng);
            y = distr(eng);
        }

        points.push_back(Point(x,y));
    }

    /*
    // optional point density for circle
    normal_distribution<double> distr2(0.5, 0.08);

    // generate random coordinates for Points
    for (int i = 0; i < N; ++i) {
        double x = distr(eng);

        double a = distr2(eng);


        double y = std::max(std::min(a, 0.9999), 0.0001) + 0.0009999999 * distr(eng);

        points.push_back(Point(x, y));
    }*/
    


    // generate random coordinates for Points
    //for (int i = 0; i < N; ++i) {
    //    double x = distr(eng);
    //    double y = distr(eng);

    //    points.push_back(Point(x, y));
    //}

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
// calls the generate Mesh functions depending on specified options (cartesian, 1D/2D, N_row, optional lloyd preprocessing, repeating boundary conditions)
template <typename CellType>
void Mesh<CellType>::generate_grid(bool cartesian, bool is_1D, int N_row, int lloyd_iterations, bool repeating, bool structure) {

    if (cartesian) {
        // generate cartesian mesh
        if (is_1D) {
            // do it in 1D
            this->generate_uniform_grid1D(Point(0, 0), N_row, 1.0/static_cast<double>(N_row), repeating);
            is_cartesian = true;
        } else {
            // do it in 2D
            this->generate_uniform_grid2D(Point(0, 0), N_row, N_row, 1.0/static_cast<double>(N_row), 1.0/static_cast<double>(N_row), repeating);
            is_cartesian = true;
        }
    } else {
        // generate voronoi mesh
        if (is_1D) {
            // do it in 1D
            vector<Point> pts = generate_seed_points(N_row, true, 0, 1, 42, true, 100, 1);
            this->generate_vmesh1D(pts, repeating);
            is_cartesian = false;
        } else {
            // do it in 2D
            vector<Point> pts = generate_seed_points(N_row * N_row, true, 0, 1, 42, true, 100, 1);
            if (lloyd_iterations != 0) {do_lloyd_iterations(&pts, lloyd_iterations);};
            int nr;
            if (structure) {nr = add_struct(&pts, 0.0001, 0.01, "struct");}

            this->generate_vmesh2D(pts, repeating, !structure);
            is_cartesian = false;

            if (structure) {
                for (int i = 0; i < nr; i++) {
                    make_cell_boundary_cell(i);
                }
            }

        }

    }
    for (int i = 0; i<this->cells.size(); i++) {
        this->cells[i].index = i;
    }
    cout << "grid generated" << endl;
    


}


// generates a uniform grid with all the neighbour relations and so on
template <typename CellType>
void Mesh<CellType>::generate_uniform_grid2D(Point start, int n_hor, int n_vert, double distx, double disty, bool repeating) {

    n_horizontal = n_hor;
    n_vertical = n_vert;

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
            if (a == 0 && repeating == false) {f0.is_boundary = true;}
            if (a == n_hor -1 && repeating == false) {f2.is_boundary = true;}
            if (b == 0 && repeating == false) {f3.is_boundary = true;}
            if (b == n_vert -1 && repeating == false) {f1.is_boundary = true;}

            if (n_vert == 1) {
                f1.is_boundary = true;
                f3.is_boundary = true;
            }

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
            if constexpr (is_same_v<CellType, DG_Q_Cell> == true) {
                cells.emplace_back(seedin, edgesin, N_basisfunc);
            } else {
                cells.emplace_back(seedin, edgesin);
            }

            // set centroid
            cells[cells.size()-1].centroid = seedin;

        }
    }

    // now that all cells exist we define the neighbour relations
    for (int i = 0; i<cells.size(); i++) {

        cells[i].volume = distx * disty;

        for (int j = 0; j<cells[i].edges.size(); j++) {
            if (cells[i].edges[j].is_boundary == false) {
                if (j == 0) {
                    cells[i].edges[j].neighbour = &cells[i - (i%n_hor) + ((i+n_hor - 1)%n_hor)];
                } else if (j == 1) {
                    cells[i].edges[j].neighbour = &cells[((((i - (i%n_hor))/n_hor)+1)%n_vert)*n_hor + (i%n_hor)];
                } else if (j == 2) {
                    cells[i].edges[j].neighbour = &cells[i - (i%n_hor) + ((i + 1)%n_hor)];
                } else if (j == 3) {
                    cells[i].edges[j].neighbour = &cells[((((i - (i%n_hor))/n_hor)+n_vert - 1)%n_vert)*n_hor + (i%n_hor)];
                }
                
            }
        }
    }
}


// generates vmesh using vmp and converts it into data usable for this mesh type
template <typename CellType>
void Mesh<CellType>::generate_vmesh2D(vector<Point> pts, bool repeating, bool point_insertion) {

    // preprocessing for repeating boundary conditions
    vector<Point> points_plus_ghost;
    int initial_pts_size = pts.size();
    if (repeating) {
        points_plus_ghost.reserve(pts.size() * 9);

        // put points into (middle/middle) block by shrinking them by a factor of 3
        for (int i = 0; i<pts.size(); i++) {
            points_plus_ghost.emplace_back((pts[i].x/3.0) + 1.0/3.0, (pts[i].y/3.0) + 1.0/3.0);
        }

        // add the same shrinked points again but shifted in all other 8 third blocks (up/middle/down, left/middle/right)
        vector<double> pos_X = {0., 1., 2., 0., 2., 0., 1., 2.};
        vector<double> pos_Y = {0., 0., 0., 1., 1., 2., 2., 2.};
        for (int i = 0; i<8; i++) {
            for (int j = 0; j<pts.size(); j++) {
                points_plus_ghost.emplace_back((pts[j].x/3.0) + pos_X[i] * 1.0/3.0, (pts[j].y/3.0) + pos_Y[i] * 1.0/3.0);
            }
        }
        
        // replace pts with pts + additional ghost cells (eg 8 times the pts all around)
        pts = points_plus_ghost;
    }

    // generate vmesh
    VoronoiMesh vmesh(pts);
    vmesh.do_point_insertion();
    //vmesh.construct_mesh();
    if (point_insertion) {
        //vmesh.do_point_insertion();
    } else {
        //vmesh.construct_mesh();
    }


    // loop through all cells (of initial pts vector) to set everything but neighbour relations
    for (int i = 0; i<initial_pts_size; i++) {

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

        // push back new cell in cells vector
        if constexpr (is_same_v<CellType, DG_Q_Cell> == true) {
            cells.emplace_back(seedin, edgesin, N_basisfunc);
        } else {
            cells.emplace_back(seedin, edgesin);
        }

        // set centroid
        cells[i].centroid = vmesh.vcells[i].get_centroid();

    }

    // loop through all cells (of initial pts size) and set neighbour relations
    for (int i = 0; i<initial_pts_size; i++) {

        cells[i].volume = vmesh.vcells[i].get_area();

        for (int j = 0; j<vmesh.vcells[i].edges.size(); j++) {

            // neighbour index as used in vmesh, using modulo here to get neighbour relations
            // correct even for repeating boundaries since index will be multiple of original index
            int neigbour_index;
            neigbour_index = (vmesh.vcells[i].edges[j].index2)%initial_pts_size;

            // exclude boundaries
            if (neigbour_index >= 0) {

                // neighbour as defined before over index now with adress
                cells[i].edges[j].neighbour = &cells[neigbour_index];

            }

        }

    }

    // if repeating boundary conditions rescale the mesh back to normal
    if (repeating) {

        // loop through all cells
        for (int i = 0; i<cells.size(); i++) {

            // redo scaling for seed and volume
            cells[i].seed.x = (cells[i].seed.x - 1.0/3.0) * 3.0;
            cells[i].seed.y = (cells[i].seed.y - 1.0/3.0) * 3.0;
            cells[i].centroid.x = (cells[i].centroid.x - 1.0/3.0) * 3.0;
            cells[i].centroid.y = (cells[i].centroid.y - 1.0/3.0) * 3.0;
            cells[i].volume = cells[i].volume * 3.0 * 3.0;

            // loop through all faces
            for (int j = 0; j<cells[i].edges.size(); j++) {

                face f = cells[i].edges[j];

                // redo scaling for face positions
                f.a.x = (f.a.x - 1.0/3.0) * 3.0;
                f.a.y = (f.a.y - 1.0/3.0) * 3.0;
                f.b.x = (f.b.x - 1.0/3.0) * 3.0;
                f.b.y = (f.b.y - 1.0/3.0) * 3.0;
                cells[i].edges[j].a.x = f.a.x;
                cells[i].edges[j].a.y = f.a.y;
                cells[i].edges[j].b.x = f.b.x;
                cells[i].edges[j].b.y = f.b.y;

                // recalculate length with correct scaling
                cells[i].edges[j].length = sqrt((f.a.x - f.b.x)*(f.a.x - f.b.x) + (f.a.y - f.b.y)*(f.a.y - f.b.y));

            }

        }

    }


}


// does the lloyd_iterations to some points
template <typename CellType>
void Mesh<CellType>::do_lloyd_iterations(vector<Point>* pts, int lloyd_iterations) {
    // preprocessing step to change pts for mesh into pts for approx centroidal vmesh
    if (lloyd_iterations != 0) {
        cout << "start lloyd_iterations" << endl;

        // calculate original mesh
        VoronoiMesh initial_vmesh(*pts);
        initial_vmesh.do_point_insertion();
        
        // do multiple iterations of lloyds algorithm
        for (int i = 0; i<lloyd_iterations; i++) {

            // calculate centroids
            vector<Point> centroids;
            centroids.reserve(initial_vmesh.vcells.size());
            for (int i = 0; i<initial_vmesh.vcells.size(); i++) {
                centroids.push_back(initial_vmesh.vcells[i].get_centroid());
            }

            // replace original mesh seeds with centroids and calculate mesh again
            initial_vmesh = VoronoiMesh(centroids);
            initial_vmesh.do_point_insertion();

        }

        // after iterations store final calculated seeds in pts
        vector<Point> centroidal_seeds;
        centroidal_seeds.reserve(initial_vmesh.vcells.size());
        for (int i = 0; i<initial_vmesh.vcells.size(); i++) {
            centroidal_seeds.push_back(initial_vmesh.vcells[i].seed);
        }
        *pts = centroidal_seeds;
        cout << "finished lloyd_iterations" << endl;
    }
}


// generates a 1D uniform grid with all the neighbour relations and so on
template <typename CellType>
void Mesh<CellType>::generate_uniform_grid1D(Point start, int n, double dist, bool repeating) {
    generate_uniform_grid2D(start, n, 1, dist, 1, repeating);
}


// generates a 1D voronoi mesh, only works if points are between 0 and 1
template <typename CellType>
void Mesh<CellType>::generate_vmesh1D(vector<Point> pts, bool repeating) {

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

        if (repeating == true) {
            f0.is_boundary = false;
            f2.is_boundary = false;
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
        if constexpr (is_same_v<CellType, DG_Q_Cell> == true) {
            cells.emplace_back(seedin, edgesin, N_basisfunc);
        } else {
            cells.emplace_back(seedin, edgesin);
        }

        cells[cells.size()-1].centroid = seedin;

    }

    // now that all cells exist we define the neighbour relations
    for (int i = 0; i < cells.size(); i++) {

        cells[i].volume = cells[i].edges[0].length * cells[i].edges[1].length;

        if (cells[i].edges[0].is_boundary == false) {
            cells[i].edges[0].neighbour = &cells[i - (i%cells.size()) + ((i+cells.size() - 1)%cells.size())];
        }
        if (cells[i].edges[2].is_boundary == false) {
            cells[i].edges[2].neighbour = &cells[i - (i%cells.size()) + ((i + 1)%cells.size())];
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

    // set Q for all indices between a and b to value
    for (int i = a; i < b; i+=step) {
        cells[i].Q = value;
    }

}


// sets the initial condition according to given analytical Q_circle
template <typename CellType>
void Mesh<CellType>::initialize_Q_circle(Point p0, double r, double Qval) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Q_Cell> == false) {
        cerr << "initalize_Q_circle called with wrong cell type, you must use Q_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    // set Q_values for initial circle using analytical solution at t = 0
    for (int i = 0; i < cells.size(); i++) {
        cells[i].Q = Qval * advecting_circle(cells[i].seed, 0, Point(0, 0), p0, r);
    }

}


// Function to initalize a gaussian for shallow water equations
template <typename CellType>
void Mesh<CellType>::initialize_SWE_gaussian(Point p0, double A, double sigma) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, SWE_Cell> == false) {
        cerr << "initialize_SWE_gaussian called with wrong cell type, you must use SWE_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < cells.size(); i++) {

        // set height according to gaussian shape
        double dx = cells[i].seed.x - p0.x;
        double dy = cells[i].seed.y - p0.y;
        cells[i].h += A * exp(- (dx * dx + dy * dy) / (2 * sigma * sigma));

        // no inital velocities
        cells[i].u = 0;
        cells[i].v = 0;

    }

}


// sets initial conditions for SWE dam break (x, diagonal, circular)
template <typename CellType>
void Mesh<CellType>::initialize_SWE_dam_break(double h1, double h2, double pos, int dam_break_type) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, SWE_Cell> == false) {
        cerr << "initialize_SWE_dam_break called with wrong cell type, you must use SWE_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    // go through all cells
    for (int i = 0; i< cells.size(); i++) {

        bool set_values = false;

        // dam break in x direction
        if (dam_break_type == 0) {
            set_values = (cells[i].seed.x <= pos);
        // dam break in y direction
        } else if (dam_break_type == 1) {
            set_values = (cells[i].seed.y <= pos);
        // dam break in diagonal direction
        } else if (dam_break_type == 2) {
            set_values = (cells[i].seed.y + cells[i].seed.x  <= 2*pos);
        // circular dam break
        } else if (dam_break_type == 3) {
            set_values = (sqrt(cells[i].seed.y*cells[i].seed.y + cells[i].seed.x*cells[i].seed.x)  <= pos);
        }

        // set values according to bool
        if (set_values) {
            cells[i].h += h1-1;
            cells[i].u = 0;
            cells[i].v = 0;
        } else {
            cells[i].h += h2-1;
            cells[i].u = 0;
            cells[i].v = 0;
        }
    }
    
}


// function to initalize sods shock tube on the mesh
template <typename CellType>
void Mesh<CellType>::initialize_euler_shock_tube() {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Euler_Cell> == false) {
        cerr << "initialize_euler_shock_tube called with wrong cell type, you must use Euler_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        if (cells[i].seed.x + cells[i].seed.y < 1) {
            cells[i].rho += 1-1;
            cells[i].u = 0;
            cells[i].v = 0;
            cells[i].E += (1/(cells[i].gamma - 1)) * 1 - 1;
        } else {
            cells[i].rho += 0.125-1;
            cells[i].u = 0;
            cells[i].v = 0;
            cells[i].E += (1/(cells[i].gamma - 1)) * 0.1 - 1;
        }
    }
}

// function to initalize a kelvin helmholtz instability on the mesh
template <typename CellType>
void Mesh<CellType>::initialize_kelvin_helmholtz() {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Euler_Cell> == false) {
        cerr << "initalize_kelvin_helmholtz called with wrong cell type, you must use Euler_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        double pi = 3.14159265358979323846;

        if (cells[i].seed.y > 0.3 && cells[i].seed.y < 0.7) {
            cells[i].rho = 0.5;
            cells[i].u = 0.3;
            cells[i].v = 0;
            cells[i].E = (1/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u);
        } else {
            cells[i].rho = 0.2;
            cells[i].u = -0.3;
            cells[i].v = 0;
            cells[i].E = (1/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u);
        }
    }
  
}



// function to initalize a rayleigh taylor instability on the mesh
template <typename CellType>
void Mesh<CellType>::initialize_rayleigh_taylor(Point g) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Euler_Cell> == false) {
        cerr << "initalize_rayleigh_taylor called with wrong cell type, you must use Euler_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        double pi = 3.14159265358979323846;



        if (cells[i].seed.y > 0.5 + 0.03*cos(pi*2*cells[i].seed.x*2)) {
            cells[i].rho = 1;
            cells[i].u = 0;
            cells[i].v = 0;
            double P = 1 + cells[i].rho * g.y * (cells[i].seed.y - 0.5);
            cells[i].E = (P/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u);
        } else {
            cells[i].rho = 0.2;
            cells[i].u = 0;
            cells[i].v = 0;
            double P = 1 + cells[i].rho * g.y * (cells[i].seed.y - 0.5);
            cells[i].E = (P/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u);
        }
    }
  
}


// function to initialize a constant flow
template <typename CellType>
void Mesh<CellType>::initialize_const_flow(Point v) {
    for (int i = 0; i<cells.size(); i++) {
        cells[i].u += v.x;
        cells[i].v += v.y;
    }
}


// function to initialize a quad shock
template <typename CellType>
void Mesh<CellType>::initialize_quad_shock() {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Euler_Cell> == false) {
        cerr << "initialize_quad_shock called with wrong cell type, you must use Euler_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        if (cells[i].seed.x >= 0.5 && cells[i].seed.y >= 0.5) {
            cells[i].rho = 1.5;
            cells[i].u = 0;
            cells[i].v = 0;
            cells[i].E = (1.5/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        } else if (cells[i].seed.x < 0.5 && cells[i].seed.y > 0.5) {
            cells[i].rho = 0.5323;
            cells[i].u = 1.206;
            cells[i].v = 0;
            cells[i].E = (0.3/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        } else if (cells[i].seed.x < 0.5 && cells[i].seed.y < 0.5) {
            cells[i].rho = 0.138;
            cells[i].u = 1.206;
            cells[i].v = 1.206;
            cells[i].E = (0.029/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        } else if (cells[i].seed.x > 0.5 && cells[i].seed.y < 0.5) {
            cells[i].rho = 0.5323;
            cells[i].u = 0;
            cells[i].v = 1.206;
            cells[i].E = (0.3/(cells[i].gamma - 1)) + 0.5*cells[i].rho*(cells[i].u*cells[i].u + cells[i].v*cells[i].v);
        }
    }

}






// Function to create an internal boundary structure (just a square at the moment)
template <typename CellType>
void Mesh<CellType>::initialize_boundary_struct(Point p0, double l_x, double l_y) {

    // loop through all cells
    for (int i = 0; i < cells.size(); i++) {

        // if cell[i] is inside square make it boundary cell (in principle any other if condition could be built here)
        if (cells[i].seed.x < p0.x + l_x && cells[i].seed.x > p0.x && cells[i].seed.y < p0.y + l_y && cells[i].seed.y > p0.y) { 
            make_cell_boundary_cell(i);
        }
    }
}


// makes cell[i] boundary cell
template <typename CellType>
void Mesh<CellType>::make_cell_boundary_cell(int i) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Q_Cell> == true) {
        // boundary cells have Q = -INFINITY such that they will not be plotted in visualization
        cells[i].Q = -INFINITY;
    }

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, SWE_Cell> == true) {
        // boundary cells have h = -INFINITY such that they will not be plotted in visualization
        cells[i].h = -10;
        cells[i].u = 0;
        cells[i].v = 0;
    }

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Euler_Cell> == true) {
        // boundary cells have h = -INFINITY such that they will not be plotted in visualization
        cells[i].rho = -INFINITY;
        cells[i].u = 0;
        cells[i].v = 0;
        cells[i].E = -INFINITY;
    }

    // loop through edges of that boundary cell
    for (int j = 0; j < cells[i].edges.size(); j++) {

        // if already boundary there is nothing to change
        if (cells[i].edges[j].is_boundary == false) {

            // set face from internal side to boundary
            cells[i].edges[j].is_boundary = true;

            // go through all faces of the edges[j].neighbour were looking at
            for (int k = 0; k < cells[i].edges[j].neighbour->edges.size(); k++) {

                // find face with neighbour == our cell, then set its boundary also true
                if (cells[i].edges[j].neighbour->edges[k].neighbour == &cells[i]) {

                    // set that boundary true
                    cells[i].edges[j].neighbour->edges[k].is_boundary = true;
                }
            }
        }
    }
}


// loads structure and places it in meshpoints
// structures are given by a list of points specifying the verticies
// verticies need to be in clockwise orientation!!!
template <typename CellType>
int Mesh<CellType>::add_struct(vector<Point>* pts, double dist_a, double safety_dist, string structname) {
    
    // load structpoints from file
    vector<Point> verticies;
    ifstream file("../src/files/" + structname + ".csv");
    string line;

    while (getline(file, line)) {
        stringstream ss(line);
        string x_string, y_string;

        if (getline(ss, x_string, ',') && getline(ss, y_string, ',')) {
            double x = stod(x_string);
            double y = stod(y_string);
            verticies.emplace_back(x, y);
        }
    }
    file.close();

    // calculate corresponding seeds in inner_seeds/outer_seeds
    vector<Point> edge_vectors;
    vector<Point> a = verticies;
    vector<Point> b;
    vector<Point> midpoints;
    for (int i = 0; i < verticies.size(); i++) {
        b.emplace_back(verticies[(i+1)%verticies.size()]);
    }
    for (int i = 0; i < verticies.size(); i++) {
        edge_vectors.emplace_back((b[i].x - a[i].x), (b[i].y - a[i].y));
        midpoints.emplace_back((a[i].x + b[i].x)/2.0, (a[i].y + b[i].y)/2.0);
    }
    vector<Point> normal_vectors;
    for (int i = 0; i<edge_vectors.size(); i++) {
        normal_vectors.emplace_back((-1*(edge_vectors[i].y)/(sqrt(edge_vectors[i].x*edge_vectors[i].x + edge_vectors[i].y*edge_vectors[i].y))),
                                    ((edge_vectors[i].x)/(sqrt(edge_vectors[i].x*edge_vectors[i].x + edge_vectors[i].y*edge_vectors[i].y))));
    }
    vector<Point> inner_points;
    vector<Point> outer_points;
    for (int i = 0; i< normal_vectors.size(); i++) {
        inner_points.emplace_back((midpoints[i].x - dist_a * normal_vectors[i].x), (midpoints[i].y - dist_a * normal_vectors[i].y));
        outer_points.emplace_back((midpoints[i].x + dist_a * normal_vectors[i].x), (midpoints[i].y + dist_a * normal_vectors[i].y));
    }


    // remove seeds inside of struct (using ray-casting-algorithm, ray along x direction) and seeds to close to struct
    vector<Point> pts_removed;
    for (int i = 0; i < (*pts).size(); i++) {
        Point pt = (*pts)[i];
        int counter = 0;
        double min_dist_to_outer_seeds = 2;
        for (int j = 0; j < verticies.size(); j++) {
            double crossing_x = a[j].x + ((pt.y - a[j].y)/(b[j].y - a[j].y)) * (b[j].x - a[j].x);
            //cout << i << " _ " << (pt.y > a[j].y && pt.y < b[j].y) << " " << (pt.y < a[j].y && pt.y > b[j].y) << endl;
            if (((pt.y > a[j].y && pt.y < b[j].y) || (pt.y < a[j].y && pt.y > b[j].y)) && (pt.x < crossing_x)) {
                counter += 1;
            }

            double dist_to_outer_seed = sqrt((pt.x - outer_points[j].x)*(pt.x - outer_points[j].x) + (pt.y - outer_points[j].y)*(pt.y - outer_points[j].y));
            if (dist_to_outer_seed < min_dist_to_outer_seeds) {
                min_dist_to_outer_seeds = dist_to_outer_seed;
            }
        }

        if (counter%2 == 0 && min_dist_to_outer_seeds > safety_dist) {
            pts_removed.push_back((*pts)[i]);
            //(*pts).erase((*pts).begin() + i);
        }
    }
    *pts = pts_removed;



    // return int = inner_seeds.size(), set *pts =  {inner_seeds, outer_seeds, pts}
    vector<Point> new_pts;
    for (int i = 0; i<inner_points.size(); i++) {
        new_pts.push_back(inner_points[i]);
    }
    for (int i = 0; i<outer_points.size(); i++) {
        new_pts.push_back(outer_points[i]);
    }
    for (int i = 0; i<(*pts).size(); i++) {
        new_pts.push_back((*pts)[i]);
    }
    *pts = new_pts;

    return inner_points.size();


}



// SAVE DATA TO FILE ------------------------------------------------------------------------------
// saves mesh into csv file readable for python script
template <typename CellType>
void Mesh<CellType>::save_mesh(int file_nr, string name, double dt) {

    // open correct file
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
    output_file << "seed.x, seed.y | a.x, a.y ; a.x ... ; | Quantities" << endl;

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

        // store sim time
        output_file << file_nr * dt << ",";

        // if cell type is Q_cell or Conway Cell save Q
        if constexpr (is_same_v<CellType, Q_Cell> || is_same_v<CellType, Conway_Cell>) {
            output_file << cells[i].Q;
        }

        // if cell type is SWE_cell save h, u, v
        if constexpr (is_same_v<CellType, SWE_Cell>) {
            output_file << cells[i].h << "," << cells[i].u << "," << cells[i].v;
        }

        // if cell type is Euler_cell save rho, u, v, E, P
        if constexpr(is_same_v<CellType, Euler_Cell>) {
            output_file << cells[i].rho << "," << cells[i].u << "," << cells[i].v << "," << cells[i].E << "," << cells[i].get_P();
        }

        // if cell type is DG_Q_Cell save all Q out of vector
        if constexpr(is_same_v<CellType, DG_Q_Cell>) {
            for (int l = 0; l < N_basisfunc - 1; l++) {
                output_file << cells[i].Q(l) << ",";
            }
            output_file << cells[i].Q(N_basisfunc - 1);
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

    // get Q_values for numerical and analytical solution at current timestep
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

    // write caluclated L1 difference for given time in file
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

    // get Q_values for numerical and analytical solution at current timestep
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

    // write caluclated L1 difference for given time in file
    output_file << t << "," << L1_error(Qs_num, Qs_ana) << endl;
    output_file.close();
}

// function to calculate and store the L1 error of an advecting circle, requires Q_cells
template <typename CellType>
void Mesh<CellType>::save_L1_swe_dam_break(double t, bool reset_file) {

    vector<double> h_num;
    vector<double> h_ana;
    h_num.reserve(cells.size());
    h_ana.reserve(cells.size());

    // get h_values for numerical and analytical solution at current timestep
    for (int i = 0; i<cells.size(); i++) {
        h_num.push_back(cells[i].h);
        h_ana.push_back(swe1D_dam_break(cells[i].seed, t));
    }

    // open new file or in append mode
    ofstream output_file;
    if (reset_file) {
        output_file = ofstream("../src/files/L1_error.csv");
    } else {
        output_file = ofstream("../src/files/L1_error.csv", ios::app);
    }

    // write caluclated L1 difference for given time in file
    output_file << t << "," << L1_error(h_num, h_ana) << endl;
    output_file.close();
}




template <typename CellType>
Mesh<CellType>::~Mesh() {}