#include <vector>
#include <string>
#include "cell_types/Cell.h"
#include "cell_types/Point.h"
#include "cell_types/Q_Cell.h"
#include "cell_types/Conway_Cell.h"
#include "vmp/VoronoiMesh.h"

#ifndef Mesh_h
#define Mesh_h

template <typename CellType>
class Mesh {

public:
    Mesh<CellType>();
    ~Mesh<CellType>();

    // cell list
    vector<CellType> cells;
    double total_Q_initial;
    bool is_cartesian = false;
    
    // generate mesh
    void generate_grid(bool cartesian, bool is_1D, int N_row, int lloyd_iterations = 0, bool repeating = false);
    
    // initial conditions
    void initialize_Q_cells(int a, int b, double value = 1, int step = 1);
    void initalize_Q_circle(Point p0 = Point(0.5, 0.5), double r = 0.1);

    // save mesh
    void save_mesh(int file_nr, string name, double dt);
    void save_Q_diff(double t, bool reset_file = false, bool is_density = true);
    void save_L1_adv_circle(double t, bool reset_file, Point v = Point(0.5, 0.3), Point p0 = Point(0.5, 0.5), double r = 0.1);
    void save_L1_adv_1Dstepfunc(double t, bool reset_file, double v = 0.5, double a0 = 0, double b0 = 0.1);


private:

    // helper functions for vmesh generation out of Voronoi_Mesh_Project (vmp)
    vector<Point> generate_seed_points(int N, bool fixed_random_seed, double min, int max, int rd_seed, bool sort_pts, int sort_precision, int sort_scheme);
    int get_sort_index(Point pt, int sort_grid_size, int sort_scheme);

    // different grid generations
    void generate_uniform_grid2D(Point start, int n_hor, int n_vert, double distx, double disty, bool repeating = false);
    void generate_vmesh2D(vector<Point> pts, int lloyd_iterations);
    void generate_uniform_grid1D(Point start, int n, double dist, bool repeating = false);
    void generate_vmesh1D(vector<Point> pts, bool repeating = false);
};

#include "Mesh.tpp"

#endif