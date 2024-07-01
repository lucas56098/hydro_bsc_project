#include <vector>
#include <string>
#include "cell_types/Cell.h"
#include "cell_types/Point.h"
#include "cell_types/Q_Cell.h"
#include "cell_types/Conway_Cell.h"
#include "cell_types/SWE_Cell.h"
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
    double n_horizontal;
    double n_vertical;
    
    // generate mesh
    void generate_grid(bool cartesian, bool is_1D, int N_row, int lloyd_iterations = 0, bool repeating = false);
    
    // initial conditions
    void initialize_Q_cells(int a, int b, double value = 1, int step = 1);
    void initalize_Q_circle(Point p0 = Point(0.5, 0.5), double r = 0.1, double Qval = 1);
    void inialize_boundary_struct(Point p0, double l_x, double l_y);
    void initalize_SWE_gaussian(Point p0 = Point(0.5, 0.5), double A = 1, double sigma = 0.5);
    void initalize_SWE_dam_break(double h1 = 2.0, double h2 = 1.0, double pos = 0.5, int dam_break_type = 0);

    // save mesh
    void save_mesh(int file_nr, string name, double dt);
    void save_Q_diff(double t, bool reset_file = false, bool is_density = true);
    void save_L1_adv_circle(double t, bool reset_file, Point v = Point(0.5, 0.3), Point p0 = Point(0.5, 0.5), double r = 0.1);
    void save_L1_adv_1Dstepfunc(double t, bool reset_file, double v = 0.5, double a0 = 0, double b0 = 0.1);
    void save_L1_swe_dam_break(double t, bool reset_file);


private:

    // helper functions
    vector<Point> generate_seed_points(int N, bool fixed_random_seed, double min, int max, int rd_seed, bool sort_pts, int sort_precision, int sort_scheme);
    int get_sort_index(Point pt, int sort_grid_size, int sort_scheme);
    void make_cell_boundary_cell(int i);

    // different grid generations
    void generate_uniform_grid2D(Point start, int n_hor, int n_vert, double distx, double disty, bool repeating = false);
    void generate_vmesh2D(vector<Point> pts, int lloyd_iterations, bool repeating = false);
    void generate_uniform_grid1D(Point start, int n, double dist, bool repeating = false);
    void generate_vmesh1D(vector<Point> pts, bool repeating = false);
};

#include "Mesh.tpp"

#endif