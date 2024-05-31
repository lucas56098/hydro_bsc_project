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
    void generate_grid(bool cartesian, bool is_1D, int N_row);
    
    // initial conditions
    void initialize_Q_cells(int a, int b, double value = 1, int step = 1);

    // save mesh
    void save_mesh(int file_nr, string name);
    void save_Q_diff(bool reset_file = false, bool is_density = true);

private:

    // helper functions for vmesh generation out of Voronoi_Mesh_Project (vmp)
    vector<Point> generate_seed_points(int N, bool fixed_random_seed, double min, int max, int rd_seed, bool sort_pts, int sort_precision, int sort_scheme);
    int get_sort_index(Point pt, int sort_grid_size, int sort_scheme);

    // different grid generations
    void generate_uniform_grid2D(Point start, int n_hor, int n_vert, double distx, double disty);
    void generate_vmesh2D(vector<Point> pts);
    void generate_uniform_grid1D(Point start, int n, double dist);
    void generate_vmesh1D(vector<Point> pts);
};

#include "Mesh.tpp"

#endif