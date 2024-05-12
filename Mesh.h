#include <vector>
#include "cell_types/Cell.h"
#include "cell_types/Point.h"
#include "cell_types/Q_Cell.h"
#include "cell_types/Conway_Cell.h"
#include "vmp/VoronoiMesh.h"

#ifndef Mesh_h
#define Mesh_h

class Mesh {

public:
    Mesh();
    Mesh(int cell_typein);
    ~Mesh();

    // internal lists and cell type
    int cell_type;  // 0 = Cell, 1 = Q_Cell, 2 = Conway_Cell
    vector<Cell> cells;
    vector<Q_Cell> q_cells;
    vector<Conway_Cell> conway_cells;
    
    // function for adaptive cell handling
    void smart_push_back(Point seedin, vector<face> edgesin);
    face& smart_cells_edge(int i, int j);
    Cell& smart_cells(int i);
    int get_edges_size(int i);
    int get_cells_size();

    // generate mesh
    void generate_uniform_grid2D(Point start, int n_hor, int n_vert, double distx, double disty);
    void generate_vmesh2D(vector<Point> pts);
    void generate_uniform_grid1D(Point start, int n, double dist);
    void generate_vmesh1D(vector<Point> pts);
    
    // initial conditions
    void initialize_cells(int a, int b, double value = 1);
    void initialize_random();

    // save mesh
    void save_mesh(int file_nr, bool cartesian = true);
    void save_Q_diff(double initial_total_Q, bool reset_file = false);

};

#endif