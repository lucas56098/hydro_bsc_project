#include <vector>
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

    // internal lists and cell type
    vector<CellType> cells;
    double total_Q_initial;
    
    // generate mesh
    void generate_uniform_grid2D(Point start, int n_hor, int n_vert, double distx, double disty);
    void generate_vmesh2D(vector<Point> pts);
    void generate_uniform_grid1D(Point start, int n, double dist);
    void generate_vmesh1D(vector<Point> pts);
    
    // initial conditions
    void initialize_Q_cells(int a, int b, double value = 1, int step = 1);
    void initialize_random();

    // save mesh
    void save_Q_mesh(int file_nr, bool cartesian = true);
    void save_Q_diff(bool reset_file = false, bool is_density = true);

};

#include "Mesh.tpp"

#endif