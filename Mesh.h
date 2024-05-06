#include <vector>
#include "Cell.h"
#include "Point.h"
#include "vmp/VoronoiMesh.h"

#ifndef Mesh_h
#define Mesh_h

class Mesh {

public:
    Mesh();
    ~Mesh();

    void generate_uniform_grid(Point start, int n_hor, int n_vert, double dist);
    void generate_from_vmesh(VoronoiMesh vmesh);
    void updateQ();
    void printQ(int n_hor, int n_vert);
    void save_mesh(int file_nr);

    vector<Cell> cells;


};

#endif