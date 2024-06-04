#include <iostream>
#include <vector>
#include <string>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "Mesh.h"

#ifndef Solver_h
#define Solver_h

template <typename CellType>
class Solver {

public:
    Solver<CellType>();
    Solver<CellType>(Mesh<CellType>* gridin);
    ~Solver<CellType>();

    // equation solvers
    void diffusion_like(double dt);
    void conway();
    void advection(double dt, Point v);
    void shallow_water_1D_cartesian(double dt);
    
private:

    // grid to solve the equations on
    Mesh<CellType>* grid;

    // helper functions
    Point get_normal_vec(Point a, Point b);

    // sub solvers
    void advection_cartesian(double dt, Point v);
    void advection_vmesh(double dt, Point v);


};

#include "Solver.tpp"

#endif