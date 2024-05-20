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

    Mesh<CellType>* grid;

    Point get_normal_vec(Point a, Point b);
    void diffusion_like_step(double dt);
    void conway_step();
    void advection_finite_difference_upwind_cartesian(double dt, Point v);
    void advection_vmesh(double dt, Point v);
};

#include "Solver.tpp"

#endif