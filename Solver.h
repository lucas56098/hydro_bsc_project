#include <iostream>
#include <vector>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "Mesh.h"

#ifndef Solver_h
#define Solver_h

class Solver {

public:
    Solver();
    Solver(Mesh* gridin);
    ~Solver();

    Mesh* grid;

    void diffusion_like_step(double dt);
    void conway_step();
    void advection_finite_difference_upwind_cartesian1D(double dt, double v);
};

#endif