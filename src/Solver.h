#include <iostream>
#include <vector>
#include <array>
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
    //void shallow_water_1D_cartesian(double dt);
    //void shallow_water_2D_cartesian(double dt, int boundary_condition = -1); // boundary condition -1 -> reflective, 1-> sink, for repeating boundary change mesh structure
    void shallow_water_2D(double dt, int boundary_cond = -1);
    
private:

    // grid to solve the equations on
    Mesh<CellType>* grid;

    // helper functions
    Point get_normal_vec(Point a, Point b);

    // sub solvers
    void advection_cartesian(double dt, Point v);
    void advection_vmesh(double dt, Point v);

    // riemann solvers
    array<double, 3> roe_solver_swe_2D(array<double, 3> huv_i_n, array<double, 3> huv_j_n, array<double, 3> f_i, array<double, 3> f_j, array<double, 3> g_i, array<double, 3> g_j, double g, Point n, double add_diffusion_coeff = 0.0);
    vector<double> hll_solver_swe_2D(vector<double> huv_i_n, vector<double> huv_j_n, vector<double> f_i, vector<double> f_j, vector<double> g_i, vector<double> g_j, double g, Point n);


};

#include "Solver.tpp"

#endif