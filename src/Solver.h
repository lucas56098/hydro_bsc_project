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
    void shallow_water(double dt, int boundary_cond = -1, double numerical_diffusion_coeff = 0, int sim_order = 1, double g = 1.0);
    
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
    //array<double, 3> hll_solver_swe_2D(array<double, 3> huv_i_n, array<double, 3> huv_j_n, array<double, 3> f_i, array<double, 3> f_j, array<double, 3> g_i, array<double, 3> g_j, double g, Point n, double add_diffusion_coeff = 0.0);
    array<double, 3> hll_solver_swe_2D(array<double, 3> huv_i_n, array<double, 3> huv_j_n, double g, Point n, double add_diffusion_coeff);

    array<double, 3> rotate2Dswe(array<double, 3> vec, double angle);
    array<array<double, 3>, 2> calc_swe_gradients(int i, int boundary_cond, double g);
    array<double, 3> linear_extrapolate_swe(array<double, 3> huv_i_n, array<array<double, 3>, 2> gradient, double dt, int i, int j, double g);

};

#include "Solver.tpp"

#endif