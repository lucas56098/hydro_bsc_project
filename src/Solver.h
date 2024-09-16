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
    void euler(double dt, int boundary_cond = -1, int sim_order = 2, Point g = Point(0, 0));

    // other
    Point calc_struct_force(int a = 200, int b = 400);

    

private:

    // grid to solve the equations on
    Mesh<CellType>* grid;

    // helper functions
    Point get_normal_vec(Point a, Point b);
    Point get_f_mid(int i, int j);

        //swe related
        array<double, 3> get_huv_j(array<double, 3> huv_i, int i, int j, int boundary_cond);
        array<double, 3> get_flux_f_swe(array<double, 3> huv, double g);
        array<double, 3> huv_to_U(array<double, 3> huv);
        array<double, 3> U_to_huv(array<double, 3> U);
        array<double, 3> rotate2Dswe(array<double, 3> vec, double angle);

        // euler related
        array<double, 4> get_puvE_j(array<double, 4> puvE_i, int i, int j, int boundary_cond);
        array<double, 4> get_flux_f_euler(array<double, 4> puvE);
        array<double, 4> rotate2Deuler(array<double, 4> vec, double angle);
        double get_P_ideal_gas(array<double, 4> puvE);
        array<double, 4> puvE_to_U(array<double, 4> puvE);
        array<double, 4> U_to_puvE(array<double, 4> U);
    

    // sub solvers
    void advection_cartesian(double dt, Point v);
    void advection_vmesh(double dt, Point v);

    // riemann solvers
    //array<double, 3> roe_solver_swe_2D(array<double, 3> huv_i_n, array<double, 3> huv_j_n, array<double, 3> f_i, array<double, 3> f_j, array<double, 3> g_i, array<double, 3> g_j, double g, Point n, double add_diffusion_coeff = 0.0);
    array<double, 3> hll_solver_swe_2D(array<double, 3> huv_i_n, array<double, 3> huv_j_n, double g, int i, int j, double add_diffusion_coeff);
    array<double, 4> hll_solver_euler_2D(array<double, 4> puvE_i_n, array<double, 4> puvE_j_n, int i, int j);

    // second order swe
    array<array<double, 3>, 2> calc_swe_gradients(int i, int boundary_cond, double g);
    array<double, 3> linear_extrapolate_swe(array<double, 3> huv_i_n, vector<array<array<double, 3>, 2>> &gradients, double dt, int i, int j, double g);
    array<double, 3> linear_extrapolate_neighbour_swe(int i, int j, array<double, 3> huv_j_n, array<double, 3> huv_i_ext, vector<array<array<double, 3>, 2>> &gradients, double dt, double g, int boundary_cond);

    // second order euler fv
    array<array<double, 4>, 2> calc_euler_gradients(int i, int boundary_cond);
    array<array<double, 4>, 2> slope_limit_maxmin(array<array<double, 4>, 2> gradientU, int i, array<double, 4> U_i, int boundary_cond);
    array<array<double, 4>, 2> slope_limit_tvd(array<array<double, 4>, 2> gradientU, int i, array<double, 4> U_i, int boundary_cond, double theta = 1);
    array<double, 4> linear_extrapolate_euler(array<double, 4> puvE_i_n, vector<array<array<double, 4>, 2>> &gradients, double dt, int i, int j);
    array<double, 4> linear_extrapolate_neighbour_euler(int i, int j, array<double, 4> puvE_j_n, array<double, 4> puvE_i_ext, vector<array<array<double, 4>, 2>> &gradients, double dt, int boundary_cond);


};

#include "Solver.tpp"

#endif