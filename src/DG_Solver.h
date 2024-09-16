#include <iostream>
#include <vector>
#include <array>
#include <string>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "Mesh.h"

#ifndef DG_Solver_h
#define DG_Solver_h

template <typename CellType>
class DG_Solver {

public:
    DG_Solver<CellType>();
    DG_Solver<CellType>(Mesh<CellType>* gridin);
    ~DG_Solver<CellType>(); 

    void advection1D(double dt, double slope_limited = 0, double u = 1);
    void advection2D(double dt, double slope_limited = 0, Point u = Point(0.5, 0.5), bool initialize = false);

private:

    // grid to solve the equations on
    Mesh<CellType>* grid;

    // helper functions for 1D scalar advection
    double minmod(double a, double b, double c);
    double get_element_average_adv1D(int index, vector<VectorXd>& new_Qs);
    vector<VectorXd> slope_limit_adv1D(vector<VectorXd> new_Qs, double slope_limited);

    // Matrices for 2D scalar advection
    MatrixXd F_hat_A;
    MatrixXd F_hat_B;
    MatrixXd F_hat_C;
    MatrixXd F_hat_D;

    // helper functions for cartesian 2D scalar advection 
    MatrixXd get_F_hat_A_adv2D();
    MatrixXd get_F_hat_B_adv2D();
    MatrixXd get_F_hat_C_adv2D();
    MatrixXd get_F_hat_D_adv2D();
    vector<VectorXd> slope_limit_adv2D(vector<VectorXd> new_Qs, double slope_limited);


};

#include "DG_Solver.tpp"

#endif