#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "Mesh.h"
#include "DG_Solver.h"

template <typename CellType>
DG_Solver<CellType>::DG_Solver() {};

template <typename CellType>
DG_Solver<CellType>::DG_Solver(Mesh<CellType>* gridin) {
    grid = gridin;
}    


// - - - - - - - 1D scalar upwind advection - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// function to calculate next step in scalar upwind 1D advection discontinous galerkin
template <typename CellType>
void DG_Solver<CellType>::advection1D(double dt, double slope_limited, double u) {

    // vector to store new Q values
    vector<VectorXd> new_Qs;
    new_Qs.reserve(grid->cells.size());

    // calculate new Q values
    // loop through all elements T_k for that
    for (int k = 0; k < grid->cells.size(); k++) {
        
        DG_Q_Cell T_k = grid->cells[k];

        // calculate M^k_ij inverse
        T_k.M_ij_k_inv = T_k.calc_M_ij_k_inv_adv1D();

        // calculate S^k_ij
        T_k.S_ij_k = T_k.calc_S_ij_k_adv1D();

        // calculate RHS_I
        VectorXd RHS_I = u * T_k.S_ij_k * T_k.Q;

        // calculate RHS_II
        VectorXd RHS_II;
    
        // total RHS_II is given by: u * (- phi^k_j(x^k_r) phi^k_i(x^k_r)) * Q^k_j  + u * (phi^k-1_j(x^k-1_r) phi^k_i(x^k_l)) * Q^k-1_j     with sum convention over j used here

        // this is the (- phi^k_j(x^k_r) phi^k_i(x^k_r)) matrix (F_hat_A)_ij
        MatrixXd F_hat_A;
        F_hat_A = - F_hat_A.Ones(grid->N_basisfunc,grid->N_basisfunc);

        // this is the (phi^k-1_j(x^k-1_r) phi^k_i(x^k_l)) matrix (F_hat_B)_ij
        MatrixXd F_hat_B(grid->N_basisfunc,grid->N_basisfunc);
        for (int i = 0; i < grid->N_basisfunc; i++) {
            for (int j = 0; j < grid->N_basisfunc; j++) {
                F_hat_B(i, j) = (i%2 == 0) ? 1 : -1;
            }
        }

        VectorXd Q_k = T_k.get_coeff();
        if (T_k.edges[0].is_boundary == false) {

            VectorXd Q_km1 = T_k.edges[0].neighbour->get_coeff();
            RHS_II = u * F_hat_A * Q_k + u * F_hat_B * Q_km1;

        } else {

            // inflow BC with Q_km1 = {0, ...}
            RHS_II = u * F_hat_A * Q_k;

        }

        // assemble RHS
        VectorXd RHS = RHS_I + RHS_II;

        // calculate new Q using forward euler
        VectorXd Q_k_np1 = Q_k + dt * T_k.M_ij_k_inv * RHS; 

        // store it in new Q vector
        new_Qs.push_back(Q_k_np1);

    }

    // optional slope limiting as a post processing step
    if (grid->N_basisfunc > 1 && slope_limited != 0) {
        /*
        vector<VectorXd> sl_new_Qs;
        sl_new_Qs.reserve(grid->cells.size());
        sl_new_Qs = new_Qs;

        // go through all elements
        for (int i = 0; i < grid->cells.size(); i++) {

            // get averaged slope on T_k (for N_p = 1 the slope is constant so already averaged)
            double dx = grid->cells[i].edges[1].b.x - grid->cells[i].edges[1].a.x;
            double dy = 0;
            for (int a = 0; a < grid->N_basisfunc; a++) {
                if (a%2 != 0) {dy += 2*new_Qs[i](a);}
            }
            double avg_slope_k = dy/dx;


            // get averaged q values for T_k-1, T_k and T_k+1
            double q_avg_m = (i >= 1) ? get_element_average_adv1D(i-1, new_Qs) : 0;
            double q_avg_k = get_element_average_adv1D(i, new_Qs);
            double q_avg_p = (i <= new_Qs.size()-2) ? get_element_average_adv1D(i+1, new_Qs) : get_element_average_adv1D(i, new_Qs);    // note that we do respect the inflow and outflow BC here


            // calc averaged slopes between T_k and T_k-1 as well as T_k+1 and T_k
            double avg_slope_m = slope_limited*(q_avg_k - q_avg_m)/dx;
            double avg_slope_p = slope_limited*(q_avg_p - q_avg_k)/dx;

            // do minmod limiting of averaged slopes
            double slope_minmod = minmod(avg_slope_k, avg_slope_m, avg_slope_p);

            // if slope_minmod == avg_slope_k just set sl_new_Qs[i] = new_Qs[i]
            if (slope_minmod != avg_slope_k) {
                sl_new_Qs[i](0) = q_avg_k;
                sl_new_Qs[i](1) = (slope_minmod*dx)/2.0;
                for (int n = 2; n < sl_new_Qs[i].size(); n++) {
                    sl_new_Qs[i](n) = 0;
                }
            }         
        }*/

        // exchange new_Qs with slope_limited_new_Qs
        //new_Qs = sl_new_Qs;
        new_Qs = slope_limit_adv1D(new_Qs, slope_limited);

    }
    
    // apply new Q values
    for (int i = 0; i < grid->cells.size(); i++) {
        grid->cells[i].Q = new_Qs[i];
    }
}


// function to slope limit the coefficients of the 1D advection
template <typename CellType>
vector<VectorXd> DG_Solver<CellType>::slope_limit_adv1D(vector<VectorXd> new_Qs, double slope_limited) {

    vector<VectorXd> sl_new_Qs;
    sl_new_Qs.reserve(grid->cells.size());
    sl_new_Qs = new_Qs;

    // go through all elements
    for (int i = 0; i < grid->cells.size(); i++) {

        // get averaged slope on T_k
        double dx = grid->cells[i].edges[1].b.x - grid->cells[i].edges[1].a.x;
        double dy = 0;
        for (int a = 0; a < grid->N_basisfunc; a++) {
            if (a%2 != 0) {dy += 2*new_Qs[i](a);}
        }
        double avg_slope_k = dy/dx;

        // get averaged q values for T_k-1, T_k and T_k+1
        double q_avg_m = (i >= 1) ? get_element_average_adv1D(i-1, new_Qs) : 0;
        double q_avg_k = get_element_average_adv1D(i, new_Qs);
        double q_avg_p = (i <= new_Qs.size()-2) ? get_element_average_adv1D(i+1, new_Qs) : get_element_average_adv1D(i, new_Qs);    // note that we do respect the inflow and outflow BC here

        // calc averaged slopes between T_k and T_k-1 as well as T_k+1 and T_k
        double avg_slope_m = slope_limited*(q_avg_k - q_avg_m)/dx;
        double avg_slope_p = slope_limited*(q_avg_p - q_avg_k)/dx;

        // do minmod limiting of averaged slopes
        double slope_minmod = minmod(avg_slope_k, avg_slope_m, avg_slope_p);

        // if slope_minmod == avg_slope_k just set sl_new_Qs[i] = new_Qs[i]
        if (slope_minmod != avg_slope_k) {
            sl_new_Qs[i](0) = q_avg_k;
            sl_new_Qs[i](1) = (slope_minmod*dx)/2.0;
            for (int n = 2; n < sl_new_Qs[i].size(); n++) {
                sl_new_Qs[i](n) = 0;
            }
        }

    }

    return sl_new_Qs;
}


// function to apply minmod limiting to three values
template <typename CellType>
double DG_Solver<CellType>::minmod(double a, double b, double c) {

    if (a/abs(a) == b/abs(b) && b/abs(b) == c/abs(c)) {

        double s = a/abs(a);
        return s * min({abs(a), abs(b), abs(c)});

    } else {
        return 0;
    }
}

// returns element average in 1D advection monomial basis
template <typename CellType>
double DG_Solver<CellType>::get_element_average_adv1D(int index, vector<VectorXd>& new_Qs) {

    double q_avg = 0;

    for (int i = 0; i < grid->N_basisfunc; i++) {
        if (i%2 == 0) {
            q_avg += (1.0/(static_cast<double>(i)+1.0)) * new_Qs[index](i);
        }
    }

    return q_avg;
}

// - - - - - - - 2D scalar upwind advection - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// assuming u.x and u.y > 0
// function to calculate next step in scalar upwind 2D advection discontinous galerkin
template <typename CellType>
void DG_Solver<CellType>::advection2D(double dt, double slope_limited, Point u, bool intitialize) {

    if (intitialize) {
        for (int k = 0; k < grid->cells.size(); k++) {
            // calculate M^k_ij inverse - also the N_row is absorbed in here
            grid->cells[k].M_ij_k_inv = grid->cells[k].calc_M_ij_k_inv_adv2D();
            // calculate S^k_ij (absorbing a u in there for later on)
            grid->cells[k].S_ij_k = grid->cells[k].calc_S_ij_k_adv2D(u);
        }
        F_hat_A = get_F_hat_A_adv2D();
        F_hat_B = get_F_hat_B_adv2D();
        F_hat_C = get_F_hat_C_adv2D();
        F_hat_D = get_F_hat_D_adv2D();
    }

    // vector to store new Q values
    vector<VectorXd> new_Qs;
    new_Qs.reserve(grid->cells.size());

    VectorXd zeroQ(grid->N_basisfunc);
    zeroQ.setZero();

    // calculate new Q values
    // loop through all elements T_k for that
    for (int k = 0; k < grid->cells.size(); k++) {

        // calculate RHS_I
        VectorXd RHS_I = grid->cells[k].S_ij_k * grid->cells[k].Q;

    
        // get neighbour element coefficients
        VectorXd Q_km1;
        if (grid->cells[k].edges[0].is_boundary) {
            Q_km1 = zeroQ;
        } else {
            Q_km1 = grid->cells[k].edges[0].neighbour->get_coeff();
        }
        
        VectorXd Q_kmN;
        if (grid->cells[k].edges[3].is_boundary) {
            Q_kmN = zeroQ;
        } else {
            Q_kmN = grid->cells[k].edges[3].neighbour->get_coeff();
        }

        //VectorXd RHS = RHS_I + flux1 + flux2 + flux3 + flux4;
        VectorXd RHS = RHS_I + (u.x * F_hat_A * Q_km1) 
                             - (u.x * F_hat_B * grid->cells[k].Q) 
                             - (u.y * F_hat_C * grid->cells[k].Q) 
                             + (u.y * F_hat_D * Q_kmN);

        // calculate new Q using forward euler
        VectorXd Q_k_np1 = grid->cells[k].Q + dt * grid->cells[k].M_ij_k_inv * RHS; 

        // store it in new Q vector
        new_Qs.push_back(Q_k_np1);

    }

    // optional slope limiting as a post processing step
    if (grid->N_basisfunc > 1 && slope_limited != 0) {
        new_Qs = slope_limit_adv2D(new_Qs, slope_limited);
    }

    // apply new Q values
    for (int i = 0; i < grid->cells.size(); i++) {
        grid->cells[i].Q = new_Qs[i];
    }

}

// returns F_hat_A for 2d cartesian upwind advection
template <typename CellType>
MatrixXd DG_Solver<CellType>::get_F_hat_A_adv2D() {

    if (grid->N_basisfunc != 1 && grid->N_basisfunc != 3 && grid->N_basisfunc != 6) {
        cerr << "Wrong number of basis functions for F_hat_A_adv2D. Please change to 1 or 3 or 6." << endl;
        EXIT_FAILURE;
    }

    MatrixXd F_hat_A(grid->N_basisfunc, grid->N_basisfunc);
    F_hat_A.setZero();
    F_hat_A(0, 0) = 1;
    if (grid->N_basisfunc >= 3) {
        F_hat_A(0, 1) = 1;
        F_hat_A(1, 0) = -1;
        F_hat_A(1, 1) = -1;
        F_hat_A(2, 2) = 1.0/3.0;
    }
    if (grid->N_basisfunc >= 6) {
        F_hat_A(2, 3) = 1.0/3.0;
        F_hat_A(3, 2) = -1.0/3.0;
        F_hat_A(3, 3) = -1.0/3.0;
        F_hat_A(4, 0) = 1;
        F_hat_A(4, 1) = 1;
        F_hat_A(0, 4) = 1;
        F_hat_A(4, 4) = 1;
        F_hat_A(1, 4) = -1;
        F_hat_A(5, 5) = 1.0/5.0;

    }
    
    return F_hat_A;
}

// returns F_hat_B for 2d cartesian upwind advection
template <typename CellType>
MatrixXd DG_Solver<CellType>::get_F_hat_B_adv2D() {

    if (grid->N_basisfunc != 1 && grid->N_basisfunc != 3 && grid->N_basisfunc != 6) {
        cerr << "Wrong number of basis functions for F_hat_B_adv2D. Please change to 1 or 3 or 6." << endl;
        EXIT_FAILURE;
    }

    MatrixXd F_hat_B(grid->N_basisfunc, grid->N_basisfunc);
    F_hat_B.setZero();
    F_hat_B(0, 0) = 1;
    if (grid->N_basisfunc >= 3) {
        F_hat_B(0, 1) = 1;
        F_hat_B(1, 0) = 1;
        F_hat_B(1, 1) = 1;
        F_hat_B(2, 2) = 1.0/3.0;
    }
    if (grid->N_basisfunc >= 6) {
        F_hat_B(2, 3) = 1.0/3.0;
        F_hat_B(3, 2) = 1.0/3.0;
        F_hat_B(3, 3) = 1.0/3.0;
        F_hat_B(4, 0) = 1;
        F_hat_B(4, 1) = 1;
        F_hat_B(0, 4) = 1;
        F_hat_B(4, 4) = 1;
        F_hat_B(1, 4) = 1;
        F_hat_B(5, 5) = 1.0/5.0;
    }
    
    return F_hat_B;
}

// returns F_hat_C for 2d cartesian upwind advection
template <typename CellType>
MatrixXd DG_Solver<CellType>::get_F_hat_C_adv2D() {

    if (grid->N_basisfunc != 1 && grid->N_basisfunc != 3 && grid->N_basisfunc != 6) {
        cerr << "Wrong number of basis functions for F_hat_C_adv2D. Please change to 1 or 3 or 6." << endl;
        EXIT_FAILURE;
    }

    MatrixXd F_hat_C(grid->N_basisfunc, grid->N_basisfunc);
    F_hat_C.setZero();
    F_hat_C(0, 0) = 1;
    if (grid->N_basisfunc >= 3) {
        F_hat_C(1, 1) = 1.0/3.0;
        F_hat_C(0, 2) = 1;
        F_hat_C(2, 0) = 1;
        F_hat_C(2, 2) = 1;
    }
    if (grid->N_basisfunc >= 6) {
        F_hat_C(1, 3) = 1.0/3.0;
        F_hat_C(3, 1) = 1.0/3.0;
        F_hat_C(3, 3) = 1.0/3.0;
        F_hat_C(5, 0) = 1;
        F_hat_C(0, 5) = 1;
        F_hat_C(5, 2) = 1;
        F_hat_C(2, 5) = 1;
        F_hat_C(5, 5) = 1;
        F_hat_C(4, 4) = 1.0/5.0;
    }
    
    return F_hat_C;
}

// returns F_hat_D for 2d cartesian upwind advection
template <typename CellType>
MatrixXd DG_Solver<CellType>::get_F_hat_D_adv2D() {

    if (grid->N_basisfunc != 1 && grid->N_basisfunc != 3 && grid->N_basisfunc != 6) {
        cerr << "Wrong number of basis functions for F_hat_D_adv2D. Please change to 1 or 3 or 6." << endl;
        EXIT_FAILURE;
    }

    MatrixXd F_hat_D(grid->N_basisfunc, grid->N_basisfunc);
    F_hat_D.setZero();
    F_hat_D(0, 0) = 1;
    if (grid->N_basisfunc >= 3) {
        F_hat_D(1, 1) = 1.0/3.0;
        F_hat_D(0, 2) = 1;
        F_hat_D(2, 0) = -1;
        F_hat_D(2, 2) = -1;
    }
    if (grid->N_basisfunc >= 6) {
        F_hat_D(1, 3) = 1.0/3.0;
        F_hat_D(3, 1) = -1.0/3.0;
        F_hat_D(3, 3) = -1.0/3.0;
        F_hat_D(5, 0) = 1;
        F_hat_D(0, 5) = 1;
        F_hat_D(5, 2) = 1;
        F_hat_D(2, 5) = -1;
        F_hat_D(5, 5) = 1;
        F_hat_D(4, 4) = 1.0/5.0;
    }
    
    return F_hat_D;
}


// slope limits 2D scalar advection with minmod while handling x and y direction independently, only works for N_bf = 1 or 3
template <typename CellType>
vector<VectorXd> DG_Solver<CellType>::slope_limit_adv2D(vector<VectorXd> new_Qs, double slope_limited) {

    if (grid->N_basisfunc != 3 && grid->N_basisfunc != 6) {
        cerr << "slope limiter does only work with 3 or 6 basis func, please change N_basisfunc" << endl;
        EXIT_FAILURE;
    }

    vector<VectorXd> sl_new_Qs;
    sl_new_Qs.reserve(grid->cells.size());
    sl_new_Qs = new_Qs;

    // go through all elements
    for (int i = 0; i < grid->cells.size(); i++) {

        // get averaged slopes on T_k (in x and y direction)
        double dx = grid->cells[i].edges[1].b.x - grid->cells[i].edges[1].a.x;
        double avg_slope_k_x = 2*new_Qs[i](1)/dx;
        double avg_slope_k_y = 2*new_Qs[i](2)/dx;


        // get averaged q values for T_k and neighbours
        double q_avg_k = new_Qs[i](0);
        double q_avg_l = (i >= 1) ? new_Qs[i-1](0) : 0; // consider inflow bc
        double q_avg_d = (i >= grid->n_horizontal) ? new_Qs[i-grid->n_horizontal](0) : 0; // consider inflow bc
        double q_avg_r = (i <= new_Qs.size()-2) ? new_Qs[i+1](0) : q_avg_k; // consider outflow bc
        double q_avg_u = (i <= new_Qs.size()-grid->n_horizontal-1) ? new_Qs[i+grid->n_horizontal](0) : q_avg_k; // consider outflow bc

        // calc averaged slopes between T_k and neighbours
        double avg_slope_l = slope_limited * (q_avg_k - q_avg_l)/dx;
        double avg_slope_r = slope_limited * (q_avg_r - q_avg_k)/dx;
        double avg_slope_u = slope_limited * (q_avg_u - q_avg_k)/dx;
        double avg_slope_d = slope_limited * (q_avg_k - q_avg_d)/dx;

        // do minmod limiting of averaged slopes
        double slope_minmod_x = minmod(avg_slope_k_x, avg_slope_l, avg_slope_r);
        double slope_minmod_y = minmod(avg_slope_k_y, avg_slope_d, avg_slope_u);

        // if slope_minmod != avg_slope then reconstruct limited with linear (x and y reconstructed linearly here)
        if (slope_minmod_x != avg_slope_k_x) {
            sl_new_Qs[i](0) = q_avg_k;
            sl_new_Qs[i](1) = (slope_minmod_x*dx)/2.0;
            if (sl_new_Qs[i].size() >= 6) {
                //sl_new_Qs[i](3) = 0;
                sl_new_Qs[i](4) = 0;
            }
        }
        if (slope_minmod_y != avg_slope_k_y) {
            sl_new_Qs[i](0) = q_avg_k;
            sl_new_Qs[i](2) = (slope_minmod_y*dx)/2.0;
            if (sl_new_Qs[i].size() >= 6) {
                //sl_new_Qs[i](3) = 0;
                sl_new_Qs[i](5) = 0;
            }
        }

    }

    return sl_new_Qs;
}


template <typename CellType>
DG_Solver<CellType>::~DG_Solver() {};