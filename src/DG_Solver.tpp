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

        // calculate M^k_ij
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
            double q_avg_m = (i >= 1) ? get_element_average_adv1D(i-1, new_Qs) : 0;                             // again for N_p = 1 the averages are quite straight forward
            double q_avg_k = get_element_average_adv1D(i, new_Qs);
            double q_avg_p = (i <= new_Qs.size()-2) ? get_element_average_adv1D(i+1, new_Qs) : get_element_average_adv1D(i, new_Qs);    // note that we do respect the inflow and outflow BC here


            // calc averaged slopes between T_k and T_k-1 as well as T_k+1 and T_k
            double avg_slope_m = slope_limited*(q_avg_k - q_avg_m)/dx;
            double avg_slope_p = slope_limited*(q_avg_p - q_avg_k)/dx;

            // do minmod limiting of averaged slopes
            double slope_minmod = minmod(avg_slope_k, avg_slope_m, avg_slope_p);

            //cout << avg_slope_k <<"   " << new_Qs[i](0) << "      " << new_Qs[i](1) << endl;

            // if slope_minmod == avg_slope_k just set sl_new_Qs[i] = new_Qs[i]
            if (slope_minmod != avg_slope_k) {
                sl_new_Qs[i](0) = q_avg_k;
                sl_new_Qs[i](1) = (slope_minmod*dx)/2.0;
                for (int n = 2; n < sl_new_Qs[i].size(); n++) {
                    sl_new_Qs[i](n) = 0;
                }
            }         
        }

        // exchange new_Qs with slope_limited_new_Qs
        new_Qs = sl_new_Qs;

    }
    /*
    // optional slope limiting as a post processing step
    // first try only for N_p = 1 <-> N_bf = 2
    if (grid->N_basisfunc == 2 && slope_limited != 0) {

        vector<VectorXd> sl_new_Qs;
        sl_new_Qs.reserve(grid->cells.size());
        sl_new_Qs = new_Qs;

        // go through all elements
        for (int i = 0; i < grid->cells.size(); i++) {

            // get averaged slope on T_k (for N_p = 1 the slope is constant so already averaged)
            double dx = grid->cells[i].edges[1].b.x - grid->cells[i].edges[1].a.x;
            double dy = 2*new_Qs[i](1);
            double avg_slope_k = dy/dx;

            // get averaged q values for T_k-1, T_k and T_k+1
            double q_avg_m = (i >= 1) ? new_Qs[i-1](0) : 0;                             // again for N_p = 1 the averages are quite straight forward
            double q_avg_k = new_Qs[i](0);
            double q_avg_p = (i <= new_Qs.size()-2) ? new_Qs[i+1](0) : new_Qs[i](0);    // note that we do respect the inflow and outflow BC here

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

        // exchange new_Qs with slope_limited_new_Qs
        new_Qs = sl_new_Qs;

    }*/

    // apply new Q values
    for (int i = 0; i < grid->cells.size(); i++) {
        grid->cells[i].Q = new_Qs[i];
    }
}

template <typename CellType>
double DG_Solver<CellType>::minmod(double a, double b, double c) {

    if (a/abs(a) == b/abs(b) && b/abs(b) == c/abs(c)) {

        double s = a/abs(a);
        return s * min({abs(a), abs(b), abs(c)});

    } else {
        return 0;
    }
}

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



template <typename CellType>
DG_Solver<CellType>::~DG_Solver() {};