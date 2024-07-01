#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "cell_types/Q_Cell.h"
#include "Mesh.h"
#include "Solver.h"

template <typename CellType>
Solver<CellType>::Solver() {}

template <typename CellType>
Solver<CellType>::Solver(Mesh<CellType>* gridin) {
    grid = gridin;
}

// PRIVATE helper functions -----------------------------------------------------------------------
// function to get the normal vector on a face definde by two points
template <typename CellType>
Point Solver<CellType>::get_normal_vec(Point a, Point b) {

    double Delta_x = b.x - a.x;
    double Delta_y = b.y - a.y;

    double norm = sqrt(Delta_x*Delta_x + Delta_y*Delta_y);

    return Point(-Delta_y/norm, Delta_x/norm);
}

// SOLVERS (update the Mesh by a single step) -----------------------------------------------------
// SOLVER: Fake Diffusion--------------------------------------------------------------------------
// Fake Diffusion Equation as first try, requires Q_cells, Q is a absolute value not a density
template <typename CellType>
void Solver<CellType>::diffusion_like(double dt) {
 
    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Q_Cell> == false) {
        cerr << "diffusion_like step called on Mesh with wrong cell type, you must use Q_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    vector<double> Qs_new;
    Qs_new.reserve(grid->cells.size());

    // go through all cells and calculate total in/outflow
    for (int i = 0; i < grid->cells.size(); i++) {

        double deltaQ = 0;

        // calculate flux for each face and add to total in/outflow
        for (int j = 0; j< grid->cells[i].edges.size(); j++) {

            // calculate flux for normal faces and boundary faces
            if (grid->cells[i].edges[j].is_boundary == false) {

                // diffusion like flux
                deltaQ += 0.1 * (grid->cells[i].edges[j].neighbour->getQ() - grid->cells[i].Q);

            } else {
                deltaQ += 0;  // no flux through boundaries at the moment
            }
        }
        // new Q = old Q + delta Q
        Qs_new.push_back(grid->cells[i].Q + deltaQ);
    }

    double total = 0;

    // after calculating now exchange Q_new for Q_old
    for (int i = 0; i < grid->cells.size(); i++) {
        grid->cells[i].Q = Qs_new[i];
        total += Qs_new[i];
    }
}


// SOLVER: Conaway's Game of Life -----------------------------------------------------------------
// Just a side joke, Evolution Step for Conway's Game of Life, requires conway_cell
template <typename CellType>
void Solver<CellType>::conway() {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Conway_Cell> == false) {
        cerr << "conway step called on Mesh with wrong cell type, you must use Conway_Cells" << endl;
        exit(EXIT_FAILURE);
    }


    vector<int> Qs_new;
    Qs_new.reserve(grid->cells.size());

    // go through all cells and decide whether they will live or die
    for (int i = 0; i < grid->cells.size(); i++) {

        int neighbour_count = 0;

        // count the living neighbours
        for (int j = 0; j< grid->cells[i].edges.size(); j++) {

            if (grid->cells[i].edges[j].is_boundary == false) {

                neighbour_count += grid->cells[i].edges[j].neighbour->getQ();

            }

        }

        // also count cell on top left
        if (grid->cells[i].edges.size()>=2) {
            if (grid->cells[i].edges[1].is_boundary == false) {
                if (grid->cells[i].edges[1].neighbour->edges.size()>=1) {
                    if (grid->cells[i].edges[1].neighbour->edges[0].is_boundary == false) {
                        neighbour_count += grid->cells[i].edges[1].neighbour->edges[0].neighbour->getQ();
                    }
                }
            }
        }

        // also count cell on top right
        if (grid->cells[i].edges.size()>=2) {
            if (grid->cells[i].edges[1].is_boundary == false) {
                if (grid->cells[i].edges[1].neighbour->edges.size()>=3) {
                    if (grid->cells[i].edges[1].neighbour->edges[2].is_boundary == false) {
                        neighbour_count += grid->cells[i].edges[1].neighbour->edges[2].neighbour->getQ();
                    }
                }
            }
        }   

        // also count cell on bottom left
        if (grid->cells[i].edges.size()>=4) {
            if (grid->cells[i].edges[3].is_boundary == false) {
                if (grid->cells[i].edges[3].neighbour->edges.size()>=1) {
                    if (grid->cells[i].edges[3].neighbour->edges[0].is_boundary == false) {
                        neighbour_count += grid->cells[i].edges[3].neighbour->edges[0].neighbour->getQ();
                    }
                }
            }
        }  

        // also count cell on bottom right
        if (grid->cells[i].edges.size()>=4) {
            if (grid->cells[i].edges[3].is_boundary == false) {
                if (grid->cells[i].edges[3].neighbour->edges.size()>=3) {
                    if (grid->cells[i].edges[3].neighbour->edges[2].is_boundary == false) {
                        neighbour_count += grid->cells[i].edges[3].neighbour->edges[2].neighbour->getQ();
                    }
                }
            }
        }    

        // rules of conways game of life
        if (neighbour_count<=1 || neighbour_count >=4) {
            Qs_new.push_back(0);
        } else if (neighbour_count == 3) {
            Qs_new.push_back(1);
        } else {
            Qs_new.push_back(grid->cells[i].Q);
        }

    }

    // after calculating now exchange Q_new for Q_old
    for (int i = 0; i < grid->cells.size(); i++) {
        grid->cells[i].Q = Qs_new[i];
    }

}


// SOLVER: ADVECTION ------------------------------------------------------------------------------
// FV Upwind Advection on cartesian or Vmesh, requires Q_Cells
template <typename CellType>
void Solver<CellType>::advection(double dt, Point v) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, Q_Cell> == false) {
        cerr << "advection step called on Mesh with wrong cell type, you must use Q_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    if (grid->is_cartesian) {
        advection_cartesian(dt, v);
    } else {
        advection_vmesh(dt, v);
    }
}

// FV Upwind Advection, requires Q_Cells and cartesian mesh
template <typename CellType>
void Solver<CellType>::advection_cartesian(double dt, Point v) {

    if (dt * abs(v.x) >= grid->cells[0].edges[1].length || dt * abs(v.y) >= grid->cells[0].edges[2].length) {
        cerr << "the CFL timestep condition isnt fullfilled, please make sure that dt < dx/v. otherwise the simulation will show unstable behaviour" << endl;
        exit(EXIT_FAILURE);
    }

    vector<double> new_Qs;
    new_Qs.reserve(grid->cells.size());

    // calculate all the new Qs
    for (int i = 0; i<grid->cells.size(); i++) {

        double Qi_np1 = 0;
        double Qi_n = grid->cells[i].Q;
        double Qim1x_n = 0;
        double Qim1y_n = 0;
        double dx = grid->cells[i].edges[1].length;
        double dy = grid->cells[i].edges[2].length;

        // find upwind Q values
        if (v.x > 0 && grid->cells[i].edges[0].is_boundary == false) {
            Qim1x_n = grid->cells[i].edges[0].neighbour->getQ();
        }
        if (v.x < 0 && grid->cells[i].edges[2].is_boundary == false) {
            Qim1x_n = grid->cells[i].edges[2].neighbour->getQ();
        }
        if (v.y > 0 && grid->cells[i].edges[3].is_boundary == false) {
            Qim1y_n = grid->cells[i].edges[3].neighbour->getQ();
        }
        if (v.y < 0 && grid->cells[i].edges[1].is_boundary == false) {
            Qim1y_n = grid->cells[i].edges[1].neighbour->getQ();
        }

        // calculate new Q
        Qi_np1 = Qi_n - abs(v.x) * (dt/dx) * (Qi_n - Qim1x_n) - abs(v.y) * (dt/dy) * (Qi_n - Qim1y_n);

        new_Qs.push_back(Qi_np1);

    }

    // apply all the new Qs
    for (int i = 0; i<grid->cells.size(); i++) {

        grid->cells[i].Q = new_Qs[i];

    }

}

// FV Upwind Advection on Vmesh, requires Q_Cells
template <typename CellType>
void Solver<CellType>::advection_vmesh(double dt, Point v) {

    vector<double> new_Qs;
    new_Qs.reserve(grid->cells.size());

    // calculate all the new Qs
    for (int i = 0; i < grid->cells.size(); i++) {

        // get different variables we'll later need
        double Q_i_np1 = 0;
        double Q_i_n = grid->cells[i].Q;
        double total_Flux = 0;
        double dV = grid->cells[i].volume;

        // calculate fluxes for all faces and add to total_flux
        for (int j = 0; j < grid->cells[i].edges.size(); j++) {

            // calculate v_eff and get l_ij
            Point n_ij;
            n_ij = get_normal_vec(grid->cells[i].edges[j].a, grid->cells[i].edges[j].b);
            double v_eff = (v.x*n_ij.x + v.y*n_ij.y);
            double l_ij = grid->cells[i].edges[j].length;
            
            // get Q_j_n
            double Q_j_n;
            if (grid->cells[i].edges[j].is_boundary) {
                Q_j_n = 0;
            } else {
                Q_j_n = grid->cells[i].edges[j].neighbour->getQ();
            }

            // add to total flux
            if (v_eff>=0) {
                total_Flux += v_eff * Q_i_n * l_ij;
            } else {
                total_Flux += v_eff * Q_j_n * l_ij;
            }


        }

        // calculate next Q
        Q_i_np1 = Q_i_n - dt/dV *(total_Flux);

        new_Qs.push_back(Q_i_np1);

    }


    // apply the calculated Qs
    for (int i = 0; i<grid->cells.size(); i++) {

        grid->cells[i].Q = static_cast<double>(new_Qs[i]);

    }


}


// SOLVER: Shallow Water Equations ----------------------------------------------------------------
// 1st order implementation of shallow water equation 2D voronoi using HLL solver
template<typename CellType>
void Solver<CellType>::shallow_water(double dt, int boundary_cond, double numerical_diffusion_coeff) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, SWE_Cell> == false) {
        cerr << "shallow water step called on Mesh with wrong cell type, you must use SWE_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    vector<array<double, 3>> Q_new;
    Q_new.reserve(grid->cells.size());

    //go through each cell to calculate the new values
    for (int i = 0; i<grid->cells.size(); i++) {

        //cout << "Cell: " << i << endl;

        double A = grid->cells[i].volume;

        array<double, 3> huv_i_np1 = {0, 0, 0};
        array<double, 3> huv_i_n = {grid->cells[i].h, grid->cells[i].u, grid->cells[i].v};
        array<double, 3> total_flux = {0, 0, 0};

        // calculate total_flux by summing over edge fluxed * edge_length
        for (int j = 0; j < grid->cells[i].edges.size(); j++) {
            
            Point n = get_normal_vec(grid->cells[i].edges[j].a, grid->cells[i].edges[j].b);
            double l_i_j = grid->cells[i].edges[j].length;
            double g = 1;

            // values needed for f and g calculations
            // values of other cell (initalized with sink boundaries for boundary cond = 1)
            array<double, 3> huv_j_n = {huv_i_n[0], boundary_cond * huv_i_n[1],  boundary_cond * huv_i_n[2]};
            // if reflective boundary change velocities accordingly
            if (boundary_cond == -1) {
                huv_j_n[1] = huv_i_n[1] - 2*(huv_i_n[1]*n.x + huv_i_n[2]*n.y)*n.x;
                huv_j_n[2] = huv_i_n[2] - 2*(huv_i_n[1]*n.x + huv_i_n[2]*n.y)*n.y;
            }

            // get actual huv values if its not a boundary
            if (grid->cells[i].edges[j].is_boundary == false) {
                huv_j_n[0] = grid->cells[i].edges[j].neighbour->get_h();
                huv_j_n[1] = grid->cells[i].edges[j].neighbour->get_u();
                huv_j_n[2] = grid->cells[i].edges[j].neighbour->get_v();
            }

            // f_i
            array<double, 3> f_i = {
                                  huv_i_n[0] * huv_i_n[1],
                                  huv_i_n[0] * huv_i_n[1] * huv_i_n[1] + (g/2.0) * huv_i_n[0] * huv_i_n[0],
                                  huv_i_n[0] * huv_i_n[2] * huv_i_n[1]
                                 };
            
            // f_j
            array<double, 3> f_j = {
                                  huv_j_n[0] * huv_j_n[1], 
                                  huv_j_n[0] * huv_j_n[1] * huv_j_n[1] + (g/2.0) * huv_j_n[0] * huv_j_n[0], 
                                  huv_j_n[0] * huv_j_n[2] * huv_j_n[1]
                                 };
            
            // g_i
            array<double, 3> g_i = {
                                  huv_i_n[0] * huv_i_n[2], 
                                  huv_i_n[0] * huv_i_n[2] * huv_i_n[1], 
                                  huv_i_n[0] * huv_i_n[2] * huv_i_n[2] + (g/2.0) * huv_i_n[0] * huv_i_n[0]
                                     };
            // g_j
            array<double, 3> g_j = {
                                  huv_j_n[0] * huv_j_n[2], 
                                  huv_j_n[0] * huv_j_n[2] * huv_j_n[1], 
                                  huv_j_n[0] * huv_j_n[2] * huv_j_n[2] + (g/2.0) * huv_j_n[0] * huv_j_n[0]
                                 };

            // Flux Estimation using Roe approximate Riemann Solver
            array<double, 3> F_eff;

            // get effective Flux using a split scheme HLL solver
            F_eff = hll_solver_swe_2D(huv_i_n, huv_j_n, f_i, f_j, g_i, g_j, g, n, numerical_diffusion_coeff);

            // add to total flux * length
            total_flux[0] += F_eff[0] * l_i_j;
            total_flux[1] += F_eff[1] * l_i_j;
            total_flux[2] += F_eff[2] * l_i_j;
        }

        // calculate new values
        huv_i_np1[0] = huv_i_n[0] - (dt/A) * total_flux[0];
        huv_i_np1[1] = (huv_i_n[0]*huv_i_n[1] - (dt/A) * total_flux[1])/huv_i_np1[0];
        huv_i_np1[2] = (huv_i_n[0]*huv_i_n[2] - (dt/A) * total_flux[2])/huv_i_np1[0];

        // store new values
        Q_new.push_back(huv_i_np1);

    }

    // then apply all the changes
    for (int i = 0; i<grid->cells.size(); i++) {
        grid->cells[i].h = Q_new[i][0];
        grid->cells[i].u = Q_new[i][1];
        grid->cells[i].v = Q_new[i][2];
    }

}

// RIEMANN SOLVERS --------------------------------------------------------------------------------
// UNPHYSICAL FOR WHATEVER REASON: Roe's Solver for Shallow Water Equation 2D unstructured
template <typename CellType>
array<double, 3> Solver<CellType>::roe_solver_swe_2D(array<double, 3> huv_i_n, array<double, 3> huv_j_n, array<double, 3> f_i, array<double, 3> f_j, array<double, 3> g_i, array<double, 3> g_j, double g, Point n, double add_diffusion_coeff) {

    cerr << "WARNING: using unphysical roe_solver. Is that intended. Maybe instead use HLL solver?" << endl;

    double u_roe = (huv_i_n[1]*sqrt(huv_i_n[0]) + huv_j_n[1]*sqrt(huv_j_n[0]))/(sqrt(huv_i_n[0]) + sqrt(huv_j_n[0]));
    double v_roe = (huv_i_n[2]*sqrt(huv_i_n[0]) + huv_j_n[2]*sqrt(huv_j_n[0]))/(sqrt(huv_i_n[0]) + sqrt(huv_j_n[0]));
    double c_roe = sqrt((g/2.0) * (huv_i_n[0] + huv_j_n[0]));

    // eigenvalues and determinant
    double lambda1 = u_roe * n.x + v_roe * n.y;
    double lambda2 = u_roe * n.x + v_roe * n.y - c_roe;
    double lambda3 = u_roe * n.x + v_roe * n.y + c_roe;

    double detA = lambda1 * lambda2 * lambda3;

    // get effective fluxes
    array<double, 3> F_eff = {(0.5) * (((f_i[0]*n.x + g_i[0]*n.y) + (f_j[0]*n.x + g_j[0]*n.y)) - detA * (huv_i_n[0] - huv_j_n[0])) - add_diffusion_coeff * (huv_i_n[0] - huv_j_n[0]), 
                            (0.5) * (((f_i[1]*n.x + g_i[1]*n.y) + (f_j[1]*n.x + g_j[1]*n.y)) - detA * (huv_i_n[0] * huv_i_n[1] - huv_j_n[0] * huv_j_n[1]))- add_diffusion_coeff * (huv_i_n[0] * huv_i_n[1] - huv_j_n[0] * huv_j_n[1]), 
                            (0.5) * (((f_i[2]*n.x + g_i[2]*n.y) + (f_j[2]*n.x + g_j[2]*n.y)) - detA * (huv_i_n[0] * huv_i_n[2] - huv_j_n[0] * huv_j_n[2]))- add_diffusion_coeff * (huv_i_n[0] * huv_i_n[2] - huv_j_n[0] * huv_j_n[2])
                            };
    return F_eff;
}



// CONSTRUCTION SITE: -----------------------------------------------------------------------------
// HLL Solver for Shallow Water Equation 2D unstructured
template <typename CellType>
array<double, 3> Solver<CellType>::hll_solver_swe_2D(array<double, 3> huv_i_n, array<double, 3> huv_j_n, array<double, 3> f_i, array<double, 3> f_j, array<double, 3> g_i, array<double, 3> g_j, double g, Point n, double add_diffusion_coeff) {
    
    // solves Riemann problem in a split scheme
    // start with x direction and flux F

    // make sure that left and right state are chosen correct for F
    array<double, 3> huv_l_n;
    array<double, 3> huv_r_n; 
    array<double, 3> f_l; 
    array<double, 3> f_r;
    if (n.x > 0) {
        huv_l_n = huv_i_n;
        huv_r_n = huv_j_n;
        f_l = f_i;
        f_r = f_j;
    } else {
        huv_l_n = huv_j_n;
        huv_r_n = huv_i_n;
        f_l = f_j;
        f_r = f_i;
    }
    
    // calculate wave speeds
    double SLF = min(huv_l_n[1] - sqrt(g*huv_l_n[0]), huv_r_n[1] - sqrt(g*huv_r_n[0]));
    double SRF = max(huv_r_n[1] + sqrt(g*huv_r_n[0]), huv_l_n[1] + sqrt(g*huv_l_n[0]));

    // HLL solver for F
    array<double, 3> F_HLL;
    if (SLF >= 0) {
        F_HLL = f_l;
    } else if (SLF < 0 && SRF > 0) {
        F_HLL[0] = (SRF * f_l[0] - SLF * f_r[0] + SLF * SRF * (huv_r_n[0] - huv_l_n[0])) / (SRF - SLF);
        F_HLL[1] = (SRF * f_l[1] - SLF * f_r[1] + SLF * SRF * (huv_r_n[0] * huv_r_n[1] - huv_l_n[0] * huv_l_n[1])) / (SRF - SLF);
        F_HLL[2] = (SRF * f_l[2] - SLF * f_r[2] + SLF * SRF * (huv_r_n[0] * huv_r_n[2] - huv_l_n[0] * huv_l_n[2])) / (SRF - SLF);
    } else if (SRF <= 0) {
        F_HLL = f_r;
    }

    // now solve Riemann problem in y direction for flux G
    // make sure that left and right state are chosen correct for G
    array<double, 3> g_l;
    array<double, 3> g_r;
    if (n.y > 0) {
        huv_l_n = huv_i_n;
        huv_r_n = huv_j_n;
        g_l = g_i;
        g_r = g_j;
    } else {
        huv_l_n = huv_j_n;
        huv_r_n = huv_i_n;
        g_l = g_j;
        g_r = g_i;
    }

    // calculate wave speeds
    double SLG = min(huv_l_n[1] - sqrt(g*huv_l_n[0]), huv_r_n[1] - sqrt(g*huv_r_n[0]));
    double SRG = max(huv_r_n[1] + sqrt(g*huv_r_n[0]), huv_l_n[1] + sqrt(g*huv_l_n[0]));

    // HLL Solver for G
    array<double, 3> G_HLL;
    if (SLG >= 0) {
        G_HLL = g_l;
    } else if (SLG < 0 && SRG > 0) {
        G_HLL[0] = (SRG * g_l[0] - SLG * g_r[0] + SLG * SRG * (huv_r_n[0] - huv_l_n[0])) / (SRG - SLG);
        G_HLL[1] = (SRG * g_l[1] - SLG * g_r[1] + SLG * SRG * (huv_r_n[0] * huv_r_n[1] - huv_l_n[0] * huv_l_n[1])) / (SRG - SLG);
        G_HLL[2] = (SRG * g_l[2] - SLG * g_r[2] + SLG * SRG * (huv_r_n[0] * huv_r_n[2] - huv_l_n[0] * huv_l_n[2])) / (SRG - SLG);
    } else if (SRG <= 0) {
        G_HLL = g_r;
    }

    // now combine split scheme + optional additional diffusion (not needed normally)
    array<double, 3> F_eff = {
        (F_HLL[0]*n.x + G_HLL[0]*n.y) - add_diffusion_coeff * (huv_l_n[0] - huv_r_n[0]),
        (F_HLL[1]*n.x + G_HLL[1]*n.y) - add_diffusion_coeff * (huv_l_n[0] * huv_l_n[1] - huv_r_n[0] * huv_r_n[1]),
        (F_HLL[2]*n.x + G_HLL[2]*n.y) - add_diffusion_coeff * (huv_l_n[0] * huv_l_n[2] - huv_r_n[0] * huv_r_n[2])
    };

    return F_eff;
}

template <typename CellType>
Solver<CellType>::~Solver() {};