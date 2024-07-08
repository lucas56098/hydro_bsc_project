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

// function to rotate a swe vector where [1] and [2] are the velocities to rotate and [0] stays const.
template <typename CellType>
array<double, 3> Solver<CellType>::rotate2Dswe(array<double, 3> vec, double angle) {

    array<double, 3> rot_vec = {
                                    vec[0],
                                    cos(angle) * vec[1] - sin(angle) * vec[2],
                                    sin(angle) * vec[1] + cos(angle) * vec[2]
                               };

    return rot_vec;

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
void Solver<CellType>::shallow_water(double dt, int boundary_cond, double numerical_diffusion_coeff, int sim_order, double g) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, SWE_Cell> == false) {
        cerr << "shallow water step called on Mesh with wrong cell type, you must use SWE_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    vector<array<double, 3>> Q_new;
    Q_new.reserve(grid->cells.size());

    vector<array<array<double, 3>, 2>> gradients;
    if (sim_order == 2) {
        // calculate gradients
        gradients.reserve(grid->cells.size());

        // go through each cell to calculate gradient
        for (int i = 0; i < grid->cells.size(); i++) {
            gradients.push_back(calc_swe_gradients(i, boundary_cond, g));
        }
    }


    //go through each cell to calculate the new values
    for (int i = 0; i<grid->cells.size(); i++) {

        double A = grid->cells[i].volume;

        array<double, 3> huv_i_np1 = {0, 0, 0};
        array<double, 3> huv_i_n = {grid->cells[i].h, grid->cells[i].u, grid->cells[i].v};
        array<double, 3> total_flux = {0, 0, 0};

        // calculate total_flux by summing over edge fluxed * edge_length
        for (int j = 0; j < grid->cells[i].edges.size(); j++) {
            
            Point n = get_normal_vec(grid->cells[i].edges[j].a, grid->cells[i].edges[j].b);
            double l_i_j = grid->cells[i].edges[j].length;

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

            array<double, 3> F_ij = {0.0, 0.0, 0.0};

            if (sim_order == 1) {
                // get Flux using a HLL solver (first order)
                F_ij = hll_solver_swe_2D(huv_i_n, huv_j_n, g, n, numerical_diffusion_coeff);
            } else if (sim_order == 2) {
                // extrapolate values
                // extrapolate U_i
                array<double, 3> huv_i_ext = linear_extrapolate_swe(huv_i_n, gradients[i], dt, i, j, g);

                // if j not boundary extrapolate U_j
                array<double, 3> huv_j_ext = {0.0, 0.0, 0.0};
                if (grid->cells[i].edges[j].is_boundary == false) {
                    // extrapolate U_j

                    // index_j is the index of the Cell_j in the cell list
                    int index_j = grid->cells[i].edges[j].neighbour->index;
                    // index_i is the neighbour index of Cell_i inside of Cell_j
                    int index_i;

                    // get index_i by looping through all edges of Cell_j
                    for (int k = 0; k < grid->cells[index_j].edges.size(); k++) {
                        if (grid->cells[index_j].edges[k].is_boundary == false) {
                            if (grid->cells[index_j].edges[k].neighbour->index == i) {
                                index_i = k;
                            }
                        }
                    }
                    
                    huv_j_ext = linear_extrapolate_swe(huv_j_n, gradients[index_j], dt, index_j, index_i, g);
                // else U_j according to boundary conditions
                } else {
                    huv_j_ext = {huv_i_ext[0], boundary_cond*huv_i_ext[1], boundary_cond*huv_i_ext[2]};

                    if (boundary_cond == -1) {
                        huv_j_ext[1] = huv_i_ext[1] - 2*(huv_i_ext[1]*n.x + huv_i_ext[2]*n.y)*n.x;
                        huv_j_ext[2] = huv_i_ext[2] - 2*(huv_i_ext[1]*n.x + huv_i_ext[2]*n.y)*n.y;
                    }
                }

                // do HLL with extrapolated values to get Flux
                F_ij = hll_solver_swe_2D(huv_i_ext, huv_j_ext, g, n, numerical_diffusion_coeff);
            } else {
                cerr << "no such simulation order allowed!, choose sim_order = 1 or = 2" << endl;
            }

            // add to total flux * length
            total_flux[0] += F_ij[0] * l_i_j;
            total_flux[1] += F_ij[1] * l_i_j;
            total_flux[2] += F_ij[2] * l_i_j;
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


// HLL Solver for Shallow Water Equation 2D unstructured
template <typename CellType>
array<double, 3> Solver<CellType>::hll_solver_swe_2D(array<double, 3> huv_i_n, array<double, 3> huv_j_n, double g, Point n, double add_diffusion_coeff) {
    
    // rotate huv_i, huv_j into n-frame
    double theta = atan2(n.y, n.x);
    array<double, 3> huv_l = rotate2Dswe(huv_i_n, theta);
    array<double, 3> huv_r = rotate2Dswe(huv_j_n, theta);

    // using huv_l, huv_r (the rotated ones) now calc F_HLLn using HLL solver
    // calc f_l and f_r
    array<double, 3> f_l = {
                                huv_l[0] * huv_l[1],
                                huv_l[0]*huv_l[1]*huv_l[1] + (g/2.0) * huv_l[0] * huv_l[0],
                                huv_l[0] * huv_l[1] * huv_l[2]
                           };

    array<double, 3> f_r = {
                                huv_r[0] * huv_r[1],
                                huv_r[0]*huv_r[1]*huv_r[1] + (g/2.0) * huv_r[0] * huv_r[0],
                                huv_r[0] * huv_r[1] * huv_r[2]
                           };

    // calculate wave speeds
    double SL = min(huv_l[1] - sqrt(g*huv_l[0]), huv_r[1] - sqrt(g*huv_r[0]));
    double SR = max(huv_r[1] + sqrt(g*huv_r[0]), huv_l[1] + sqrt(g*huv_l[0]));

    // HLL solver for F
    array<double, 3> F_HLL_n;
    if (SL >= 0) {
        F_HLL_n = f_l;
    } else if (SL < 0 && SR > 0) {
        F_HLL_n[0] = (SR * f_l[0] - SL * f_r[0] + SL * SR * (huv_r[0] - huv_l[0])) / (SR - SL);
        F_HLL_n[1] = (SR * f_l[1] - SL * f_r[1] + SL * SR * (huv_r[0] * huv_r[1] - huv_l[0] * huv_l[1])) / (SR - SL);
        F_HLL_n[2] = (SR * f_l[2] - SL * f_r[2] + SL * SR * (huv_r[0] * huv_r[2] - huv_l[0] * huv_l[2])) / (SR - SL);
    } else if (SR <= 0) {
        F_HLL_n = f_r;
    }

    // get F_HLL by rotating F_HLLn back into lab frame
    array<double, 3> F_HLL = rotate2Dswe(F_HLL_n, -theta);
    
    // optional additional diffusion
    F_HLL[0] += - add_diffusion_coeff * (huv_l[0] - huv_r[0]);
    F_HLL[1] += - add_diffusion_coeff * (huv_l[0] * huv_l[1] - huv_r[0] * huv_r[1]);
    F_HLL[2] += - add_diffusion_coeff * (huv_l[0] * huv_l[2] - huv_r[0] * huv_r[2]);

    return F_HLL;
}


// GRADIENTS AND EXTRAPOLATION --------------------------------------------------------------------
template <typename CellType>
array<array<double, 3>, 2> Solver<CellType>::calc_swe_gradients(int i, int boundary_cond, double g) {

    double A = grid->cells[i].volume;
    array<double, 3> huv_i_n = {grid->cells[i].h, grid->cells[i].u, grid->cells[i].v};

    array<array<double, 3>, 2> gradientU;
    gradientU[0] = {0.0, 0.0, 0.0};
    gradientU[1] = {0.0, 0.0, 0.0};

    // go through all edges
    for (int j = 0; j < grid->cells[i].edges.size(); j++) {
            
        Point n = get_normal_vec(grid->cells[i].edges[j].a, grid->cells[i].edges[j].b);
        double l_i_j = grid->cells[i].edges[j].length;

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

        // get c_ij and |r_i_j|
        Point f_mid = Point((grid->cells[i].edges[j].a.x + grid->cells[i].edges[j].b.x)/2.0, (grid->cells[i].edges[j].a.y + grid->cells[i].edges[j].b.y)/2.0);
        Point r_i = grid->cells[i].seed;
        Point mid = Point(f_mid.x - r_i.x, f_mid.y - r_i.y);


        double r_ij = 2 * (mid.x * n.x + mid.y * n.y);
        Point c_ij = Point(
                            mid.x - 0.5*r_ij*n.x,
                            mid.y - 0.5*r_ij*n.y
                          );


        // add to sum of gradient
        gradientU[0][0] += (l_i_j/A) * (((huv_i_n[0] + huv_j_n[0])/2.0)*n.x + (huv_j_n[0] - huv_i_n[0]) * (c_ij.x/r_ij));
        gradientU[0][1] += (l_i_j/A) * (((huv_i_n[0]*huv_i_n[1] + huv_j_n[0]*huv_j_n[1])/2.0)*n.x + (huv_j_n[0]*huv_j_n[1] - huv_i_n[0]*huv_i_n[1]) * (c_ij.x/r_ij));
        gradientU[0][2] += (l_i_j/A) * (((huv_i_n[0]*huv_i_n[2] + huv_j_n[0]*huv_j_n[2])/2.0)*n.x + (huv_j_n[0]*huv_j_n[2] - huv_i_n[0]*huv_i_n[2]) * (c_ij.x/r_ij));

        gradientU[1][0] += (l_i_j/A) * (((huv_i_n[0] + huv_j_n[0])/2.0)*n.y + (huv_j_n[0] - huv_i_n[0]) * (c_ij.y/r_ij));
        gradientU[1][1] += (l_i_j/A) * (((huv_i_n[0]*huv_i_n[1] + huv_j_n[0]*huv_j_n[1])/2.0)*n.y + (huv_j_n[0]*huv_j_n[1] - huv_i_n[0]*huv_i_n[1]) * (c_ij.y/r_ij));
        gradientU[1][2] += (l_i_j/A) * (((huv_i_n[0]*huv_i_n[2] + huv_j_n[0]*huv_j_n[2])/2.0)*n.y + (huv_j_n[0]*huv_j_n[2] - huv_i_n[0]*huv_i_n[2]) * (c_ij.y/r_ij));

    }
    
    // slope limiting
    // initalize a_i
    array<double, 3> a_i = {1.0, 1.0, 1.0};
    array<double, 3> U_i = {
                            huv_i_n[0],
                            huv_i_n[0] * huv_i_n[1],
                            huv_i_n[0] * huv_i_n[2]
                         };

    // go through all edges to calculate correct a_i
    for (int j=0; j<grid->cells[i].edges.size(); j++) {
        // calculate delta_ij
        array<double, 3> delta_ij;
        Point f_mid = Point((grid->cells[i].edges[j].a.x + grid->cells[i].edges[j].b.x)/2.0, (grid->cells[i].edges[j].a.y + grid->cells[i].edges[j].b.y)/2.0);
        Point centroid = grid->cells[i].centroid;
        Point n = get_normal_vec(grid->cells[i].edges[j].a, grid->cells[i].edges[j].b);

        for (int k = 0; k < 3; k++) {
            delta_ij[k] = gradientU[0][k] * (f_mid.x - centroid.x) + gradientU[1][k] * (f_mid.y - centroid.y);
        }

        // calcualte psi_ij
        array<double, 3> psi_ij;

        for (int k = 0; k < 3; k++) {
            if (delta_ij[k] > 0) {
                // (psi_max - psi_i)/delta_psi_ij
                double psi_max = U_i[k];
                // loop through all edges to maximise psi
                for (int l = 0; l < grid->cells[i].edges.size(); l++) {
                    
                    // values of other cell (initalized with sink boundaries for boundary cond = 1)
                    array<double, 3> huv_j_n = {huv_i_n[0], boundary_cond * huv_i_n[1],  boundary_cond * huv_i_n[2]};
                    // if reflective boundary change velocities accordingly
                    if (boundary_cond == -1) {
                        huv_j_n[1] = huv_i_n[1] - 2*(huv_i_n[1]*n.x + huv_i_n[2]*n.y)*n.x;
                        huv_j_n[2] = huv_i_n[2] - 2*(huv_i_n[1]*n.x + huv_i_n[2]*n.y)*n.y;
                    }
                    // get actual huv values if its not a boundary
                    if (grid->cells[i].edges[l].is_boundary == false) {
                        huv_j_n[0] = grid->cells[i].edges[l].neighbour->get_h();
                        huv_j_n[1] = grid->cells[i].edges[l].neighbour->get_u();
                        huv_j_n[2] = grid->cells[i].edges[l].neighbour->get_v();
                    }

                    array<double, 3> U_j = {
                                                huv_j_n[0],
                                                huv_j_n[0] * huv_j_n[1],
                                                huv_j_n[0] * huv_j_n[2],
                                           };

                    // actually maximise
                    if (U_j[k] > psi_max) {psi_max = U_j[k];}

                } 

                psi_ij[k] = (psi_max - U_i[k])/delta_ij[k];


            } else if (delta_ij[k] < 0) {
                // (psi_min - psi_i)/delta_psi_ij
                double psi_min = U_i[k];

                // loop through all edges to minimze psi
                for (int l = 0; l < grid->cells[i].edges.size(); l++) {
                    
                    // values of other cell (initalized with sink boundaries for boundary cond = 1)
                    array<double, 3> huv_j_n = {huv_i_n[0], boundary_cond * huv_i_n[1],  boundary_cond * huv_i_n[2]};
                    // if reflective boundary change velocities accordingly
                    if (boundary_cond == -1) {
                        huv_j_n[1] = huv_i_n[1] - 2*(huv_i_n[1]*n.x + huv_i_n[2]*n.y)*n.x;
                        huv_j_n[2] = huv_i_n[2] - 2*(huv_i_n[1]*n.x + huv_i_n[2]*n.y)*n.y;
                    }
                    // get actual huv values if its not a boundary
                    if (grid->cells[i].edges[l].is_boundary == false) {
                        huv_j_n[0] = grid->cells[i].edges[l].neighbour->get_h();
                        huv_j_n[1] = grid->cells[i].edges[l].neighbour->get_u();
                        huv_j_n[2] = grid->cells[i].edges[l].neighbour->get_v();
                    }

                    array<double, 3> U_j = {
                                                huv_j_n[0],
                                                huv_j_n[0] * huv_j_n[1],
                                                huv_j_n[0] * huv_j_n[2],
                                           };

                    // actually minimize
                    if (U_j[k] < psi_min) {psi_min = U_j[k];}

                } 

                psi_ij[k] = (psi_min - U_i[k])/delta_ij[k];

            } else if (delta_ij[k] == 0) {
                psi_ij[k] = 1;
            }
        }

        // update a_i (e.g. minimizing a_i = min_{over j}(1, psi_ij))
        for (int k = 0; k < 3; k++) {
            if (psi_ij[k] < a_i[k]) {a_i[k] = psi_ij[k];}
        }
    }
    
    // apply slope limiters
    gradientU[0][0] = a_i [0] * gradientU[0][0];
    gradientU[1][0] = a_i [0] * gradientU[1][0];
    gradientU[0][1] = a_i [1] * gradientU[0][1];
    gradientU[1][1] = a_i [1] * gradientU[1][1];
    gradientU[0][2] = a_i [2] * gradientU[0][2];
    gradientU[1][2] = a_i [2] * gradientU[1][2];
    
    // return gradient which now can be used to extrapolate the values
    return gradientU;
}

template <typename CellType>
array<double, 3> Solver<CellType>::linear_extrapolate_swe(array<double, 3> huv_i_n, array<array<double, 3>, 2> gradient, double dt, int i, int j, double g) {

    // U_i
    double h = huv_i_n[0];
    double u = huv_i_n[1];
    double v = huv_i_n[2];
    array<double, 3> U_i = {
                                h,
                                h*u,
                                h*v
                           };

    // d
    Point f_mid = Point((grid->cells[i].edges[j].a.x + grid->cells[i].edges[j].b.x)/2.0, (grid->cells[i].edges[j].a.y + grid->cells[i].edges[j].b.y)/2.0);
    Point centroid = grid->cells[i].centroid;
    Point d = Point(f_mid.x - centroid.x, f_mid.y - centroid.y);

    // U_i_ext
    array<double, 3> U_i_ext = {
                                    U_i[0] + (d.x*gradient[0][0] - (dt/2.0)*(gradient[0][1])) + (d.y*gradient[1][0] - (dt/2.0)*(gradient[1][2])),
                                    U_i[1] + (d.x*gradient[0][1] - (dt/2.0)*((g*h - u*u)*gradient[0][0] + (2*u)*gradient[0][1])) + (d.y*gradient[1][1] - (dt/2.0)*((-1*u*v)*gradient[1][0] + (v)*gradient[1][1] + (u)*gradient[1][2])),
                                    U_i[2] + (d.x*gradient[0][2] - (dt/2.0)*((-1*u*v)*gradient[0][0] + (v)*gradient[0][1] + (u)*gradient[0][2])) + (d.y*gradient[1][0] - (dt/2.0)*((g*h - v*v)*gradient[1][0] + (2*v)* gradient[1][2]))
                               };

    // huv_i_ext;
    array<double, 3> huv_i_ext;
    huv_i_ext[0] = U_i_ext[0];
    huv_i_ext[1] = U_i_ext[1]/U_i_ext[0];
    huv_i_ext[2] = U_i_ext[2]/U_i_ext[0];

    return huv_i_ext;
}



template <typename CellType>
Solver<CellType>::~Solver() {};