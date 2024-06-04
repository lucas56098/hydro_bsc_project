#include <iostream>
#include <vector>
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


// Evolution Step for Conway's Game of Life, requires conway_cell
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
template<typename CellType>
void Solver<CellType>::shallow_water_1D_cartesian(double dt) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, SWE_Cell> == false) {
        cerr << "shallow water step called on Mesh with wrong cell type, you must use SWE_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    // first calculate all new values and store in vectors

    vector<double> h_new;
    vector<double> u_new;
    h_new.reserve(grid->cells.size());
    u_new.reserve(grid->cells.size());

    // loop through all cells to calcualte values
    for (int i = 0; i<grid->cells.size(); i++) {

        // get everything we need
        double l_x = grid->cells[i].edges[1].length;
        double g = 9.81;
        double h_i_n = grid->cells[i].h;
        double h_ip1_n = grid->cells[i].h; // will be later changed unless its a boundary
        double h_im1_n = grid->cells[i].h; // will be later changed unless its a boundary
        double u_i_n = grid->cells[i].u;
        double u_ip1_n = 0;                // will be later changed unless its a boundary
        double u_im1_n = 0;                // will be later changed unless its a boundary

        if (grid->cells[i].edges[0].is_boundary == false) {
            h_im1_n = grid->cells[i].edges[0].neighbour->get_h();
            u_im1_n = grid->cells[i].edges[0].neighbour->get_u();
        }
        if (grid->cells[i].edges[2].is_boundary == false) {
            h_ip1_n = grid->cells[i].edges[2].neighbour->get_h();
            u_ip1_n = grid->cells[i].edges[2].neighbour->get_u();
        }

        // new values using FV with Lax Friedrichs Flux estimation
        //double h_i_np1 = h_i_n - ((dt)/(2.0*l_x)) * ( ((h_im1_n*u_im1_n - h_i_n*u_i_n) - ((l_x)/(dt))*(h_i_n - h_im1_n)) - ((h_ip1_n*u_ip1_n - h_i_n*u_i_n) - ((l_x)/(dt))*(h_i_n - h_ip1_n)) );
        //double h_i_np1 = h_i_n - ((dt)/(2.0*l_x)) * ( ((h_im1_n*u_im1_n - h_i_n*u_i_n) + ((l_x)/(dt))*(h_i_n - h_im1_n)) - ((h_ip1_n*u_ip1_n - h_i_n*u_i_n) - ((l_x)/(dt))*(h_i_n - h_ip1_n)) );
        //double u_i_np1 = ((h_i_n * u_i_n)/(h_i_np1)) - ((dt)/(2.0*l_x*h_i_np1)) *  ( ((h_im1_n * u_im1_n * u_im1_n + (g/2.0)* h_im1_n * h_im1_n) - (h_i_n * u_i_n * u_i_n + (g/2.0)* h_i_n * h_i_n) - ((l_x)/(dt))*(h_i_n * u_i_n - h_im1_n * u_im1_n)) 
        //                                                                           - ((h_ip1_n * u_ip1_n * u_ip1_n + (g/2.0)* h_ip1_n * h_ip1_n) - (h_i_n * u_i_n * u_i_n + (g/2.0)* h_i_n * h_i_n) - ((l_x)/(dt))*(h_i_n * u_i_n - h_ip1_n * u_ip1_n)) );
        //double u_i_np1 = ((h_i_n * u_i_n)/(h_i_np1)) - ((dt)/(2.0*l_x*h_i_np1)) *  ( ((h_im1_n * u_im1_n * u_im1_n + (g/2.0)* h_im1_n * h_im1_n) - (h_i_n * u_i_n * u_i_n + (g/2.0)* h_i_n * h_i_n) + ((l_x)/(dt))*(h_i_n * u_i_n - h_im1_n * u_im1_n)) 
        //                                                                           - ((h_ip1_n * u_ip1_n * u_ip1_n + (g/2.0)* h_ip1_n * h_ip1_n) - (h_i_n * u_i_n * u_i_n + (g/2.0)* h_i_n * h_i_n) - ((l_x)/(dt))*(h_i_n * u_i_n - h_ip1_n * u_ip1_n)) );
    //    double h_i_np1 = h_i_n - ((dt)/(2.0*l_x)) * ( (-1*(h_im1_n*u_im1_n - h_i_n*u_i_n) + ((l_x)/(dt))*(h_i_n - h_im1_n)) - ( -1*(h_ip1_n*u_ip1_n - h_i_n*u_i_n) - ((l_x)/(dt))*(h_i_n - h_ip1_n)) );
  //      double u_i_np1 = ((h_i_n * u_i_n)/(h_i_np1)) + ((dt)/(2.0*l_x*h_i_np1)) *  ( ((h_im1_n * u_im1_n * u_im1_n + (g/2.0)* h_im1_n * h_im1_n) - (h_i_n * u_i_n * u_i_n + (g/2.0)* h_i_n * h_i_n) - ((l_x)/(dt))*(h_i_n * u_i_n - h_im1_n * u_im1_n)) 
//                                                                                   - ((h_ip1_n * u_ip1_n * u_ip1_n + (g/2.0)* h_ip1_n * h_ip1_n) - (h_i_n * u_i_n * u_i_n + (g/2.0)* h_i_n * h_i_n) + ((l_x)/(dt))*(h_i_n * u_i_n - h_ip1_n * u_ip1_n)) );

        double h_i_np1 = h_i_n - (dt/(2*l_x)) * ( (h_ip1_n*u_ip1_n - h_im1_n*u_im1_n) + (l_x/dt)*(2*h_i_n - h_im1_n - h_ip1_n) );
        double u_i_np1 = ((h_i_n*u_i_n)/h_i_np1) - (dt/(2*l_x*h_i_np1)) * ( (h_ip1_n*u_ip1_n*u_ip1_n + (g/2.0)*h_ip1_n*h_ip1_n) - ((h_im1_n*u_im1_n*u_im1_n + (g/2.0)*h_im1_n*h_im1_n)) + (l_x/dt)*(2*h_i_n*u_i_n - h_ip1_n*u_ip1_n - h_im1_n*u_im1_n) );

        h_new.push_back(h_i_np1);
        u_new.push_back(u_i_np1);
    }
    
    // then apply all the changes
    for (int i = 0; i<grid->cells.size(); i++) {
        grid->cells[i].h = h_new[i];
        grid->cells[i].u = u_new[i];
    }

}

template <typename CellType>
Solver<CellType>::~Solver() {};