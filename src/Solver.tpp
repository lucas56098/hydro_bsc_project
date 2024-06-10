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
        double u_ip1_n = - grid->cells[i].u;  // will be later changed unless its a boundary
        double u_im1_n = - grid->cells[i].u;  // will be later changed unless its a boundary

        if (grid->cells[i].edges[0].is_boundary == false) {
            h_im1_n = grid->cells[i].edges[0].neighbour->get_h();
            u_im1_n = grid->cells[i].edges[0].neighbour->get_u();
        }
        if (grid->cells[i].edges[2].is_boundary == false) {
            h_ip1_n = grid->cells[i].edges[2].neighbour->get_h();
            u_ip1_n = grid->cells[i].edges[2].neighbour->get_u();
        }

        // new values using FV with Lax Friedrichs Flux estimation
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

// implementation of shallow water equation 2D cartesian using lax friedrichs flux
template<typename CellType>
void Solver<CellType>::shallow_water_2D_cartesian(double dt, int boundary_condition) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, SWE_Cell> == false) {
        cerr << "shallow water step called on Mesh with wrong cell type, you must use SWE_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    vector<double> h_new;
    vector<double> u_new;
    vector<double> v_new;
    h_new.reserve(grid->cells.size());
    u_new.reserve(grid->cells.size());
    v_new.reserve(grid->cells.size());

    // define d_x, d_y
    double dx = grid->cells[0].edges[1].length;
    double dy = grid->cells[0].edges[0].length;
    double g_half = 9.81/2.0;

    //go through each cell to calculate the new values
    for (int i = 0; i < grid->cells.size(); i++) {

        // Q_i_np1
        double h_i_np1 = 0;
        double u_i_np1 = 0;
        double v_i_np1 = 0;

        // Q_i_n
        double h_i_n = grid->cells[i].h;
        double u_i_n = grid->cells[i].u;
        double v_i_n = grid->cells[i].v;

        // conditions as if it were a boundary (will be changed again if not)
        // reflective boundary at the moment (will be generalized later on)
        // Q_im1_n (left cell)
        double h_im1_n = grid->cells[i].h;
        double u_im1_n = boundary_condition * grid->cells[i].u; //minus here
        double v_im1_n = grid->cells[i].v;

        // Q_ip1_n (right cell)
        double h_ip1_n = grid->cells[i].h;
        double u_ip1_n = boundary_condition * grid->cells[i].u; //minus here
        double v_ip1_n = grid->cells[i].v;

        // Q_ipN_n (upper cell)
        double h_ipN_n = grid->cells[i].h;
        double u_ipN_n = grid->cells[i].u;
        double v_ipN_n = boundary_condition * grid->cells[i].v; //minus here

        // Q_imN_n (lower cell)
        double h_imN_n = grid->cells[i].h;
        double u_imN_n = grid->cells[i].u;
        double v_imN_n = boundary_condition * grid->cells[i].v; //minus here


        if (grid->cells[i].edges[0].is_boundary == false) {
            // Q_im1_n (left cell)
            h_im1_n = grid->cells[i].edges[0].neighbour->get_h();
            u_im1_n = grid->cells[i].edges[0].neighbour->get_u();
            v_im1_n = grid->cells[i].edges[0].neighbour->get_v();
        }
        if (grid->cells[i].edges[1].is_boundary == false) {
            // Q_ipN_n (upper cell)
            h_ipN_n = grid->cells[i].edges[1].neighbour->get_h();
            u_ipN_n = grid->cells[i].edges[1].neighbour->get_u();
            v_ipN_n = grid->cells[i].edges[1].neighbour->get_v();

        }
        if (grid->cells[i].edges[2].is_boundary == false) {
            // Q_ip1_n (right cell)
            h_ip1_n = grid->cells[i].edges[2].neighbour->get_h();
            u_ip1_n = grid->cells[i].edges[2].neighbour->get_u();
            v_ip1_n = grid->cells[i].edges[2].neighbour->get_v();
        }
        if (grid->cells[i].edges[3].is_boundary == false) {
            // Q_imN_n (lower cell)
            h_imN_n = grid->cells[i].edges[3].neighbour->get_h();
            u_imN_n = grid->cells[i].edges[3].neighbour->get_u();
            v_imN_n = grid->cells[i].edges[3].neighbour->get_v();
        }

        // define F_im1, F_i, F_ip1
        // F_im1
        double F_im1_h = h_im1_n * u_im1_n;
        double F_im1_u = h_im1_n * u_im1_n * u_im1_n + g_half * h_im1_n * h_im1_n;
        double F_im1_v = h_im1_n * u_im1_n * v_im1_n;
        // F_i
        double F_i_h = h_i_n * u_i_n;
        double F_i_u = h_i_n * u_i_n * u_i_n + g_half * h_i_n * h_i_n;
        double F_i_v = h_i_n * u_i_n * v_i_n;
        // F_ip1
        double F_ip1_h = h_ip1_n * u_ip1_n;
        double F_ip1_u = h_ip1_n * u_ip1_n * u_ip1_n + g_half * h_ip1_n * h_ip1_n;
        double F_ip1_v = h_ip1_n * u_ip1_n * v_ip1_n;

        // define G_im1, G_i, G_ip1
        // G_im1
        double G_imN_h = h_imN_n * v_imN_n;
        double G_imN_u = h_imN_n * v_imN_n * u_imN_n;
        double G_imN_v = h_imN_n * v_imN_n * v_imN_n + g_half * h_imN_n * h_imN_n;
        // G_i
        double G_i_h = h_i_n * v_i_n;
        double G_i_u = h_i_n * v_i_n * u_i_n;
        double G_i_v = h_i_n * v_i_n * v_i_n + g_half * h_i_n * h_i_n;
        // G_ip1
        double G_ipN_h = h_ipN_n * v_ipN_n;
        double G_ipN_u = h_ipN_n * v_ipN_n * u_ipN_n;
        double G_ipN_v = h_ipN_n * v_ipN_n * v_ipN_n + g_half * h_ipN_n * h_ipN_n;

        // calculate F_m, F_p using lax friedrichs flux approximation
        // F_m
        double F_m_h = 0.5 * (F_im1_h + F_i_h) - (1.0/100.0) * ((dx)/(2.0*dt)) * (h_i_n - h_im1_n); 
        double F_m_u = 0.5 * (F_im1_u + F_i_u) - (1.0/100.0) * ((dx)/(2.0*dt)) * (h_i_n * u_i_n - h_im1_n * u_im1_n);
        double F_m_v = 0.5 * (F_im1_v + F_i_v) - (1.0/100.0) * ((dx)/(2.0*dt)) * (h_i_n * v_i_n - h_im1_n * v_im1_n);
        //double F_m_h = (0.5 * (F_im1_h + F_i_h) >= 0) ? F_im1_h : F_i_h;
        //double F_m_u = (0.5 * (F_im1_u + F_i_u) >= 0) ? F_im1_u : F_i_u;
        //double F_m_v = (0.5 * (F_im1_v + F_i_v) >= 0) ? F_im1_v : F_i_v;


        // F_p
        double F_p_h = 0.5 * (F_i_h + F_ip1_h) - (1.0/100.0) * ((dx)/(2.0*dt)) * (h_ip1_n - h_i_n);
        double F_p_u = 0.5 * (F_i_u + F_ip1_u) - (1.0/100.0) * ((dx)/(2.0*dt)) * (h_ip1_n * u_ip1_n - h_i_n * u_i_n);
        double F_p_v = 0.5 * (F_i_v + F_ip1_v) - (1.0/100.0) * ((dx)/(2.0*dt)) * (h_ip1_n * v_ip1_n - h_i_n * v_i_n);
        //double F_p_h = (0.5 * (F_ip1_h + F_i_h) <= 0) ? F_ip1_h : F_i_h;
        //double F_p_u = (0.5 * (F_ip1_u + F_i_u) <= 0) ? F_ip1_u : F_i_u;
        //double F_p_v = (0.5 * (F_ip1_v + F_i_v) <= 0) ? F_ip1_v : F_i_v;

        // calculate G_m, G_p using lax friedrichs flux approximation
        // G_m
        double G_m_h = 0.5 * (G_imN_h + G_i_h) - (1.0/100.0) * ((dy)/(2.0*dt)) * (h_i_n - h_imN_n);
        double G_m_u = 0.5 * (G_imN_u + G_i_u) - (1.0/100.0) * ((dy)/(2.0*dt)) * (h_i_n * u_i_n - h_imN_n * u_imN_n);
        double G_m_v = 0.5 * (G_imN_v + G_i_v) - (1.0/100.0) * ((dy)/(2.0*dt)) * (h_i_n * v_i_n - h_imN_n * v_imN_n);
        //double G_m_h = (0.5 * (G_imN_h + G_i_h) >= 0) ? G_imN_h : G_i_h;
        //double G_m_u = (0.5 * (G_imN_u + G_i_u) >= 0) ? G_imN_u : G_i_u;
        //double G_m_v = (0.5 * (G_imN_v + G_i_v) >= 0) ? G_imN_v : G_i_v;


        // G_p
        double G_p_h = 0.5 * (G_i_h + G_ipN_h) - (1.0/100.0) * ((dy)/(2.0*dt)) * (h_ipN_n - h_i_n);
        double G_p_u = 0.5 * (G_i_u + G_ipN_u) - (1.0/100.0) * ((dy)/(2.0*dt)) * (h_ipN_n * u_ipN_n - h_i_n * u_i_n);
        double G_p_v = 0.5 * (G_i_v + G_ipN_v) - (1.0/100.0) * ((dy)/(2.0*dt)) * (h_ipN_n * v_ipN_n - h_i_n * v_i_n);
        //double G_p_h = (0.5 * (G_imN_h + G_i_h) <= 0) ? G_ipN_h : G_i_h;
        //double G_p_u = (0.5 * (G_imN_u + G_i_u) <= 0) ? G_ipN_u : G_i_u;
        //double G_p_v = (0.5 * (G_imN_v + G_i_v) <= 0) ? G_ipN_v : G_i_v;

        // calculate new state vector Q_i_np1 using fluxes
        h_i_np1 = h_i_n + ((dt)/(dx)) * (F_m_h - F_p_h) + ((dt)/(dy)) * (G_m_h - G_p_h);
        u_i_np1 = (h_i_n * u_i_n + ((dt)/(dx)) * (F_m_u - F_p_u) + ((dt)/(dy)) * (G_m_u - G_p_u))/h_i_np1;
        v_i_np1 = (h_i_n * v_i_n + ((dt)/(dx)) * (F_m_v - F_p_v) + ((dt)/(dy)) * (G_m_v - G_p_v))/h_i_np1;

        // append new state vector to vectors
        h_new.push_back(h_i_np1);
        u_new.push_back(u_i_np1);
        v_new.push_back(v_i_np1);

    }

    // then apply all the changes
    for (int i = 0; i < grid->cells.size(); i++) {
        grid->cells[i].h = h_new[i];
        grid->cells[i].u = u_new[i];
        grid->cells[i].v = v_new[i];
    }

}


// DOESNT WORK!!! -- STILL CONSTRUCTION SITE -- implementation of shallow water equation 2D voronoi
template<typename CellType>
void Solver<CellType>::shallow_water_2D(double dt) {

    // make sure that the cell type is correct
    if constexpr (is_same_v<CellType, SWE_Cell> == false) {
        cerr << "shallow water step called on Mesh with wrong cell type, you must use SWE_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    vector<double> h_new;
    vector<double> u_new;
    vector<double> v_new;
    h_new.reserve(grid->cells.size());
    u_new.reserve(grid->cells.size());
    v_new.reserve(grid->cells.size());

    //go through each cell to calculate the new values
    for (int i = 0; i<grid->cells.size(); i++) {

        double A = grid->cells[i].volume;

        double h_i_np1 = 0;
        double u_i_np1 = 0;
        double v_i_np1 = 0;

        double h_i_n = grid->cells[i].h;
        double u_i_n = grid->cells[i].u;
        double v_i_n = grid->cells[i].v;

        double total_flux_h = 0;
        double total_flux_u = 0;
        double total_flux_v = 0;

        // calculate total_flux by summing over edge fluxed * edge_length
        for (int j = 0; j < grid->cells[i].edges.size(); j++) {
            
            // make sure that no flux through boundaries 
            if (grid->cells[i].edges[j].is_boundary == false) {

                Point n = get_normal_vec(grid->cells[i].edges[j].a, grid->cells[i].edges[j].b);
                double l_i_j = grid->cells[i].edges[j].length;
                double g = 9.81;

                // values needed for f and g calculations
                // values of other cell
                double h_j_n = grid->cells[i].edges[j].neighbour->get_h();
                double u_j_n = grid->cells[i].edges[j].neighbour->get_u();
                double v_j_n = grid->cells[i].edges[j].neighbour->get_v();

                // f_i
                double f_i_h = h_i_n * u_i_n;
                double f_i_u = h_i_n * u_i_n * u_i_n + (g/2.0) * h_i_n * h_i_n;
                double f_i_v = h_i_n * v_i_n * u_i_n;
                
                // f_j
                double f_j_h = h_j_n * u_j_n;
                double f_j_u = h_j_n * u_j_n * u_j_n + (g/2.0) * h_j_n * h_j_n;
                double f_j_v = h_j_n * v_j_n * u_j_n;
                
                // g_i
                double g_i_h = h_i_n * v_i_n;
                double g_i_u = h_i_n * v_i_n * u_i_n;
                double g_i_v = h_i_n * v_i_n * v_i_n + (g/2.0) * h_i_n * h_i_n;
                
                // g_j
                double g_j_h = h_j_n * v_j_n;
                double g_j_u = h_j_n * v_j_n * u_j_n;
                double g_j_v = h_j_n * v_j_n * v_j_n + (g/2.0) * h_j_n * h_j_n;

                // distance_i_j
                Point seed_i = grid->cells[i].seed;
                Point seed_j = grid->cells[i].edges[j].neighbour->seed;
                double d_i_j = sqrt((seed_i.x - seed_j.x)*(seed_i.x - seed_j.x) + (seed_i.y - seed_j.y)*(seed_i.y - seed_j.y));


                // calculate f using lax_friedrichs_scheme
                double f_i_j_h = (0.5) * (f_i_h + f_j_h) - (d_i_j/(2.0*dt)) * (h_j_n - h_i_n);
                double f_i_j_u = (0.5) * (f_i_u + f_j_u) - (d_i_j/(2.0*dt)) * (h_j_n * u_j_n - h_i_n * u_i_n);
                double f_i_j_v = (0.5) * (f_i_v + f_j_v) - (d_i_j/(2.0*dt)) * (h_j_n * v_j_n - h_i_n * v_i_n);
            
                // calculate g using lax_friedrichs_scheme
                double g_i_j_h = (0.5) * (g_i_h + g_j_h) - (d_i_j/(2.0*dt)) * (h_j_n - h_i_n);
                double g_i_j_u = (0.5) * (g_i_u + g_j_u) - (d_i_j/(2.0*dt)) * (h_j_n * u_j_n - h_i_n * u_i_n);
                double g_i_j_v = (0.5) * (g_i_v + g_j_v) - (d_i_j/(2.0*dt)) * (h_j_n * v_j_n - h_i_n * v_i_n);

                // add flux (corrected for direction) times lenght to total_flux
                total_flux_h += (f_i_j_h*n.x + g_i_j_h*n.y) * l_i_j;
                total_flux_u += (f_i_j_u*n.x + g_i_j_u*n.y) * l_i_j;
                total_flux_v += (f_i_j_v*n.x + g_i_j_v*n.y) * l_i_j;

            } 

        }

        // calculate new values
        h_i_np1 = h_i_n + (dt/A) * total_flux_h;
        u_i_np1 = (h_i_n*u_i_n + (dt/A) * total_flux_u)/h_i_np1;
        v_i_np1 = (h_i_n*v_i_n + (dt/A) * total_flux_v)/h_i_np1;

        // store new values
        h_new.push_back(h_i_np1);
        u_new.push_back(u_i_np1);
        v_new.push_back(v_i_np1);

    }

    // then apply all the changes
    for (int i = 0; i<grid->cells.size(); i++) {
        grid->cells[i].h = h_new[i];
        grid->cells[i].u = u_new[i];
        grid->cells[i].v = v_new[i];
    }

}


template <typename CellType>
Solver<CellType>::~Solver() {};