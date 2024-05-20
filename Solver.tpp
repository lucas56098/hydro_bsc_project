#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
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

template <typename CellType>
Point Solver<CellType>::get_normal_vec(Point a, Point b) {

    double Delta_x = b.x - a.x;
    double Delta_y = b.y - a.y;

    double norm = sqrt(Delta_x*Delta_x + Delta_y*Delta_y);

    return Point(-Delta_y/norm, Delta_x/norm);
}


// requires Q_Cells
template <typename CellType>
void Solver<CellType>::diffusion_like_step(double dt) {
 
    if (grid->cells[0].cell_type != 1) {
        cerr << "diffusion_like_step called on Mesh with cell_type !=1, you must use Q_Cells" << endl;
        //exit(EXIT_FAILURE);
    } if (0.1 * dt >1){ //(/grid->q_cells[0].edges[0].length * grid->q_cells[0].edges[0].length) > 1) {
        cerr << "WARNING: c * dt/l^2 > 1, this time step likely wont be numerical stable, please reduce dt" << endl;
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
                //deltaQ += 0.1 * (grid->cells[i].edges[j].neighbour->getQ() - grid->cells[i].Q) * dt;
                deltaQ += 0.1 * (grid->cells[i].edges[j].neighbour->getQ() - grid->cells[i].Q);


            } else {

                deltaQ += 0;  // no flux through boundaries at the moment

            }

        }

        //cout << i << " : " << deltaQ << endl;

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

// requires Conway_Cells
template <typename CellType>
void Solver<CellType>::conway_step() {

    if (grid->cells[0].cell_type != 2) {
        cerr << "conway_step called on Mesh with cell_type !=2, you must use Conway_Cells" << endl;
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


        if (grid->cells[i].edges.size()>=2) {
            if (grid->cells[i].edges[1].is_boundary == false) {
                if (grid->cells[i].edges[1].neighbour->edges.size()>=1) {
                    if (grid->cells[i].edges[1].neighbour->edges[0].is_boundary == false) {
                        neighbour_count += grid->cells[i].edges[1].neighbour->edges[0].neighbour->getQ();
                    }
                }
            }
        }

        if (grid->cells[i].edges.size()>=2) {
            if (grid->cells[i].edges[1].is_boundary == false) {
                if (grid->cells[i].edges[1].neighbour->edges.size()>=3) {
                    if (grid->cells[i].edges[1].neighbour->edges[2].is_boundary == false) {
                        neighbour_count += grid->cells[i].edges[1].neighbour->edges[2].neighbour->getQ();
                    }
                }
            }
        }   

        if (grid->cells[i].edges.size()>=4) {
            if (grid->cells[i].edges[3].is_boundary == false) {
                if (grid->cells[i].edges[3].neighbour->edges.size()>=1) {
                    if (grid->cells[i].edges[3].neighbour->edges[0].is_boundary == false) {
                        neighbour_count += grid->cells[i].edges[3].neighbour->edges[0].neighbour->getQ();
                    }
                }
            }
        }  

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

// requires Q_Cells and 1D cartesian mesh, v_x and v_y >= 0! since upwind scheme doesnt work otherwhise
template <typename CellType>
void Solver<CellType>::advection_finite_difference_upwind_cartesian(double dt, Point v) {

    if (v.x < 0 || v.y < 0) {
        cerr << "v.x or v.y is < 0. This is not allowed for current version of cartesian upwind advection scheme" << endl;
        exit(EXIT_FAILURE);
    }

    vector<double> new_Qs;
    new_Qs.reserve(grid->cells.size());

    // calculate all the new Qs
    for (int i = 0; i<grid->cells.size(); i++) {

        double Qi_np1 = 0;
        double Qi_n = grid->cells[i].Q;
        double dx = grid->cells[i].edges[1].length;
        double dy = grid->cells[i].edges[2].length;


        if (dt * v.x >= dx || dt * v.y >= dy) {
            cerr << "the CFL timestep condition isnt fullfilled, please make sure that dt < dx/v. otherwise the simulation will show unstable behaviour" << endl;
            exit(EXIT_FAILURE);
        }

        // Qi_np1 = Qi_n - v * (dt/dx) * (Qi_n - 0);

        // double Qim1_n = grid->q_cells[i].edges[0].neighbour->getQ();

        // Qi_np1 = Qi_n - v * (dt/dx) * (Qi_n - Qim1_n);

        if (grid->cells[i].edges[0].is_boundary && grid->cells[i].edges[3].is_boundary) {
            Qi_np1 = Qi_n - v.x * (dt/dx) * (Qi_n - 0) - v.y * (dt/dy) * (Qi_n - 0);
        } else if (grid->cells[i].edges[0].is_boundary) {
            double Qim1y_n = grid->cells[i].edges[3].neighbour->getQ();
            Qi_np1 = Qi_n - v.x * (dt/dx) * (Qi_n - 0) - v.y * (dt/dy) * (Qi_n - Qim1y_n);
        } else if (grid->cells[i].edges[3].is_boundary) {
            double Qim1x_n = grid->cells[i].edges[0].neighbour->getQ();
            Qi_np1 = Qi_n - v.x * (dt/dx) * (Qi_n - Qim1x_n) - v.y * (dt/dy) * (Qi_n - 0);
        } else {
            double Qim1x_n = grid->cells[i].edges[0].neighbour->getQ();
            double Qim1y_n = grid->cells[i].edges[3].neighbour->getQ();
            Qi_np1 = Qi_n - v.x * (dt/dx) * (Qi_n - Qim1x_n) - v.y * (dt/dy) * (Qi_n - Qim1y_n);

        }

        new_Qs.push_back(Qi_np1);

    }

    // apply all the new Qs
    for (int i = 0; i<grid->cells.size(); i++) {

        grid->cells[i].Q = new_Qs[i];

    }

}


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

template <typename CellType>
Solver<CellType>::~Solver() {};