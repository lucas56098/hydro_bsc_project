#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "cell_types/Q_Cell.h"
#include "Mesh.h"
#include "Solver.h"

Solver::Solver() {}

Solver::Solver(Mesh* gridin) {

    grid = gridin;

}

Point Solver::get_normal_vec(Point a, Point b) {

    double Delta_x = b.x - a.x;
    double Delta_y = b.y - a.y;

    double norm = sqrt(Delta_x*Delta_x + Delta_y*Delta_y);

    return Point(-Delta_y/norm, Delta_x/norm);
}


// requires Q_Cells
void Solver::diffusion_like_step(double dt) {
 

    if (grid->cell_type != 1) {
        cerr << "diffusion_like_step called on Mesh with cell_type !=1, you must use Q_Cells" << endl;
        exit(EXIT_FAILURE);
    } if (0.1 * dt >1){ //(/grid->q_cells[0].edges[0].length * grid->q_cells[0].edges[0].length) > 1) {
        cerr << "WARNING: c * dt/l^2 > 1, this time step likely wont be numerical stable, please reduce dt" << endl;
        exit(EXIT_FAILURE);
    }

    vector<double> Qs_new;

    // go through all cells and calculate total in/outflow
    for (int i = 0; i < grid->q_cells.size(); i++) {

        double deltaQ = 0;

        // calculate flux for each face and add to total in/outflow
        for (int j = 0; j< grid->q_cells[i].edges.size(); j++) {

            // calculate flux for normal faces and boundary faces
            if (grid->q_cells[i].edges[j].is_boundary == false) {


                // diffusion like flux
                deltaQ += 0.1 * (grid->q_cells[i].edges[j].neighbour->getQ() - grid->q_cells[i].Q) * dt;


            } else {

                deltaQ += 0;  // no flux through boundaries at the moment

            }

        }


        // new Q = old Q + delta Q
        Qs_new.push_back(grid->q_cells[i].Q + deltaQ);

    }

    // after calculating now exchange Q_new for Q_old
    for (int i = 0; i < grid->q_cells.size(); i++) {
        grid->q_cells[i].Q = Qs_new[i];
    }

}

// requires Conway_Cells
void Solver::conway_step() {

    if (grid->cell_type != 2) {
        cerr << "conway_step called on Mesh with cell_type !=2, you must use Conway_Cells" << endl;
        exit(EXIT_FAILURE);
    }

    vector<int> Qs_new;

    // go through all cells and decide whether they will live or die
    for (int i = 0; i < grid->conway_cells.size(); i++) {

        int neighbour_count = 0;

        // count the living neighbours
        for (int j = 0; j< grid->conway_cells[i].edges.size(); j++) {

            if (grid->conway_cells[i].edges[j].is_boundary == false) {

                neighbour_count += grid->conway_cells[i].edges[j].neighbour->getQ();

            }

        }


        if (grid->conway_cells[i].edges.size()>=2) {
            if (grid->conway_cells[i].edges[1].is_boundary == false) {
                if (grid->conway_cells[i].edges[1].neighbour->edges.size()>=1) {
                    if (grid->conway_cells[i].edges[1].neighbour->edges[0].is_boundary == false) {
                        neighbour_count += grid->conway_cells[i].edges[1].neighbour->edges[0].neighbour->getQ();
                    }
                }
            }
        }

        if (grid->conway_cells[i].edges.size()>=2) {
            if (grid->conway_cells[i].edges[1].is_boundary == false) {
                if (grid->conway_cells[i].edges[1].neighbour->edges.size()>=3) {
                    if (grid->conway_cells[i].edges[1].neighbour->edges[2].is_boundary == false) {
                        neighbour_count += grid->conway_cells[i].edges[1].neighbour->edges[2].neighbour->getQ();
                    }
                }
            }
        }   

        if (grid->conway_cells[i].edges.size()>=4) {
            if (grid->conway_cells[i].edges[3].is_boundary == false) {
                if (grid->conway_cells[i].edges[3].neighbour->edges.size()>=1) {
                    if (grid->conway_cells[i].edges[3].neighbour->edges[0].is_boundary == false) {
                        neighbour_count += grid->conway_cells[i].edges[3].neighbour->edges[0].neighbour->getQ();
                    }
                }
            }
        }  

        if (grid->conway_cells[i].edges.size()>=4) {
            if (grid->conway_cells[i].edges[3].is_boundary == false) {
                if (grid->conway_cells[i].edges[3].neighbour->edges.size()>=3) {
                    if (grid->conway_cells[i].edges[3].neighbour->edges[2].is_boundary == false) {
                        neighbour_count += grid->conway_cells[i].edges[3].neighbour->edges[2].neighbour->getQ();
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
            Qs_new.push_back(grid->conway_cells[i].Q);
        }

    }

    // after calculating now exchange Q_new for Q_old
    for (int i = 0; i < grid->conway_cells.size(); i++) {
        grid->conway_cells[i].Q = Qs_new[i];
    }

}

// requires Q_Cells and 1D cartesian mesh, v_x and v_y >= 0! since upwind scheme doesnt work otherwhise
void Solver::advection_finite_difference_upwind_cartesian(double dt, Point v) {

    if (v.x < 0 || v.y < 0) {
        cerr << "v.x or v.y is < 0. This is not allowed for current version of cartesian upwind advection scheme" << endl;
        exit(EXIT_FAILURE);
    }

    vector<double> new_Qs;

    // calculate all the new Qs
    for (int i = 0; i<grid->q_cells.size(); i++) {

        double Qi_np1 = 0;
        double Qi_n = grid->q_cells[i].Q;
        double dx = grid->q_cells[i].edges[1].length;
        double dy = grid->q_cells[i].edges[2].length;


        if (dt * v.x >= dx || dt * v.y >= dy) {
            cerr << "the CFL timestep condition isnt fullfilled, please make sure that dt < dx/v. otherwise the simulation will show unstable behaviour" << endl;
            exit(EXIT_FAILURE);
        }

        // Qi_np1 = Qi_n - v * (dt/dx) * (Qi_n - 0);

        // double Qim1_n = grid->q_cells[i].edges[0].neighbour->getQ();

        // Qi_np1 = Qi_n - v * (dt/dx) * (Qi_n - Qim1_n);

        if (grid->q_cells[i].edges[0].is_boundary && grid->q_cells[i].edges[3].is_boundary) {
            Qi_np1 = Qi_n - v.x * (dt/dx) * (Qi_n - 0) - v.y * (dt/dy) * (Qi_n - 0);
        } else if (grid->q_cells[i].edges[0].is_boundary) {
            double Qim1y_n = grid->q_cells[i].edges[3].neighbour->getQ();
            Qi_np1 = Qi_n - v.x * (dt/dx) * (Qi_n - 0) - v.y * (dt/dy) * (Qi_n - Qim1y_n);
        } else if (grid->q_cells[i].edges[3].is_boundary) {
            double Qim1x_n = grid->q_cells[i].edges[0].neighbour->getQ();
            Qi_np1 = Qi_n - v.x * (dt/dx) * (Qi_n - Qim1x_n) - v.y * (dt/dy) * (Qi_n - 0);
        } else {
            double Qim1x_n = grid->q_cells[i].edges[0].neighbour->getQ();
            double Qim1y_n = grid->q_cells[i].edges[3].neighbour->getQ();
            Qi_np1 = Qi_n - v.x * (dt/dx) * (Qi_n - Qim1x_n) - v.y * (dt/dy) * (Qi_n - Qim1y_n);

        }

        new_Qs.push_back(Qi_np1);

    }

    // apply all the new Qs
    for (int i = 0; i<grid->q_cells.size(); i++) {

        grid->q_cells[i].Q = new_Qs[i];

    }

}



void Solver::advection_vmesh(double dt, Point v) {

    vector<double> new_Qs;

    // calculate all the new Qs
    for (int i = 0; i < grid->q_cells.size(); i++) {

        double Q_i_np1 = 0;
        double Q_i_n = grid->q_cells[i].Q;
        double total_Flux = 0;

        for (int j = 0; j < grid->q_cells[i].edges.size(); j++) {

            Point n_ij;
            n_ij = get_normal_vec(grid->q_cells[i].edges[j].a, grid->q_cells[i].edges[j].b);
            double c_ij = max(-1 * (v.x*n_ij.x + v.y*n_ij.y), 0.0);
            double l_ij = grid->q_cells[i].edges[j].length;
            double d_ij = 1;
            double Q_j_n;

            if (grid->q_cells[i].edges[j].is_boundary) {

                if (grid->q_cells[i].edges[j].a.x == 0) {
                    d_ij = 2*grid->q_cells[i].seed.x;
                } else if (grid->q_cells[i].edges[j].a.y == 0) {
                    d_ij = 2*grid->q_cells[i].seed.y;
                } else if (grid->q_cells[i].edges[j].a.x == 1) {
                    d_ij = 2*(1-grid->q_cells[i].seed.x);
                } else if (grid->q_cells[i].edges[j].a.y == 1) {
                    d_ij = 2*(1-grid->q_cells[i].seed.y);
                }

                Q_j_n = 0;

            } else {

                double d_ij = sqrt((grid->q_cells[i].seed.x - grid->q_cells[i].edges[j].neighbour->seed.x)*
                               (grid->q_cells[i].seed.x - grid->q_cells[i].edges[j].neighbour->seed.x) + 
                               (grid->q_cells[i].seed.y - grid->q_cells[i].edges[j].neighbour->seed.y)*
                               (grid->q_cells[i].seed.y - grid->q_cells[i].edges[j].neighbour->seed.y));

                Q_j_n = grid->q_cells[i].edges[j].neighbour->getQ();

            }

            total_Flux += c_ij * (1/d_ij) * (Q_i_n - Q_j_n);
            //total_Flux += c_ij * (l_ij/d_ij) * (Q_i_n - Q_j_n);

        }

        Q_i_np1 = Q_i_n - dt*total_Flux;

        new_Qs.push_back(Q_i_np1);

    }


    // apply the calculated Qs
    for (int i = 0; i<grid->q_cells.size(); i++) {

        grid->q_cells[i].Q = new_Qs[i];

    }


}


Solver::~Solver() {};