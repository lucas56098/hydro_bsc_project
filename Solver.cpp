#include <iostream>
#include <vector>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "cell_types/Q_Cell.h"
#include "Mesh.h"
#include "Solver.h"

Solver::Solver() {}

Solver::Solver(Mesh* gridin) {

    grid = gridin;

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

// requires Q_Cells and 1D cartesian mesh
void Solver::advection_finite_difference_upwind_cartesian1D(double dt, double v) {

    vector<double> new_Qs;

    // calculate all the new Qs
    for (int i = 0; i<grid->q_cells.size(); i++) {

        double Qi_np1 = 0;
        double Qi_n = grid->q_cells[i].Q;
        double dx =grid->q_cells[i].edges[1].length;

        if (dt >= dx/v) {
            cerr << "the CFL timestep condition isnt fullfilled, please make sure that dt < dx/v. otherwise the simulation will show unstable behaviour" << endl;
            exit(EXIT_FAILURE);
        }

        if (grid->q_cells[i].edges[0].is_boundary) {

            Qi_np1 = Qi_n - v * (dt/dx) * (Qi_n - 0);

        } else {

            double Qim1_n = grid->q_cells[i].edges[0].neighbour->getQ();

            // update equation
            Qi_np1 = Qi_n - v * (dt/dx) * (Qi_n - Qim1_n);

        }

        new_Qs.push_back(Qi_np1);

    }

    // apply all the new Qs
    for (int i = 0; i<grid->q_cells.size(); i++) {

        grid->q_cells[i].Q = new_Qs[i];

    }

}

Solver::~Solver() {}