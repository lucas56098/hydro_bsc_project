#include <iostream>
#include <vector>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "Mesh.h"
#include "Solver.h"
#include "vmp/VoronoiMesh.h"
using namespace std;

// MAIN :  -------------------------------------------------------------------------------------------------------
int main () {
    

    // SPECIFICATIONS -------------------------------------
    // grid
    bool cartesian = false;
    string file_name = cartesian ? "cmesh" : "vmesh";
    bool is_1D = false;
    int N_row = 50; // total cells = N^dimension

    // simulaiton
    int sim_steps = 2000;
    double total_sim_time = 0.7;
    double dt = static_cast<double>(total_sim_time)/static_cast<double>(sim_steps);
    int save_iter = 100;

    // GRID GENERATION ------------------------------------
    Mesh<Q_Cell> grid;
    grid.generate_grid(cartesian, is_1D, N_row);

    // INITIAL CONDITIONS ---------------
    grid.initialize_Q_cells(0, 50*30, 1, 1);
    grid.save_Q_diff(true, true);

    // SIMULATION -----------------------------------------
    Solver<Q_Cell> solver(&grid);

    for (int i = 0; i<sim_steps; i++) {
        
        // save meshfiles
        if (i%save_iter == 0) {
            grid.save_mesh(i, file_name);
            grid.save_Q_diff(false, true);
            cout << i << endl;
        }

        // different update steps for the solvers
        //solver.diffusion_like(dt);
        //solver.conway();
        solver.advection(dt, Point(0.5, 0.3));
    }

    cout << "done" << endl;
    return 0;    
}


