#include <iostream>
#include <vector>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "Mesh.h"
#include "Solver.h"
#include "vmp/VoronoiMesh.h"
#include "utilities/Functions.h"
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
    int sim_steps = 600;
    double total_sim_time = 2;
    double dt = static_cast<double>(total_sim_time)/static_cast<double>(sim_steps);
    int save_iter = 10;

    // GRID GENERATION ------------------------------------
    Mesh<Q_Cell> grid;
    grid.generate_grid(cartesian, is_1D, N_row, 50, true);

    // INITIAL CONDITIONS ---------------
    //grid.initialize_Q_cells(0, 5, 1, 1);
    grid.initalize_Q_circle();
    grid.save_Q_diff(0, true, true);
    //grid.save_L1_adv_circle(0, true, Point(0, 0.5));
    //grid.save_L1_adv_1Dstepfunc(0, true, 0.5, 0, 0.1);



    // SIMULATION -----------------------------------------
    Solver<Q_Cell> solver(&grid);

    for (int i = 0; i<sim_steps; i++) {
        
        // save meshfiles
        if (i%save_iter == 0) {
            grid.save_mesh(i, file_name, dt);
            grid.save_Q_diff(i * dt, false, true);
            //grid.save_L1_adv_circle(i * dt, false, Point(0.5, 0));
            //grid.save_L1_adv_1Dstepfunc(i*dt, false, 0.5, 0, 0.1);
            cout << i << endl;
        }

        // different update steps for the solvers
        //solver.diffusion_like(dt);
        //solver.conway();
        solver.advection(dt, Point(0.5, 0.5));
    }

    cout << "done" << endl;


    return 0;    
}


