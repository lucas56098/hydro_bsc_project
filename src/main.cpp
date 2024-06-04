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
    bool cartesian = true;
    string file_name = cartesian ? "cmesh" : "vmesh";
    bool is_1D = true;
    int N_row = 10000; // total cells = N^dimension

    // simulaiton
    int sim_steps = 100000;
    double total_sim_time = 0.6;
    double dt = static_cast<double>(total_sim_time)/static_cast<double>(sim_steps);
    int save_iter = 1000;

    // GRID GENERATION ------------------------------------
    Mesh<SWE_Cell> grid;
    grid.generate_grid(cartesian, is_1D, N_row, 50, false);

    // INITIAL CONDITIONS ---------------
    for (int i = 0; i < 1000; i++) {
        grid.cells[4500+i].h = 2;
        //grid.cells[0+i].u = 1;
    }

    //grid.initialize_Q_cells(0, 5, 1, 1);
    //grid.initalize_Q_circle();
    //grid.save_Q_diff(0, true, true);
    //grid.save_L1_adv_circle(0, true, Point(0, 0.5));
    //grid.save_L1_adv_1Dstepfunc(0, true, 0.5, 0, 0.1);



    // SIMULATION -----------------------------------------
    Solver<SWE_Cell> solver(&grid);

    for (int i = 0; i<sim_steps; i++) {
        
        // save meshfiles
        if (i%save_iter == 0) {
            grid.save_mesh(i, file_name, dt);
            //grid.save_Q_diff(i * dt, false, true);
            //grid.save_L1_adv_circle(i * dt, false, Point(0.5, 0));
            //grid.save_L1_adv_1Dstepfunc(i*dt, false, 0.5, 0, 0.1);
            cout << i << endl;
        }

        // different update steps for the solvers
        //solver.diffusion_like(dt);
        //solver.conway();
        //solver.advection(dt, Point(0.5, 0.5));
        solver.shallow_water_1D_cartesian(dt);
    }

    cout << "done" << endl;


    return 0;    
}


