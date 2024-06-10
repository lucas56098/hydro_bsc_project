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
    bool is_1D = false;
    int N_row = 50*5; // total cells = N^dimension //50*5 

    // simulaiton
    int sim_steps = 5000; //5k is good
    double total_sim_time = 0.2; //0.2 is good
    double dt = static_cast<double>(total_sim_time)/static_cast<double>(sim_steps);
    int save_iter = 50; // 50 maybe

    // GRID GENERATION ------------------------------------
    Mesh<SWE_Cell> grid;
    grid.generate_grid(cartesian, is_1D, N_row, 20, false);

    // INITIAL CONDITIONS ---------------
    for (int i = 100; i < 150; i++) {
        for (int j = 100; j < 150; j++) {
            grid.cells[N_row * i + j].u = 2;
        }
    }
    //grid.cells[50*25+25].h = 1.005;
    //grid.initalize_SWE_gaussian(Point(0.2, 0.2), 1, 0.05);
    //grid.initalize_SWE_gaussian(Point(0.4, 0.6), 1, 0.05);
    //grid.initalize_SWE_gaussian(Point(0.6, 0.1), 1, 0.05);


    //grid.initialize_Q_cells(0, 5, 1, 1);
    //grid.initalize_Q_circle(Point(0.5, 0.5), 0.1);
    //grid.save_Q_diff(0, true, true);
    //grid.save_L1_adv_circle(0, true, Point(0.5/sqrt(2), 0.5/sqrt(2)));
    //grid.save_L1_adv_1Dstepfunc(0, true, 0.5, 0, 0.1);
    
    grid.inialize_boundary_struct(Point(0.7, 0.25), 0.2, 0.5);

    // SIMULATION -----------------------------------------
    Solver<SWE_Cell> solver(&grid);

    for (int i = 0; i<sim_steps; i++) {
        
        // save meshfiles
        if (i%save_iter == 0) {
            grid.save_mesh(i, file_name, dt);
            //grid.save_Q_diff(i * dt, false, true);
            //grid.save_L1_adv_circle(i * dt, false, Point(0.5/sqrt(2), 0.5/sqrt(2)));
            //grid.save_L1_adv_1Dstepfunc(i*dt, false, 0.5, 0, 0.1);
            cout << i << endl;
        }

        // different update steps for the solvers
        //solver.diffusion_like(dt);
        //solver.conway();
        //solver.advection(dt, Point(0.5/sqrt(2), 0.5/sqrt(2)));
        //solver.shallow_water_1D_cartesian(dt);
        //solver.shallow_water_2D(dt);
        solver.shallow_water_2D_cartesian(dt);
    }

    cout << "done" << endl;


    return 0;    
}


