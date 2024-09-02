#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <sstream>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "cell_types/DG_Q_Cell.h"
#include "Mesh.h"
#include "Solver.h"
#include "DG_Solver.h"
#include "vmp/VoronoiMesh.h"
#include "utilities/Functions.h"
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;

// MAIN :  -------------------------------------------------------------------------------------------------------
int main () {

    // SPECIFICATIONS -------------------------------------
    // grid
    bool cartesian = true;
    string file_name = cartesian ? "cmesh_Nbf1_res400_" : "vmesh";
    bool is_1D = true; 
    bool is_repeating = false;
    int N_row = 400; // total cells = N^dimension

    // simulation
    int sim_steps = 4000;
    double total_sim_time = 0.35;
    double dt = static_cast<double>(total_sim_time)/static_cast<double>(sim_steps);
    int save_iter = 40;

    // GRID GENERATION ------------------------------------
    Mesh<DG_Q_Cell> grid(1);
    grid.generate_grid(cartesian, is_1D, N_row, 15, is_repeating, false);


    // INITIAL CONDITIONS ---------------
    // initial conditions for advection - - - - - - - - - - 
    //grid.initalize_Q_circle(Point(0.5, 0.5), 0.1);
    //grid.save_Q_diff(0, true, true);
    //grid.save_L1_adv_circle(0, true, Point(0.5/sqrt(2), 0.5/sqrt(2)));
    //grid.save_L1_adv_1Dstepfunc(0, true, 0.5, 0, 0.1);
    
    // initial conditions for SWE - - - - - - - - - - - - - 
    //grid.initialize_SWE_dam_break(2.0, 1.0, 0.5, 0);
    //grid.initalize_SWE_gaussian(Point(0.5, 0.5), 0.5, 0.2);
    //grid.initalize_SWE_gaussian(Point(0.75, 0.25), 1, 0.05);
    //grid.initalize_SWE_gaussian(Point(0.2, 0.4), 1, 0.05);
    //grid.save_L1_swe_dam_break(0, true);

    // initial conditions for euler - - - - - - - - - - - -
    //grid.initialize_euler_shock_tube();
    //grid.initialize_rayleigh_taylor(Point(0, -1));
    //grid.initialize_quad_shock();
    //grid.initialize_kelvin_helmholtz();
    //grid.cells[155].rho = 10;
    //grid.cells[130].rho = 2;
    //for (int i = 0; i<grid.cells.size(); i++) {
    //    grid.cells[i].u = 0.5;
    //    grid.cells[i].v = 0;
    //}
    //grid.initialize_const_flow(Point(0.5, 0));

    // initial conditions for discontinous galerkin advection - - - - - - - -
    for (int i = 0; i < grid.cells.size(); i++) {
        double x = grid.cells[i].seed.x;

        //grid.cells[i].Q(0) = 1;

        grid.cells[i].Q(0) = exp(-(10 * (x - 0.25)) * (10 * (x - 0.25)));
        //grid.cells[i].Q(1) = exp(-(10 * (x - 0.25)) * (10 * (x - 0.25))) * (x - 0.25) * -200/(2*N_row);
        //grid.cells[i].Q(2) = exp(-(10 * (x - 0.25)) * (10 * (x - 0.25))) * (2000 * (x - 0.25) * (x - 0.25) - 200)/(2*N_row)/(2*N_row);
        
    }
    

    // initialize boundary condition - - - - - - - - - - - - 
    //grid.initialize_boundary_struct(Point(0.6, 0.3), 0.02, 0.4);
    //grid.initialize_boundary_struct(Point(0.0, 0.99), 1.0, 0.01);
    //grid.initialize_boundary_struct(Point(0.0, 0.0), 1.0, 0.01);

    // SIMULATION -----------------------------------------
    //Solver<Euler_Cell> solver(&grid);
    DG_Solver<DG_Q_Cell> dgsolver(&grid);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    dgsolver.advection1D(dt, 0, 1);
    // start timer
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i<sim_steps+1; i++) {
        
        // save meshfiles
        if (i%save_iter == 0) {
        
            grid.save_mesh(i, file_name, dt);

            // Q conservation - - - - - - - - - - - 
            //grid.save_Q_diff(i * dt, false, true);
            
            // L1 ERRORS - - - - - - - - - - - - - -
            //grid.save_L1_adv_circle(i * dt, false, Point(0.5/sqrt(2), 0.5/sqrt(2)));
            //grid.save_L1_adv_1Dstepfunc(i*dt, false, 0.5, 0, 0.1);
            //grid.save_L1_swe_dam_break(i*dt, false);

            // print update on simulation progress
            auto now = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = now - start;
            double eta = (elapsed.count() / (i + 1)) * (sim_steps - i - 1);
            cout << "\rStep: " << i << "/" << sim_steps << ", Time: [" << format_time(elapsed.count()) << "<" << format_time(eta) << "]" << flush;
        }

        // different update steps for the solvers - - - - - - - - - - - - - - -
        //solver.diffusion_like(dt);
        //solver.conway();
        //solver.advection(dt, Point(0.5/sqrt(2), 0.5/sqrt(2)));
        //solver.shallow_water(dt, -1, 0, 2);
        //solver.euler(dt, -1, 2, Point(0, -1));
        if (i%1 == 0) {
            dgsolver.advection1D(dt, 1, 1);
        } else {
            dgsolver.advection1D(dt, 0, 1);
        }
    }

    cout << "done" << endl;

    return 0;    
}


