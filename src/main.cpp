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
    string file_name = cartesian ? "finalplot/cmesh_Nbf1_steps200_" : "finalplot/vmesh_";
    bool is_1D = false; 
    bool is_repeating = false;
    bool structure = false;
    int N_row = 250; // total cells = N^dimension
    double sl_beta = 0;

    // simulation
    int sim_steps = N_row*4;//*2*3*30;//*5*5*5;
    double total_sim_time = 0.1;
    double dt = static_cast<double>(total_sim_time)/static_cast<double>(sim_steps);
    int save_iter = 25;

    // GRID GENERATION ------------------------------------
    //Mesh<DG_Q_Cell> grid(1);
    Mesh<DG_Q_Cell> grid(1);
    grid.generate_grid(cartesian, is_1D, N_row, 15, is_repeating, structure);


    // INITIAL CONDITIONS ---------------
    // initial conditions for advection - - - - - - - - - - 
    //grid.initialize_Q_cells(0, 10*2*5);
    //for (int i = 0; i < grid.cells.size(); i++) {
    //    grid.cells[i].Q = exp(-1*(grid.cells[i].seed.x - 0.1)*(grid.cells[i].seed.x - 0.1)/(0.002));
    //}
    //grid.initialize_Q_circle(Point(0.5, 0.5), 0.1);
    //grid.save_Q_diff(0, true, true);
    //grid.save_L1_adv_circle(0, true, Point(0.5/sqrt(2), 0.5/sqrt(2)));
    //grid.save_L1_adv_circle(0, true, Point(0.5, 0));
    //grid.save_L1_adv_1Dstepfunc(0, true, 0.5, 0, 0.1);
    //grid.save_L1_adv_1Dstepfunc(0, true, 0.5, 0, 0.2);
    
    
    // initial conditions for SWE - - - - - - - - - - - - - 
    //grid.initialize_SWE_dam_break(2.0, 1.0, 0.2, 0);
    //grid.initialize_SWE_gaussian(Point(0.25, 0.5), 1, 0.05);
    //grid.initialize_SWE_gaussian(Point(0.25, 0.5), 1, 0.05);
    //grid.initialize_SWE_gaussian(Point(0.75, 0.5), 1, 0.05);
    //grid.initialize_SWE_gaussian(Point(0.1, 0.2), 1, 0.02);
    //grid.save_L1_swe_dam_break(0, true);

    // initial conditions for euler - - - - - - - - - - - -
    //grid.initialize_euler_shock_tube();
    //grid.initialize_rayleigh_taylor(Point(0, -1));
    //grid.initialize_quad_shock();
    //grid.initialize_kelvin_helmholtz();
    //grid.cells[155].rho = 10;
    //grid.cells[130].rho = 2;
    /*double pi = 3.14159265358979323846;
    Point g = Point(0, 0);
    for (int i = 0; i<grid.cells.size(); i++) {
        if (grid.cells[i].seed.y > 0.5 + 0.01*sin(5*2*pi*grid.cells[i].seed.x)) {
            grid.cells[i].rho = 0.1;
            grid.cells[i].u = 0.5;
            grid.cells[i].v = 0.5;
            double P = 1 + grid.cells[i].rho * g.y * (grid.cells[i].seed.y - 0.5);
            grid.cells[i].E = (P/(grid.cells[i].gamma - 1)) + 0.5*grid.cells[i].rho*(grid.cells[i].u*grid.cells[i].u);
        } else {
            grid.cells[i].rho = 5;
            grid.cells[i].u = 0;
            grid.cells[i].v = 0;
            double P = 1 + grid.cells[i].rho * g.y * (grid.cells[i].seed.y - 0.5);
            grid.cells[i].E = (P/(grid.cells[i].gamma - 1)) + 0.5*grid.cells[i].rho*(grid.cells[i].u*grid.cells[i].u);
        }
    }*/
    //grid.initialize_const_flow(Point(0.3, 0));

    // initial conditions for discontinous galerkin advection 1D - - - - - - - -
    //for (int i = 0; i < grid.cells.size()/5; i++) {
    //    double x = grid.cells[i].seed.x;
    //    grid.cells[i].Q(0) = 1;
        //grid.cells[i].Q(0) = exp(-(10 * (x - 0.25)) * (10 * (x - 0.25)));
        //grid.cells[i].Q(1) = exp(-(10 * (x - 0.25)) * (10 * (x - 0.25))) * (x - 0.25) * -200/(2*N_row);
        //grid.cells[i].Q(2) = exp(-(10 * (x - 0.25)) * (10 * (x - 0.25))) * (2000 * (x - 0.25) * (x - 0.25) - 200)/(2*N_row)/(2*N_row);
        
    //}
    

    // initial conditions for discontinous galerkin advection 2D - - - - - - - -
    grid.DG_2D_initialize_gaussian_function(Point(0.5, 0.5), 0.1, 0.1);
    //grid.DG_2D_initialize_step_function();
    

    // initialize boundary condition - - - - - - - - - - - - 
    //grid.initialize_boundary_struct(Point(0.6, 0.3), 0.2, 0.4);
    //grid.initialize_boundary_struct(Point(0.0, 0.98), 1.0, 0.02);
    //grid.initialize_boundary_struct(Point(0.0, 0.0), 1.0, 0.12);
    //grid.initialize_boundary_struct(Point(0.0, 0.0), 0.05, 1);
    //grid.initialize_boundary_struct(Point(0.5, 0.0), 0.5, 1);

    // SIMULATION -----------------------------------------
//    Solver<Q_Cell> solver(&grid);
    DG_Solver<DG_Q_Cell> dgsolver(&grid);
    // initialization step (sets M, S and all the other matricies once)
    dgsolver.advection2D(0, 0, Point(0.5, 0.5), true);
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
            //grid.save_L1_adv_circle(i * dt, false, Point(0.5, 0));
            //grid.save_L1_adv_1Dstepfunc(i*dt, false, 0.5, 0, 0.2);
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
        //solver.advection(dt, Point(0.5, 0));
        //solver.shallow_water(dt, -1, 0, 2);
        //solver.euler(dt, -1, 2, Point(0, 0));
        //dgsolver.advection1D(dt, 0, 1);
        dgsolver.advection2D(dt, sl_beta, Point(0.5, 0.5));
    }
    auto final = chrono::high_resolution_clock::now();
    chrono::duration<double> total_time = final - start;

    cout << "\n Total time: " << total_time.count() << endl;

    cout << "done" << endl;

    return 0;    
}


