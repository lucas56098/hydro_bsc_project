#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "Mesh.h"
#include "cell_types/Cell.h"
#include "cell_types/Point.h"
#include "cell_types/Q_Cell.h"
#include "cell_types/Conway_Cell.h"
#include "vmp/VoronoiMesh.h"
#include <algorithm>
#include <random>


Mesh::Mesh() {}

Mesh::Mesh(int cell_typein) {

    cell_type = cell_typein;

}
/*  ADAPTION TO DIFFERENT CELL TYPES: -------------------------------------------------------------

    the following functions are written to better handle that the mesh can
    hold different derived Cell types, these functions do the operations for 
    the matching cell type. whenever a new cell type is added these functions
    need to be adapted as well as adding a vector<NewCellType> in the header file
    and adding an option in save_mesh to save the right quantities
    
*/

// push back for different derived cell types
void Mesh::smart_push_back(Point seedin, vector<face> edgesin) {

    if (cell_type == 0) {
        cells.push_back(Cell(seedin, edgesin));
    } else if (cell_type == 1) {
        q_cells.push_back(Q_Cell(seedin, edgesin, 0));
    } else if (cell_type == 2) {
        conway_cells.push_back(Conway_Cell(seedin, edgesin, 0));
    }

}

// get_edge for different derived cell types
face& Mesh::smart_cells_edge(int i, int j) {

    if (cell_type == 1) {
        return q_cells[i].edges[j];
    } else if (cell_type == 2) {
        return conway_cells[i].edges[j];
    } else {
        return cells[i].edges[j];
    }

}

// get cell for different derived cell types
Cell& Mesh::smart_cells(int i) {

    if (cell_type == 1) {
        return q_cells[i];
    } else if (cell_type == 2) {
        return conway_cells[i];
    } else {
        return cells[i];
    }

}

// get_edges.size() for different derived cell types
int Mesh::get_edges_size(int i) {

    if (cell_type == 1) {
        return q_cells[i].edges.size();
    } else if (cell_type == 2) {
        return conway_cells[i].edges.size();
    } else {
        return cells[i].edges.size();
    }

}

// get_cells.size() for different derived cell types
int Mesh::get_cells_size() {
    if (cell_type == 1) {
        return q_cells.size();
    } else if (cell_type == 2) {
        return conway_cells.size();
    } else {
        return cells.size();
    }
}


// GRID GENERATION: -------------------------------------------------------------------------------
// generates a uniform grid with all the neighbour relations and so on
void Mesh::generate_uniform_grid2D(Point start, int n_hor, int n_vert, double distx, double disty) {

    for (int b = 0; b<n_vert; b++) {
        for (int a = 0; a<n_hor; a++) {

            // set the correct seed
            Point seedin(start.x + 0.5*distx + a*distx, start.y + 0.5*disty + b*disty);
            
            // start to define the faces
            vector<face> edgesin;
            face f0; face f1; face f2; face f3;

            // set a and b for the faces
            f0.a = Point(start.x + a*distx, start.y + b*disty);
            f0.b = Point(start.x + a*distx, start.y + b*disty + disty);
            f1.a = Point(start.x + a*distx, start.y + b*disty + disty);
            f1.b = Point(start.x + a*distx + distx, start.y + b*disty + disty);
            f2.a = Point(start.x + a*distx + distx, start.y + b*disty + disty);
            f2.b = Point(start.x + a*distx + distx, start.y + b*disty);
            f3.a = Point(start.x + a*distx + distx, start.y + b*disty);
            f3.b = Point(start.x + a*distx, start.y + b*disty);

            // set correct boundary flags
            f0.is_boundary = false; f1.is_boundary = false; f2.is_boundary = false; f3.is_boundary = false;
            if (a == 0) {f0.is_boundary = true;}
            if (a == n_hor -1) {f2.is_boundary = true;}
            if (b == 0) {f3.is_boundary = true;}
            if (b == n_vert -1) {f1.is_boundary = true;}

            // push faces in edges vector
            edgesin.push_back(f0);
            edgesin.push_back(f1);
            edgesin.push_back(f2);
            edgesin.push_back(f3);

            // set face length to dist

            edgesin[0].length = disty;
            edgesin[1].length = distx;
            edgesin[2].length = disty;
            edgesin[3].length = distx;

            // push back new cell in cells vector
            //cells.push_back(Q_Cell(seedin, edgesin, 0));

            smart_push_back(seedin, edgesin);

        }
    }

    // now that all cells exist we define the neighbour relations
    for (int i = 0; i<get_cells_size(); i++) {

        smart_cells(i).volume = distx * disty;

        for (int j = 0; j<get_edges_size(i); j++) {
            if (smart_cells_edge(i, j).is_boundary == false) {
                if (j == 0) {
                    //cells[i].edges[j].neighbour = &cells[i-1];
                    smart_cells_edge(i, j).neighbour = &smart_cells(i-1);//&q_cells[i-1];
                } else if (j == 1) {
                    //cells[i].edges[j].neighbour = &cells[i+n_hor];
                    smart_cells_edge(i, j).neighbour = &smart_cells(i+n_hor);
                } else if (j == 2) {
                    //cells[i].edges[j].neighbour = &cells[i+1];
                    smart_cells_edge(i, j).neighbour = &smart_cells(i+1);
                } else if (j == 3) {
                    //cells[i].edges[j].neighbour = &cells[i-n_hor];
                    smart_cells_edge(i, j).neighbour = &smart_cells(i-n_hor);
                }
                
            }
        }
    }
}

// generates vmesh using vmp and converts it into data usable for this mesh type
void Mesh::generate_vmesh2D(vector<Point> pts) {

    VoronoiMesh vmesh(pts);
    vmesh.do_point_insertion();

    // loop through all cells to set everything but neighbour relations
    for (int i = 0; i<vmesh.vcells.size(); i++) {

        // set seed and define edge vector
        Point seedin(vmesh.vcells[i].seed.x, vmesh.vcells[i].seed.y);
        vector<face> edgesin;

        for (int j = 0; j<vmesh.vcells[i].edges.size(); j++) {

            face f;

            // set face quantities
            f.a = vmesh.vcells[i].verticies[((vmesh.vcells[i].edges.size()-1) + j)%vmesh.vcells[i].edges.size()];
            f.b = vmesh.vcells[i].verticies[j];
            f.length = sqrt((f.a.x - f.b.x)*(f.a.x - f.b.x) + (f.a.y - f.b.y)*(f.a.y - f.b.y));
            f.is_boundary = (vmesh.vcells[i].edges[j].index1 < 0);

            edgesin.push_back(f);

        }

        smart_push_back(seedin, edgesin);
        //cells.push_back(Q_Cell(seedin, edgesin, 0));

    }

    // loop through everything and set neighbour relations
    for (int i = 0; i<vmesh.vcells.size(); i++) {

        smart_cells(i).volume = vmesh.vcells[i].get_area();

        for (int j = 0; j<vmesh.vcells[i].edges.size(); j++) {

            // neighbour index as used in vmesh
            int neigbour_index;
            neigbour_index = vmesh.vcells[i].edges[j].index2;

            // exclude boundaries
            if (neigbour_index >= 0) {

                // neighbour as defined before over index now with adress
                //cells[i].edges[j].neighbour = &cells[neigbour_index];
                smart_cells_edge(i, j).neighbour = &smart_cells(neigbour_index);

            }

        }

    }


}

// generates a 1D uniform grid with all the neighbour relations and so on
void Mesh::generate_uniform_grid1D(Point start, int n, double dist) {
    generate_uniform_grid2D(start, n, 1, dist, 1);
}

// points x values have to be between 0 and 1 for 1D mesh!!
void Mesh::generate_vmesh1D(vector<Point> pts) {

    for (int i = 0; i < pts.size(); i++) {
        pts[i].y = 0.5;
    }

    vector<Point> sorted_pts = pts;
    sort(sorted_pts.begin(), sorted_pts.end(), [](const Point& a, const Point& b) {
        return a.x < b.x;
    });

    for (int i = 0; i < sorted_pts.size(); i++) {

        // set the correct seed
        Point seedin = sorted_pts[i];

        // start to define the faces
        vector<face> edgesin;
        face f0; face f1; face f2; face f3;

        // set a and b for the faces and boundary flags
        f0.is_boundary = false;
        f1.is_boundary = true;
        f2.is_boundary = false;
        f3.is_boundary = true;

        double distl;
        double distr;

        if (i == 0) {
            f0.is_boundary = true;
            f0.a = Point(0, 0);
            f0.b = Point(0, 1);
            f1.a = Point(0, 1);
            f3.b = Point(0,0);
            distl = seedin.x;
        } else {
            f0.a = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 0);
            f0.b = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 1);
            f1.a = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 1);
            f3.b = Point((sorted_pts[i-1].x + sorted_pts[i].x)/2.0, 0);
            distl = sorted_pts[i].x - sorted_pts[i-1].x;
        }
        
        if (i == sorted_pts.size()-1) {
            f2.is_boundary = true;
            f1.b = Point(1, 1);
            f2.a = Point(1, 1);
            f2.b = Point(1,0);
            f3.a = Point(1,0);
            distr = 1 - seedin.x;
        } else {
            f1.b = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 1);
            f2.a = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 1);
            f2.b = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 0);
            f3.a = Point((sorted_pts[i].x + sorted_pts[i+1].x)/2.0, 0);
            distr = sorted_pts[i+1].x - sorted_pts[i].x;
        }

        edgesin.push_back(f0);
        edgesin.push_back(f1);
        edgesin.push_back(f2);
        edgesin.push_back(f3);

        edgesin[0].length = 1;
        edgesin[2].length = 1;
        edgesin[1].length = distl + distr;
        edgesin[3].length = distl + distr;

        // push back new cell in cells vector
        smart_push_back(seedin, edgesin);

    }

    // now that all cells exist we define the neighbour relations
    for (int i = 0; i < get_cells_size(); i++) {

        smart_cells(i).volume = smart_cells_edge(i, 0).length * smart_cells_edge(i, 1).length;

        if (i != 0) {
            smart_cells_edge(i, 0).neighbour = &smart_cells(i-1);
        }
        if (i != get_cells_size()-1) {
            smart_cells_edge(i, 2).neighbour = &smart_cells(i+1);
        }

    }

}


// SET INITIAL CONDITIONS -------------------------------------------------------------------------
// sets the inital value for the first N cells of the cell vector to value
void Mesh::initialize_cells(int a, int b, double value, int step) {

    if(get_cells_size()< b) {
        cerr << "ERROR: initalize_cells(a, b, value), tried to initalize more cells then there are! b = " << b  << " > cells.size() =" << get_cells_size() << endl;
        exit(EXIT_FAILURE);
    }

    if (cell_type == 1) {
        for (int i = a; i < b; i+=step) {
            q_cells[i].Q = value;
        }
    } else if (cell_type == 2) {
        for (int i = a; i < b; i+=step) {
            conway_cells[i].Q = static_cast<int>(value);
            
        }
        if (static_cast<double>(static_cast<int>(value)) - value != 0) {
            cout <<  "Warning: while setting initial conditions specified double value = " << value << " was casted to " << static_cast<int>(value) << endl;
        }
    }

}

// for Conwway_Cells only! Initalizes grid with random integers 0,1
void Mesh::initialize_random() {

    if (cell_type !=2) {
        cerr << "initialize_random() function does only work with Conway_Cells. Please use them" << endl;
        exit(EXIT_FAILURE);
    }

    // prepare random device
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, 1);

    // fill the Q values with random integers 0, 1
    for (int i = 0; i < get_cells_size(); i++) {
        conway_cells[i].Q = dist(gen) * dist(gen);
    }
}

// SAVE DATA TO FILE ------------------------------------------------------------------------------
// saves mesh into csv file readable for python script
void Mesh::save_mesh(int file_nr, bool cartesian) {

    string filename;

    if (cartesian) {
        filename = "../files/cmesh" + to_string(file_nr) + ".csv"; 
    } else {
        filename = "../files/vmesh" + to_string(file_nr) + ".csv"; 
    }

    ofstream output_file(filename);
    
    // get maximum edge number for column correction later on
    int max_edge_nr = 0;
    for (int i = 0; i<get_cells_size(); i++) {
        if (get_edges_size(i)>=max_edge_nr) {
            max_edge_nr = get_edges_size(i);
        }
    }

    // save the mesh in the following format: ax1, ay1, bx1, by1 ; ax2, ... ; ..; | Q
    output_file << "seed.x, seed.y | a.x, a.y, b.x, b.y ; a.x ... ; | Quantities" << endl;

    for (int i = 0; i<get_cells_size(); i++) {

        output_file << smart_cells(i).seed.x << "," << smart_cells(i).seed.y << "|";

        for (int j = 0; j<get_edges_size(i); j++) {

            output_file << smart_cells_edge(i,j).a.x << ','
                        << smart_cells_edge(i,j).a.y << ','
                        << smart_cells_edge(i,j).b.x << ','
                        << smart_cells_edge(i,j).b.y
                        << ";";
        }

        // correct for empty columns such that Q is always at the same column
        for (int a = 0; a<(max_edge_nr-get_edges_size(i)); a++) {
            output_file << ";";
        }

        // here are the different storing options for the different cell types
        if (cell_type == 1) {

            output_file << "|" << q_cells[i].Q << endl;

        } else if (cell_type == 2){
            output_file << "|" << conway_cells[i].Q << endl;
        } else {
            output_file << endl;
        }

    }

    output_file.close();

}

// function to save the change in the summed up Q value over the whole grid
void Mesh::save_Q_diff(double initial_total_Q, bool reset_file) {

    double total_Q = 0;

    // calculate total_Q summed up over mesh
    for (int i = 0; i < get_cells_size(); i++) {
        total_Q += smart_cells(i).getQ();
    }
    
    // open new file or in append mode
    ofstream output_file;
    if (reset_file) {
        output_file = ofstream("../files/total_Q_diff.csv");
    } else {
        output_file = ofstream("../files/total_Q_diff.csv", ios::app);
    }

    // write caluclated difference in file
    output_file << total_Q-initial_total_Q << ",";
    

}




Mesh::~Mesh() {}