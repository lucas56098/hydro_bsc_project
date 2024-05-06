#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "Mesh.h"
#include "Cell.h"
#include "Point.h"
#include "vmp/VoronoiMesh.h"


Mesh::Mesh() {}

// generates a uniform grid with all the neighbour relations and so on
void Mesh::generate_uniform_grid(Point start, int n_hor, int n_vert, double dist) {

    for (int b = 0; b<n_vert; b++) {
        for (int a = 0; a<n_hor; a++) {

            // set the correct seed
            Point seedin(start.x + 0.5*dist + a*dist, start.y + 0.5*dist + b*dist);
            
            // start to define the faces
            vector<face> edgesin;
            face f0; face f1; face f2; face f3;

            // set a and b for the faces
            f0.a = Point(start.x + a*dist, start.y + b*dist);
            f0.b = Point(start.x + a*dist, start.y + b*dist + dist);
            f1.a = Point(start.x + a*dist, start.y + b*dist + dist);
            f1.b = Point(start.x + a*dist + dist, start.y + b*dist + dist);
            f2.a = Point(start.x + a*dist + dist, start.y + b*dist + dist);
            f2.b = Point(start.x + a*dist + dist, start.y + b*dist);
            f3.a = Point(start.x + a*dist + dist, start.y + b*dist);
            f3.b = Point(start.x + a*dist, start.y + b*dist);

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
            for (int i = 0; i<4; i++) {
                edgesin[i].length = dist;
            }

            // push back new cell in cells vector
            cells.push_back(Cell(seedin, edgesin));

        }
    }

    // now that all cells exist we define the neighbour relations
    for (int i = 0; i<cells.size(); i++) {
        for (int j = 0; j<cells[i].edges.size(); j++) {
            if (cells[i].edges[j].is_boundary == false) {
                if (j == 0) {
                    cells[i].edges[j].neighbour = &cells[i-1];
                } else if (j == 1) {
                    cells[i].edges[j].neighbour = &cells[i+n_hor];
                } else if (j == 2) {
                    cells[i].edges[j].neighbour = &cells[i+1];
                } else if (j == 3) {
                    cells[i].edges[j].neighbour = &cells[i-n_hor];
                }
                
            }
        }
    }
}

// converts voronoi mesh data into data usable for this mesh
void Mesh::generate_from_vmesh(VoronoiMesh vmesh) {

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

        cells.push_back(Cell(seedin, edgesin));

    }

    // loop through everything and set neighbour relations
    for (int i = 0; i<vmesh.vcells.size(); i++) {
        for (int j = 0; j<vmesh.vcells[i].edges.size(); j++) {

            // neighbour index as used in vmesh
            int neigbour_index;
            neigbour_index = vmesh.vcells[i].edges[j].index2;

            // exclude boundaries
            if (neigbour_index >= 0) {

                // neighbour as defined before over index now with adress
                cells[i].edges[j].neighbour = &cells[neigbour_index];

            }

        }

    }


}

// function to update Q in a general way
void Mesh::updateQ() {
    
    // calculate all the fluxes
    for (int i = 0; i<cells.size(); i++) {
        cells[i].calc_deltaQs();
    }

    // apply all the fluxes
    for (int i = 0; i<cells.size(); i++) {
        cells[i].apply_deltaQs();
    }

}


void Mesh::save_mesh(int file_nr) {

    string filename = "../files/mesh" + to_string(file_nr) + ".csv"; 

    ofstream output_file(filename);
    
    // get maximum edge number for column correction later on
    int max_edge_nr = 0;
    for (int i = 0; i<cells.size(); i++) {
        if (cells[i].edges.size()>=max_edge_nr) {
            max_edge_nr = cells[i].edges.size();
        }
    }

    // save the mesh in the following format: ax1, ay1, bx1, by1 ; ax2, ... ; ..; | Q
    for (int i = 0; i<cells.size(); i++) {

        for (int j = 0; j<cells[i].edges.size(); j++) {

            output_file << cells[i].edges[j].a.x << ','
                        << cells[i].edges[j].a.y << ','
                        << cells[i].edges[j].b.x << ','
                        << cells[i].edges[j].b.y
                        << ";";
        }

        // correct for empty columns such that Q is always at the same column
        for (int a = 0; a<(max_edge_nr-cells[i].edges.size()); a++) {
            output_file << ";";
        }

        output_file << cells[i].Q << endl;

    }

    output_file.close();

}

Mesh::~Mesh() {}