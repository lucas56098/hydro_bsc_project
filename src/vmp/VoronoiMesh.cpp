#include "VoronoiMesh.h"
#include <fstream>
#include <string>
#include <iostream>
#include <set>

VoronoiMesh::VoronoiMesh(vector<Point> points) {
    pts = points;
    total_steps = 0;
    total_frame_counter = 0;
    //vcells.reserve(pts.size());
}

VoronoiMesh::~VoronoiMesh() {}

// construct all cells using Halfplane Intersection Algorithm
void VoronoiMesh::construct_mesh() {

    vector<int> indices;

    for (int i = 0; i < pts.size(); i++) {
        indices.push_back(i);
    }

    for (int i = 0; i < pts.size(); i++) {

        // construct individual cell and add to vcells vector        
        VoronoiCell vcell(pts[i], i);
        vcell.construct_cell(pts, indices);
        vcells.push_back(vcell);
    }

}

// find cell in which the point is in
int VoronoiMesh::find_cell_index(Point point) {
    
    // guess as the cell im in the last generated cell
    VoronoiCell current_cell = vcells.back();
    double new_cell_index = current_cell.index;
    bool found_cell = true;

    // search for cell
    do {

        // current distance
        double current_cell_dist = sqrt((point.x - current_cell.seed.x)*(point.x - current_cell.seed.x) + (point.y - current_cell.seed.y)*(point.y - current_cell.seed.y));
        double new_cell_dist = current_cell_dist;
        found_cell = true;

        // try to reduce distance step
        for (int i = 0; i<current_cell.edges.size(); i++) {

            int index = current_cell.edges[i].index2;

            if (index >= 0) {
                double dist = sqrt((point.x - pts[index].x)*(point.x - pts[index].x) + 
                                    (point.y - pts[index].y)*(point.y - pts[index].y));
                if (dist < new_cell_dist) {
                    new_cell_dist = dist;
                    new_cell_index = index;
                    found_cell = false;
                }
            }
        }

        current_cell = vcells[new_cell_index];

        total_steps += 1;
    } while (!found_cell);


    return new_cell_index;
}

// function to determine the smallest positive intersection
void VoronoiMesh::find_smallest_pos_intersect(Halfplane &current_hp, int &current_cell_index, VoronoiCell &new_cell, Point &last_vertex, int &last_cell_index, Point &vertex, Halfplane &edge_hp) {
    
    // calculate distance from last_vertex to the midpoint of the current_cell
    double last_vertex_dist_to_midpoint = (last_vertex.x - current_hp.midpoint.x)*current_hp.hp_vec.x
                                        + (last_vertex.y - current_hp.midpoint.y)*current_hp.hp_vec.y;

    // get intersections of current halfplane with all edges of the current cell
    vector<intersection> intersections;
    for (int i = 0; i<vcells[current_cell_index].edges.size(); i++) {
        new_cell.intersect_two_halfplanes(current_hp, vcells[current_cell_index].edges[i], intersections);
    }

    double dist = 42;
    double epsilon = 0.000000000000001;
    bool need_to_check_for_degeneracy = false;

    // find the intersection with smallest but positive relative distance
    for (int i = 0; i<intersections.size(); i++) {

        // calculate relative distance between intersection and last vertex
        double rel_dist = intersections[i].dist_to_midpoint - last_vertex_dist_to_midpoint;

        // reduce dist if rel_dist is smaller
        if (rel_dist > 0 && rel_dist < dist && !(last_cell_index == (*intersections[i].intersecting_with).index2)) {

            if (rel_dist < dist + epsilon && rel_dist > dist - epsilon) {
                need_to_check_for_degeneracy = true;
            }

            dist = rel_dist;
            vertex = intersections[i].intersect_pt;
            edge_hp = *intersections[i].intersecting_with;

        }

    }

    // if degenerate case possible do further checks
    if (need_to_check_for_degeneracy) {
        vector<Halfplane> deg_hp_list;

        // recalculate the distances and store the ones same to dist (the minimal distance with tolearnce)
        for (int i = 0; i < intersections.size(); i++) {
            double rel_dist = intersections[i].dist_to_midpoint - last_vertex_dist_to_midpoint;

            if (rel_dist < dist + epsilon && rel_dist > dist - epsilon) {
                deg_hp_list.push_back(*intersections[i].intersecting_with);
            }
        }

        // out of the degenerate halfplanes find the one with maximum signed angle
        double max_angle = 0;
        for (int i = 0; i < deg_hp_list.size(); i++) {

            // calculate signed angle
            double angle = new_cell.get_signed_angle(current_hp.hp_vec, deg_hp_list[i].hp_vec);

            if (angle > max_angle) {
                max_angle = angle;
                edge_hp = deg_hp_list[i];
            }

        }

    }



    
}

// function to determine the index of an edge in the edge list of its voronoi cell
int VoronoiMesh::get_edge_index_in_cell(int &edge_index, VoronoiCell &vcell) {
    
    int index = -42;

    // loop through edges to find index
    for (int i = 0; i < vcell.edges.size(); i++) {
        if (edge_index == vcell.edges[i].index2) {
            index = i;
        }
    }
    if (index == -42) {
        cout << "failed to find edge in vcell.edge list" << endl;
    }

    return index;
}

// function to do the test if this is the cell where we can leave the boundary again
bool VoronoiMesh::intersection_between_start_stop(int index, Point &new_seed, int &new_seed_index, int &next_cell_index, VoronoiCell &new_cell) {
    // condition here will be: if the test_hp has its boundary intersection between the intersections
    // of the start_hp with boundary and the end_hp with boundary return true. otherwise return false


    // here the start, stop intersection and the hp with boundary intersection will be stored
    vector<intersection> intersections;
    
    // for clarity name all the halfplanes
    Halfplane start_hp = vcells[next_cell_index].edges[(index)%vcells[next_cell_index].edges.size()];
    Halfplane boundary_hp = vcells[next_cell_index].edges[(index+1)%vcells[next_cell_index].edges.size()];
    Halfplane end_hp = vcells[next_cell_index].edges[(index+2)%vcells[next_cell_index].edges.size()];
    Halfplane test_hp = Halfplane(new_seed, vcells[next_cell_index].seed, new_seed_index, next_cell_index);

    // get the intersections
    new_cell.intersect_two_halfplanes(boundary_hp, start_hp, intersections);
    new_cell.intersect_two_halfplanes(boundary_hp, end_hp, intersections);
    new_cell.intersect_two_halfplanes(boundary_hp, test_hp, intersections);

    // again for clarity name the intersections
    double start = intersections[0].dist_to_midpoint;
    double end = intersections[1].dist_to_midpoint;
    double test = intersections[2].dist_to_midpoint;

    // return true if test is in between (<= should handle degeneracys here i think)
    if (test > start && test <= end) {
        return true;
    }
    
    return false;
}

// function to insert a cell into an existing mesh, first generates new_cell and then clips all the cells around as needed
void VoronoiMesh::insert_cell(Point new_seed, int new_seed_index) {

    // find cell the new seed is in
    int cell_im_in_index = find_cell_index(new_seed);
    int current_cell_index = cell_im_in_index;
    
    // generate new_cell and initial halfplane
    VoronoiCell new_cell(new_seed, new_seed_index);
    Halfplane current_hp(new_seed, vcells[cell_im_in_index].seed, new_seed_index, cell_im_in_index);
    Halfplane first_hp = current_hp;

    // set initial last_vertex
    Point last_vertex = current_hp.midpoint;
    int last_cell_index = -42;

    // start going around
    int counter = 0;
    do {

    Point vertex;
    Halfplane edge_hp;

    find_smallest_pos_intersect(current_hp,current_cell_index, new_cell, last_vertex,last_cell_index, vertex, edge_hp);

    // store the found edge and vertex in cell
    new_cell.edges.push_back(current_hp);
    new_cell.verticies.push_back(vertex);

    // if the algorithm hits an edge
    if (edge_hp.boundary) {

        // push back new edge and get its index
        new_cell.edges.push_back(edge_hp);
        int index = get_edge_index_in_cell(edge_hp.index2, vcells[current_cell_index]);

        // go around edge until its time to leave again
        int counter2 = 0;
        bool break_condition = false;
        do {

            // if edge[index+1] is a boundary
            if (vcells[current_cell_index].edges[(index+1)%vcells[current_cell_index].edges.size()].boundary) {

                // push back edge[index + 1]
                new_cell.edges.push_back(vcells[current_cell_index].edges[(index+1)%vcells[current_cell_index].edges.size()]);

                // intersect the two boundary edges to get vertex
                vector<intersection> intersections;
                new_cell.intersect_two_halfplanes(vcells[current_cell_index].edges[index%vcells[current_cell_index].edges.size()], 
                                                    vcells[current_cell_index].edges[(index+1)%vcells[current_cell_index].edges.size()],
                                                    intersections);
                Point new_vertex = intersections[0].intersect_pt;

                // push back vertex
                new_cell.verticies.push_back(new_vertex);


            // if edge[index + 1] is not a boundary
            } else {
                
                // go into next cell
                int next_cell_index = vcells[current_cell_index].edges[(index+1)%vcells[current_cell_index].edges.size()].index2;
                index = get_edge_index_in_cell(current_cell_index, vcells[next_cell_index]);
                current_cell_index = next_cell_index;

            }

            // find out wether its time to leave
            break_condition = intersection_between_start_stop(index, new_seed, new_seed_index, current_cell_index, new_cell);

            if (break_condition) {
                // time to leave the boundary and continue in normal fashion

                // get vertex to restart
                vector<intersection> intersections;
                Halfplane new_hp(new_seed, vcells[current_cell_index].seed, new_seed_index, current_cell_index);
                new_cell.intersect_two_halfplanes(vcells[current_cell_index].edges[(index+1)%vcells[current_cell_index].edges.size()], new_hp, intersections);
                Point restart_vertex = intersections[0].intersect_pt;
                last_vertex = restart_vertex;
                new_cell.verticies.push_back(restart_vertex);

                // update other variables to get ready to leave
                last_cell_index = vcells[current_cell_index].edges[(index+1)%vcells[current_cell_index].edges.size()].index2;
                current_hp = new_hp;

                counter +=1;

            } else {

                // if breack condition is not fulfilled just increase index by one
                index += 1;
                counter2 +=1;

            }

        } while (!(break_condition) && counter2 < 1000);

    } else {

        // some updates of variables before redoing all steps
        last_vertex = vertex;
        last_cell_index = current_cell_index;
        current_cell_index = edge_hp.index2;
        current_hp = Halfplane(new_seed, vcells[current_cell_index].seed, new_seed_index, current_cell_index);
    
        // counting to break do_while when in infinite loop
        counter += 1;
        if (counter >= 10000) {

                    // if counter exeeds limit just generate the cell with the half plane intersection algoithm -> way slower in that case but more robust.
                    cout << "Failed to generate cell using pt_insertion. At: " << pts.size() <<  ". try construct_cell" << endl;
                    VoronoiCell alternative_cell(new_seed, new_seed_index);
                    vector<int> pts_indices;
                    for (int i = 0; i<pts.size(); i++) {
                        pts_indices.push_back(i);
                    }
                    alternative_cell.construct_cell(pts, pts_indices);
                    new_cell = alternative_cell;
                    cout << "cell generated using slow construct_cell algorithm. continue on other pts with point insertion" << endl;
                }
    }

    } while (!(first_hp.index2 == current_hp.index2) && counter < 10000);

    // vector memory management
    new_cell.halfplanes.shrink_to_fit();
    new_cell.edges.shrink_to_fit();
    new_cell.verticies.shrink_to_fit();

    // store new_cell in vcells and new point in pts
    vcells.push_back(new_cell);
    pts.push_back(new_seed);


    // clipping all the neighbour cells
    // loop through all edges of new cell
    for (int i = 0; i<new_cell.edges.size(); i++) {
 
        Halfplane edge = new_cell.edges[i];

        // start and end named in perspective of cell to adapt
        Point v_end = new_cell.verticies[(i-1 + new_cell.verticies.size())%new_cell.verticies.size()];
        Point v_start = new_cell.verticies[i];

        // only adapt cell if not boundary
        if (edge.index2 >= 0) {

            // get cell to adapt, start index and end index and edge to insert
            VoronoiCell cell_to_adapt = vcells[edge.index2];
            int edge_start_index = get_edge_index_in_cell(new_cell.edges[(i+1)%new_cell.edges.size()].index2, cell_to_adapt);
            int edge_end_index = get_edge_index_in_cell(new_cell.edges[(i-1 + new_cell.edges.size())%new_cell.edges.size()].index2, cell_to_adapt);
            Halfplane edge_to_insert(pts[edge.index2], pts[edge.index1], edge.index2, edge.index1);

            // simple case if the end of the vector is not crossed
            if (edge_start_index < edge_end_index) {

                // adapt edges
                cell_to_adapt.edges.erase(cell_to_adapt.edges.begin() + edge_start_index + 1, cell_to_adapt.edges.begin() + edge_end_index);
                cell_to_adapt.edges.insert(cell_to_adapt.edges.begin() + edge_start_index + 1, edge_to_insert);

                // adapt verticies
                cell_to_adapt.verticies.erase(cell_to_adapt.verticies.begin() + edge_start_index, cell_to_adapt.verticies.begin() + edge_end_index);
                cell_to_adapt.verticies.insert(cell_to_adapt.verticies.begin() + edge_start_index, v_start);
                cell_to_adapt.verticies.insert(cell_to_adapt.verticies.begin() + edge_start_index + 1, v_end);

            // more complex case if the end of the vector is crossed during clipping. conceptually doing exactly the same but kind of an index mess
            } else if (edge_start_index > edge_end_index) {

                // erase edges
                if (edge_start_index < cell_to_adapt.edges.size()-1) {
                    cell_to_adapt.edges.erase(cell_to_adapt.edges.begin() + edge_start_index + 1, cell_to_adapt.edges.end());
                }
                cell_to_adapt.edges.erase(cell_to_adapt.edges.begin(), cell_to_adapt.edges.begin() + edge_end_index);

                // erase verticies
                cell_to_adapt.verticies.erase(cell_to_adapt.verticies.begin() + edge_start_index, cell_to_adapt.verticies.end());

                // adapt verticies and edges
                if (edge_end_index > 0) {
                    cell_to_adapt.verticies.erase(cell_to_adapt.verticies.begin(), cell_to_adapt.verticies.begin() + edge_end_index );
                    cell_to_adapt.verticies.insert(cell_to_adapt.verticies.end(), v_start);
                    cell_to_adapt.verticies.insert(cell_to_adapt.verticies.begin(), v_end);
                    cell_to_adapt.edges.insert(cell_to_adapt.edges.begin(), edge_to_insert);
                // also adapt verticies and edges but nothing to delete now
                } else {
                    cell_to_adapt.verticies.insert(cell_to_adapt.verticies.end(), v_start);
                    cell_to_adapt.verticies.insert(cell_to_adapt.verticies.end(), v_end);
                    cell_to_adapt.edges.push_back(edge_to_insert);
                }       
            } else {
                cout << "while adapting cell: start and end index are the same. that shouldnt happen" << endl;
            }

            // restore property that index 0 is the one with seed closest to cell_to_adapt.seed
            // for that calculate distances
            double new_dist = sqrt((edge_to_insert.midpoint.x-cell_to_adapt.seed.x)*(edge_to_insert.midpoint.x-cell_to_adapt.seed.x)
                                    + (edge_to_insert.midpoint.y-cell_to_adapt.seed.y)*(edge_to_insert.midpoint.y-cell_to_adapt.seed.y));
            
            double old_dist = sqrt((cell_to_adapt.edges[0].midpoint.x-cell_to_adapt.seed.x)*(cell_to_adapt.edges[0].midpoint.x-cell_to_adapt.seed.x)
                                    + (cell_to_adapt.edges[0].midpoint.y-cell_to_adapt.seed.y)*(cell_to_adapt.edges[0].midpoint.y-cell_to_adapt.seed.y));

            // resort edge vector if new dist < old dist
            if (new_dist < old_dist) {

                vector<Halfplane> new_edge_list;
                vector<Point> new_verticies_list;

                // get index of inserted edge in cell to adapt
                int index;
                for (int j = 0; j < cell_to_adapt.edges.size(); j++) {
                    if (edge_to_insert.index2 == cell_to_adapt.edges[j].index2) {
                        index = j;
                    }
                }

                // resort the vector
                for (int j = index; j<cell_to_adapt.edges.size()+index; j++) {
                    new_edge_list.push_back(cell_to_adapt.edges[j%cell_to_adapt.edges.size()]);
                    new_verticies_list.push_back(cell_to_adapt.verticies[j%cell_to_adapt.verticies.size()]);
                }
                cell_to_adapt.edges = new_edge_list;
                cell_to_adapt.verticies = new_verticies_list;

            }

            // replace vcell with clipped vcell
            vcells[edge.index2] = cell_to_adapt;
        
        }

        // uncomment this and remove any file saving in main to get animation with clipping
        //save_mesh_to_files(total_frame_counter);
        //total_frame_counter += 1;

    }
   
}


// perform point insertion algorithm on pts
void VoronoiMesh::do_point_insertion() {

    vector<Point> all_pts = pts;
    pts.clear();

    // do the first few points with old algorithm
    for (int i=0; i<3; i++) {
        pts.push_back(all_pts[i]);
    }
    construct_mesh();

    // do point insertion
    for (int i = 3; i<all_pts.size(); i++) {
        if (i%100000 == 0) {
            cout << "progress: " << i << "/" << all_pts.size() << "  -> " << static_cast<double>(i)/static_cast<double>(all_pts.size())*100 << " % " << endl;
        }
        insert_cell(all_pts[i], i);
    }
}


// save the mesh to files (seedfile, edgefile, vertexfile)
void VoronoiMesh::save_mesh_to_files(int nr) {

    // save seeds to file
    string name = "files/seed_list" + to_string(nr) + ".csv";

    ofstream seed_list(name);

    seed_list << "seed_x,seed_y";

    for (int i = 0; i < vcells.size(); i++) {
        seed_list << "\n" << vcells[i].seed.x << "," << vcells[i].seed.y;
    }
    seed_list.close();

    // save vertices to file
    name = "files/vertex_list" + to_string(nr) + ".csv";

    ofstream vertex_list(name);

    vertex_list << "vertex_x,vertex_y";

    // loop through cells
    for (int i = 0; i<vcells.size(); i++) {
        
        // loop through verticies for vertex list
        for (int j = 0; j<vcells[i].verticies.size(); j++) {
            double x = vcells[i].verticies[j].x;
            double y = vcells[i].verticies[j].y;
            vertex_list << "\n" << x << "," << y;
        }
    }
    
    vertex_list.close();

    // save edges to file
    name = "files/edge_list" + to_string(nr) + ".csv";

    ofstream edge_list(name);

    edge_list << "edge1_x, edge1_y, edge2_x, edge2_y";

    // loop through cells
    for (int i = 0; i<vcells.size(); i++) {
        
        int nr = vcells[i].verticies.size();

        // loop through edges
        for (int j = 0; j<nr; j++) {
            double x1 = vcells[i].verticies[j].x;
            double y1 = vcells[i].verticies[j].y;
            double x2 = vcells[i].verticies[(j+1)%nr].x;
            double y2 = vcells[i].verticies[(j+1)%nr].y;

            edge_list << "\n" << x1 << "," << y1 << "," << x2 << "," << y2;
        }
    }

    edge_list.close();

}

// check equidistance
bool VoronoiMesh::check_equidistance() {
    bool correct_mesh = true;

    // check every cell individually
    for (int i = 0; i < vcells.size(); i++) {

        // check cell for its conditions
        bool correct_cell = vcells[i].check_equidistance_condition(pts);

        // if a single cell is false the mesh is also false
        if (!correct_cell) {
            correct_mesh = false;
        }
    }

    return correct_mesh;
}

// add up areas of all cells
double VoronoiMesh::check_area() {

    double total_area = 0;

    for (int i = 0; i < vcells.size(); i++) {

        total_area += vcells[i].get_area();
    }

    return total_area;

}

// check that all neighbours know each other
bool VoronoiMesh::check_neighbours() {

    bool all_neighbours_known = false;

    int nr_of_known_neighbours = 0;
    int nr_of_checked_edges = 0;
    

    // go through each cell
    for (int i = 0; i < vcells.size(); i++) {

        // through each edge (first edge)
        for (int j = 0; j < vcells[i].edges.size(); j++) {

            // exclude boundaries
            if (vcells[i].edges[j].boundary == false) {

                nr_of_checked_edges += 1;

                // check edges (seconde edge) of neigbour cell
                for (int k = 0; k < vcells[vcells[i].edges[j].index2].edges.size(); k++) {

                    // if there is a second edge corresponding to the first edge add nr of known neigbours +1
                    if (vcells[vcells[i].edges[j].index2].edges[k].boundary == false && 
                        vcells[vcells[i].edges[j].index2].edges[k].index2 == vcells[i].edges[j].index1) {

                            nr_of_known_neighbours +=1;

                        }

                }

            }
            
        }

    }

    // check wether or not every neighbour knows its neighbour
    if (nr_of_checked_edges == nr_of_known_neighbours) {
        all_neighbours_known = true;
    }

    return all_neighbours_known;
}

// check some conditions for mesh
bool VoronoiMesh::check_mesh() {

    cout << "checking mesh..." << endl;

    bool correct_mesh = true;

    // first check : equidistance
    bool equidist = true;
    if (!check_equidistance()) {
        correct_mesh = false;
        equidist = false;
    }
    cout << "equidistance condition: " << boolalpha << equidist << endl;
    
    // second check : area
    double total_area = check_area();
    cout << "total area = " << total_area << "+" << total_area - 1 << endl;
    if (total_area > 1 + 0.000001 || total_area < 1 - 0.000001) {
        correct_mesh = false;
    }

    // third check : neighbours
    bool neighbour = true;
    if (!check_neighbours()) {
        correct_mesh = false;
        neighbour = false;
    }
    cout << "neighbour condition: " << boolalpha << neighbour << endl;

    // return total check outcome
    return correct_mesh;
}

// function to optimize the memory (kinda only makes sense after generating the mesh and not while)
void VoronoiMesh::optimize_mesh_memory() {

    pts.shrink_to_fit(); 
    vcells.shrink_to_fit();

    for (int i = 0; i<vcells.size(); i++) {

        vcells[i].edges.shrink_to_fit();
        vcells[i].halfplanes.shrink_to_fit();
        vcells[i].verticies.shrink_to_fit();

    }

}

long long VoronoiMesh::calculate_mesh_memory(bool use_capacity) {

    long long total_size;
    if (use_capacity) {
    
            total_size = sizeof(VoronoiMesh) + sizeof(Point)*pts.capacity();
    
    } else {
        
        total_size = sizeof(VoronoiMesh) + sizeof(Point)*pts.size();
    
    }

    for (int i = 0; i<vcells.size(); i++) {

        total_size += vcells[i].calculate_cell_memory(use_capacity);

    }

    return total_size;

}