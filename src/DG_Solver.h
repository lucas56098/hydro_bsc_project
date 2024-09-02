#include <iostream>
#include <vector>
#include <array>
#include <string>
#include "cell_types/Point.h"
#include "cell_types/Cell.h"
#include "Mesh.h"

#ifndef DG_Solver_h
#define DG_Solver_h

template <typename CellType>
class DG_Solver {

public:
    DG_Solver<CellType>();
    DG_Solver<CellType>(Mesh<CellType>* gridin);
    ~DG_Solver<CellType>(); 

    void advection1D(double dt, double slope_limited = 0, double u = 1);
    double minmod(double a, double b, double c);
    double get_element_average_adv1D(int index, vector<VectorXd>& new_Qs);

private:

    // grid to solve the equations on
    Mesh<CellType>* grid;

};

#include "DG_Solver.tpp"

#endif