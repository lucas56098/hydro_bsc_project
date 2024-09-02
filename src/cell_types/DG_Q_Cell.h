#ifndef DG_Q_Cell_h
#define DG_Q_Cell_h
#include "Cell.h"
#include "../Eigen/Dense"
#include <iostream>

using namespace Eigen;


class DG_Q_Cell : public Cell {
public:
    DG_Q_Cell();
    DG_Q_Cell(Point seedin, std::vector<face> edgesin, int N_basisfunc);
    ~DG_Q_Cell();

    // coefficient vectors of q(x,t) in phi basis representation
    VectorXd Q;
    MatrixXd M_ij_k_inv;
    MatrixXd S_ij_k;

    VectorXd get_coeff() const {
        return Q;
    }


    MatrixXd calc_M_ij_k_inv_adv1D();
    MatrixXd calc_S_ij_k_adv1D();

};

#endif
