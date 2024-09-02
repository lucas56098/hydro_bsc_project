#include "DG_Q_Cell.h"

DG_Q_Cell::DG_Q_Cell() {}

DG_Q_Cell::DG_Q_Cell(Point seedin, std::vector<face> edgesin, int N_basisfunc) {

    seed = seedin;
    edges = edgesin;
    Q = VectorXd(N_basisfunc);
    
}

// sets M_ij_k for 1D monomial basis functions
MatrixXd DG_Q_Cell::calc_M_ij_k_inv_adv1D() {

    MatrixXd M_ij_k(Q.size(), Q.size());

    double A_k = (edges[1].b.x - edges[1].a.x)/2.0;

    // set values for M_ij_k
    for (int i = 1; i < Q.size()+1; i++) {
        for (int j = 1; j < Q.size()+1; j++) {
            if ((i+j-1)%2 == 0) {
                M_ij_k(i-1, j-1) = 0;
            } else {
                M_ij_k(i-1, j-1) = (2*A_k)/(i+j-1.0);
            }
        }
    }

    return M_ij_k.inverse();

}


// sets S_ij_k for 1D monomial basis functions
MatrixXd DG_Q_Cell::calc_S_ij_k_adv1D() {

    MatrixXd calc_S_ij_k(Q.size(), Q.size());

    // set values for S_ij_k
    for (int i = 1; i < Q.size()+1; i++) {
        for (int j = 1; j < Q.size()+1; j++) {
            if ((i+j-2)%2 == 0) {
                calc_S_ij_k(i-1, j-1) = 0;
            } else {
                calc_S_ij_k(i-1, j-1) = (2*(i-1))/(i+j-2.0);
            }
        }
    }

    return calc_S_ij_k;

}


DG_Q_Cell::~DG_Q_Cell() {}