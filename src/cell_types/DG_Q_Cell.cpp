#include "DG_Q_Cell.h"

DG_Q_Cell::DG_Q_Cell() {}

DG_Q_Cell::DG_Q_Cell(Point seedin, std::vector<edge> edgesin, int N_basisfunc) {

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


// sets M_ij_k for 2D legendre polynomials (at the moment just N_p = 0 and N_p = 1)
MatrixXd DG_Q_Cell::calc_M_ij_k_inv_adv2D() {

    if (Q.size() != 1 && Q.size() != 3 && Q.size() != 6) {
        cerr << "calc_M_ij_k_inv_adv2D does not work for " << Q.size() << " basis functions. Set it to 1 or 3 or 6!" << endl;
        EXIT_FAILURE;
    }

    MatrixXd M_ij_k_inverse(Q.size(), Q.size());
    M_ij_k_inverse.setZero();

    // i am absorbing the N_row into the M_ij matrix here
    double N_row = 1.0/(edges[1].b.x - edges[1].a.x);

    // set values for M_ij_k_inv
    M_ij_k_inverse(0, 0) = N_row;
    if (Q.size() >= 3) {
        M_ij_k_inverse(1, 1) = 3 * N_row;
        M_ij_k_inverse(2, 2) = 3 * N_row;
    }
    if (Q.size() >= 6) {
        M_ij_k_inverse(3, 3) = 9 * N_row;
        M_ij_k_inverse(4, 4) = 5 * N_row;
        M_ij_k_inverse(5, 5) = 5 * N_row;
    }

    return M_ij_k_inverse;
}

// sets S_ij_k for 2D legendre polynomials (at the moment just for N_p = 0 and N_p = 1)
MatrixXd DG_Q_Cell::calc_S_ij_k_adv2D(Point u) {

    if (Q.size() != 1 && Q.size() != 3 && Q.size() != 6) {
        cerr << "calc_M_ij_k_inv_adv2D does not work for " << Q.size() << " basis functions. Set it to 1 or 3 or 6!" << endl;
        EXIT_FAILURE;
    }

    MatrixXd S_ij(Q.size(), Q.size());
    S_ij.setZero();

    // set values for S_ij_k
    if (Q.size() >= 3) {
        S_ij(1, 0) = 2*u.x;
        S_ij(2, 0) = 2*u.y;
    }
    if (Q.size() >= 6) {
        S_ij(3, 2) = (2.0/3.0) * u.x;
        S_ij(3, 1) = (2.0/3.0) * u.y;
        S_ij(4, 1) = 2 * u.x;
        S_ij(5, 2) = 2 * u.y;
    }

    return S_ij;
}

DG_Q_Cell::~DG_Q_Cell() {}