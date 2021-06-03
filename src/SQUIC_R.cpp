// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#define ARMA_64BIT_WORD 1
#include "RcppArmadillo.h"
#include <stdio.h>
#include <assert.h>
#include <algorithm>


// SQUIC Library iterface
extern "C"
{
    void SQUIC_CPP(
        int mode,
        long p,
        long n, double *Y,
        double lambda,
        long *M_rinx, long *M_cptr, double *M_val, long M_nnz,
        int max_iter, double inv_tol, double term_tol, int verbose,
        long *&X_rinx, long *&X_cptr, double *&X_val, long &X_nnz,
        long *&W_rinx, long *&W_cptr, double *&W_val, long &W_nnz,
        int &info_num_iter,
        double *&info_times,      //length must be 6: [time_total,time_impcov,time_optimz,time_factor,time_aprinv,time_updte]
        double *&info_objective, // length must be size max_iter
        double &info_logdetX,
        double &info_trSX);
}

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
List SQUIC_R(arma::mat &Y, double lambda, int max_iter, double inv_tol, double term_tol, int verbose, int mode, arma::sp_mat &M, arma::sp_mat &X0, arma::sp_mat &W0) {

    // Get the key size parameters
    long p =  Y.n_rows;
    long n = Y.n_cols;

    // Make matrix pointers
    long *M_rinx ;
    long *M_cptr;
    double *M_val;
    long M_nnz = 0;

    long *X_rinx ;
    long *X_cptr;
    double *X_val;
    long X_nnz = 0;

    long *W_rinx ;
    long *W_cptr;
    double *W_val;
    long W_nnz = 0;


    //COPY M MATRIX
    if (M.n_nonzero > 0) {

        M.sync();
        auto M_row_indices = arma::access::rwp(M.row_indices);
        auto M_col_ptrs = arma::access::rwp(M.col_ptrs);
        auto M_values = arma::access::rwp(M.values);
        M_nnz = M.n_nonzero;
        {
            M_rinx = new long[M_nnz];
            M_cptr = new long[p + 1];
            M_val = new double[M_nnz];

            for (long i = 0; i < M_nnz; ++i) {
                M_rinx[i] = M_row_indices[i];
                M_val[i] = M_values[i];
            }

            for (long i = 0; i < p + 1; ++i) {
                M_cptr[i] = M_col_ptrs[i];
            }
        }
    }

    // COPY X0 MATRIX
    X0.sync();
    auto X_row_indices = arma::access::rwp(X0.row_indices);
    auto X_col_ptrs = arma::access::rwp(X0.col_ptrs);
    auto X_values = arma::access::rwp(X0.values);
    X_nnz = X0.n_nonzero;
    {
        X_rinx = new long[X_nnz];
        X_cptr = new long[p + 1];
        X_val = new double[X_nnz];

        for (long i = 0; i < X_nnz; ++i) {
            X_rinx[i] = X_row_indices[i];
            X_val[i] = X_values[i];
        }

        for (long i = 0; i < p + 1; ++i) {
            X_cptr[i] = X_col_ptrs[i];
        }
    }

    // COPY W0 MATRIX
    W0.sync();
    auto W_row_indices = arma::access::rwp(W0.row_indices);
    auto W_col_ptrs = arma::access::rwp(W0.col_ptrs);
    auto W_values = arma::access::rwp(W0.values);
    W_nnz = W0.n_nonzero;
    {
        W_rinx = new long[W_nnz];
        W_cptr = new long[p + 1];
        W_val = new double[W_nnz];

        for (long i = 0; i < W_nnz; ++i) {
            W_rinx[i] = W_row_indices[i];
            W_val[i] = W_values[i];
        }

        for (long i = 0; i < p + 1; ++i) {
            W_cptr[i] = W_col_ptrs[i];
        }
    }

    // Default Result Valuess
    int    info_num_iter = 0;
    double info_logdetX  = 0.0;
    double info_trSX     = 0.0;
    double *info_times_buffer = new double[6]();
    double *info_objective_buffer = new double[std::max(1, max_iter)]();

    // Run SQUIC
    SQUIC_CPP(
        mode,
        p,
        n, Y.memptr(),
        lambda,
        M_rinx, M_cptr, M_val, M_nnz,
        max_iter, inv_tol, term_tol, verbose,
        X_rinx, X_cptr, X_val, X_nnz,
        W_rinx, W_cptr, W_val, W_nnz,
        info_num_iter,
        info_times_buffer,
        info_objective_buffer,
        info_logdetX,
        info_trSX);


    // Copy data it standard format
    // In order to access the internal arrays of the SpMat class call .sync()
    arma::SpMat<double> iC(p, p);
    arma::SpMat<double> C(p, p);
    iC.sync();
    C.sync();

    // Making space for the elements
    iC.mem_resize(X_nnz);
    C.mem_resize(W_nnz);

    // Copying elements
    std::copy(X_rinx, X_rinx + X_nnz, arma::access::rwp(iC.row_indices));
    std::copy(X_cptr, X_cptr + p + 1, arma::access::rwp(iC.col_ptrs));
    std::copy(X_val, X_val + X_nnz, arma::access::rwp(iC.values));
    arma::access::rw(iC.n_rows) = p;
    arma::access::rw(iC.n_cols) = p;
    arma::access::rw(iC.n_nonzero) = X_nnz;

    std::copy(W_rinx, W_rinx + W_nnz, arma::access::rwp(C.row_indices));
    std::copy(W_cptr, W_cptr + p + 1, arma::access::rwp(C.col_ptrs));
    std::copy(W_val, W_val + W_nnz, arma::access::rwp(C.values));
    arma::access::rw(C.n_rows) = p;
    arma::access::rw(C.n_cols) = p;
    arma::access::rw(C.n_nonzero) = W_nnz;

    Rcpp::List output;

    if (max_iter == 0) { // Special max_iter==0: SQUIC only compute the sparse sample covariance S
        output = Rcpp::List::create(
                     Named("S") = C,
                     Named("info_time_total") = info_times_buffer[0],
                     Named("info_time_sample_cov") = info_times_buffer[1]);
    } else { // Regular case return all values

        //Copy info_objective_buffer keeping only info_num_iter elements
        arma::Col<double> info_objective(info_objective_buffer, info_num_iter);
        output = Rcpp::List::create(
                     Named("X") = iC,
                     Named("W") = C,
                     Named("info_time_total") = info_times_buffer[0],
                     Named("info_time_sample_cov") = info_times_buffer[1],
                     Named("info_time_optimize") = info_times_buffer[2],
                     Named("info_time_factor") = info_times_buffer[3],
                     Named("info_time_approximate_inv") = info_times_buffer[4],
                     Named("info_time_coordinate_upd") = info_times_buffer[5],
                     Named("info_objective") = info_objective,
                     Named("info_logdetX") = info_logdetX,
                     Named("info_trSX") = info_trSX
                 );
    }

    // Delete Buffers
    if (M.n_nonzero > 0) {
        delete[] M_rinx;
        delete[] M_cptr;
        delete[] M_val;
    }

    delete[] X_rinx;
    delete[] X_cptr;
    delete[] X_val;

    delete[] W_rinx;
    delete[] W_cptr;
    delete[] W_val;

    delete[] info_times_buffer;
    delete[] info_objective_buffer;

    return output;
}