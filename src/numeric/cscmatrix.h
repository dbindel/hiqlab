/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef CSCMATRIX_H
#define CSCMATRIX_H

#include <vector>
#include <algorithm>

#include "qcomplex.h"


class SpVecStructAssembler;
class SpVecAssembler;


/** Compressed sparse column matrix.
 */
class CSCMatrix {
public:

    /** Construct CSC matrix from basic arrays.
     */
    CSCMatrix(int real_flag = 0);
    CSCMatrix(int m, int n, int nnz, int real_flag = 0);

    template<class idx, class real>
    CSCMatrix(idx* jc, idx* ir, real* pr, real* pi, int m, int n);

    CSCMatrix(const CSCMatrix& matrix);

    int     get_m()   { return m;  }
    int     get_n()   { return n;  }
    int*    get_jc()  { return &(jc[0]); }
    int*    get_ir()  { return &(ir[0]); }
    double* get_Ax()  { return &(pr[0]); }
    double* get_Az()  { return real_flag ? NULL : &(pi[0]); }
    int     get_nnz() { return ir.size(); }
    int     is_real() { return real_flag;  }

    void    set_dims(int xm, int xn);
    void    set_nnz(int xnnz);

    /** Rescale the matrix */
    void scale(double* du, double* df);

    /** Form ax = A*x. */
    void apply(double* xx, double* xz, double* axx, double* axz);
    void apply(dcomplex* x, dcomplex* ax);
    void apply(double* x, double* ax);

    void dump(const char* fname);
    void load(const char* fname);
    void load(int* ii, int* jj, double* rdata, double* idata, int nnz);

    CSCMatrix* mul(CSCMatrix* B, CSCMatrix* result = NULL);
    CSCMatrix* mul(CSCMatrix* Y, CSCMatrix* X, CSCMatrix* result);
    CSCMatrix* add(double alpha, double beta, CSCMatrix* B, 
                   CSCMatrix* result = NULL);
    
    CSCMatrix* transpose();
    CSCMatrix* submatrix(int i1, int i2, int j1, int j2);
    friend CSCMatrix* vcat(CSCMatrix* A, CSCMatrix* B);
    friend CSCMatrix* hcat(CSCMatrix* A, CSCMatrix* B);

private:
    int m, n, real_flag;
    std::vector<int> jc;
    std::vector<int> ir;
    std::vector<double> pr;
    std::vector<double> pi;

    template<class idx, class real>
    void copy_data(idx* xjc, idx* xir, real* xpr, real* xpi, int xn);

    int nnz(int j) { return jc[j+1]-jc[j]; }

    int mul_col_struct(SpVecStructAssembler& Axs, CSCMatrix& X, int k);
    int mul_col_struct(SpVecStructAssembler& Axs, 
                       std::vector<int>& irx, int nnzx);

    void mul_col_real(SpVecAssembler& Axr, CSCMatrix& X, int k);
    void mul_col_real(SpVecAssembler& Axr, 
                      SpVecAssembler& xr, SpVecAssembler& xi,
                      bool real_x, std::vector<int>& irx, int nnzx);

    void mul_col_imag(SpVecAssembler& Axr, CSCMatrix& X, int k);
    void mul_col_imag(SpVecAssembler& Axr, 
                      SpVecAssembler& xr, SpVecAssembler& xi,
                      bool real_x, std::vector<int>& irx, int nnzx);

    int nnz_mul(CSCMatrix& X);
    int nnz_mul(CSCMatrix& Y, CSCMatrix& X);

    void mul(CSCMatrix& X, CSCMatrix& AX, bool build_index);
    void mul(CSCMatrix& Y, CSCMatrix& X, CSCMatrix& AX, bool build_index);

    int nnz_add(CSCMatrix& B);
    void add(double alpha, double beta,
             CSCMatrix& B, CSCMatrix& C, bool build_index);
};


template<class idx, class real>
    CSCMatrix::CSCMatrix(idx* xjc, idx* xir, real* xpr, real* xpi, 
                         int xm, int xn)
{
    copy_data(xjc, xir, xpr, xpi, xn);
    m = xm;
}

template<class idx, class real>
    void CSCMatrix::copy_data(idx* xjc, idx* xir, real* xpr, real* xpi, int xn)
{
    using std::copy;
    int nnz = xjc[xn];
    real_flag = (xpi == NULL);
    set_dims(xn, xn);
    set_nnz(xjc[n]);
    copy(xjc, xjc+n+1, jc.begin());
    copy(xir, xir+nnz, ir.begin());
    copy(xpr, xpr+nnz, pr.begin());
    if (xpi)
        copy(xpi, xpi+nnz, pi.begin());
}


#endif /* CSCMATRIX_H */
