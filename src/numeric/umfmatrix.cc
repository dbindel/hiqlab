/* HiQLab
 * Copyright (c): Regents of the University of California
 */

extern "C" {
  #include "umfpack.h"
}

#include <cstring>
#include <cstdio>
#include <cassert>
#include <vector>

#include "cscmatrix.h"
#include "umfmatrix.h"
#include "qcomplex.h"

using std::vector;


#define ME UMFMatrix


ME::ME(CSCMatrix* cscmat) :
    n(cscmat->get_n()),
    jc(cscmat->get_jc()),
    ir(cscmat->get_ir()),
    Ax(cscmat->get_Ax()),
    Az(cscmat->get_Az())
{
    umfpack_defaults();
}


ME::ME(const ME& matrix) :
    n(matrix.n), 
    jc(matrix.jc), 
    ir(matrix.ir), 
    Ax(matrix.Ax), 
    Az(matrix.Az)
{
    umfpack_defaults();
}


ME& ME::operator=(const ME& matrix)
{
    if (this == &matrix)
        return *this;
    umfpack_free();
    n = matrix.n;
    jc = matrix.jc;
    ir = matrix.ir;
    Ax = matrix.Ax;
    Az = matrix.Az;
    umfpack_defaults();
    return *this;
}


ME::~ME()
{
    umfpack_free();
}


int ME::factor()
{
    umfpack_free();
    int status;
    if (Az) {
        status = (umfpack_zi_symbolic(n, n, jc, ir, Ax, Az, 
                                      &Symbolic, Control, Info) != UMFPACK_OK ||
                  umfpack_zi_numeric(jc, ir, Ax, Az, Symbolic, 
                                     &Numeric, Control, Info) != UMFPACK_OK);
    } else {
        status = (umfpack_di_symbolic(n, n, jc, ir, Ax, 
                                      &Symbolic, Control, Info) != UMFPACK_OK ||
                  umfpack_di_numeric(jc, ir, Ax, Symbolic, 
                                     &Numeric, Control, Info) != UMFPACK_OK);
    }
    return status;
}


int ME::solve(double* xx, double* xz, double* bx, double* bz)
{
    int status = 0;
    if (!Numeric)
        status = factor();

    if (Az) {
        status = (status || 
                  umfpack_zi_solve(UMFPACK_A, jc, ir, Ax, Az, xx, xz, bx, bz,
                                   Numeric, Control, Info)) != UMFPACK_OK;
    } else {
        status = (status ||
                  umfpack_di_solve(UMFPACK_A, jc, ir, Ax, xx, bx, 
                                   Numeric, Control, Info) != UMFPACK_OK ||
                  umfpack_di_solve(UMFPACK_A, jc, ir, Ax, xz, bz, 
                                   Numeric, Control, Info) != UMFPACK_OK);
    }
    return status;
}


int ME::solve(dcomplex* x, dcomplex* b)
{
    vector<double> tmp(4*n);

    for (int i = 0; i < n; ++i) {
        tmp[i  ] = real(b[i]);
        tmp[i+n] = imag(b[i]);
    }

    int status = solve(&(tmp[2*n]), &(tmp[3*n]), &(tmp[0]), &(tmp[n]));

    for (int i = 0; i < n; ++i)
        x[i] = dcomplex(tmp[i+2*n], tmp[i+3*n]);
    return status;
}


int ME::solve(double* x, double* b)
{
    assert(Az == NULL);
    int status = 0;
    if (!Numeric)
        status = factor();

    status = (status || 
              umfpack_di_solve(UMFPACK_A, jc, ir, Ax, x, b, 
                               Numeric, Control, Info) != UMFPACK_OK);
    return status;
}


void ME::umfpack_defaults()
{
    if (Az)
        umfpack_zi_defaults(Control);
    else
        umfpack_di_defaults(Control);
    Numeric  = NULL;
    Symbolic = NULL;
}


void ME::umfpack_free()
{
    if (Az) {
        if (Numeric)    umfpack_zi_free_numeric(&Numeric);
        if (Symbolic)   umfpack_zi_free_symbolic(&Symbolic);
    } else {
        if (Numeric)    umfpack_di_free_numeric(&Numeric);
        if (Symbolic)   umfpack_di_free_symbolic(&Symbolic);
    }
}
