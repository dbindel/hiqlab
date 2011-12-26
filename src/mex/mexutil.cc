/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <mex.h>
#include <algorithm>
#include "cscmatrix.h"


mxArray* CSCMatrix_to_mx(CSCMatrix* K)
{
    using std::copy;
    int n   = K->get_n();
    int nnz = K->get_nnz();
    mxArray* Kmat = 
        mxCreateSparse(n, n, nnz, K->get_Az() ? mxCOMPLEX : mxREAL);
    copy(mxGetJc(Kmat), mxGetJc(Kmat)+n+1, K->get_jc());
    copy(mxGetIr(Kmat), mxGetIr(Kmat)+nnz, K->get_ir());
    copy(mxGetPr(Kmat), mxGetPr(Kmat)+nnz, K->get_Ax());
    if (K->get_Az())
        copy(mxGetPi(Kmat), mxGetPi(Kmat)+nnz, K->get_Az());
    return Kmat;
}


CSCMatrix* mx_to_CSCMatrix(const mxArray* Kmat)
{
    using std::copy;
    int n   = mxGetN(Kmat);
    int nnz = mxGetJc(Kmat)[n];
    CSCMatrix* K = new CSCMatrix(n, n, nnz, mxGetPi(Kmat) == NULL);
    copy(K->get_jc(), K->get_jc()+n+1, mxGetJc(Kmat));
    copy(K->get_ir(), K->get_ir()+nnz, mxGetIr(Kmat));
    copy(K->get_Ax(), K->get_Ax()+nnz, mxGetPr(Kmat));
    if (K->get_Az())
        copy(K->get_Az(), K->get_Az()+nnz, mxGetPi(Kmat));
    return K;
}



