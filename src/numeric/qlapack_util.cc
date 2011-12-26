#include "qlapack_util.h"
#include <vector>

// -- Declarations of the LAPACK functions
using std::complex;

extern "C"
int zgemm_(const char& transA,
           const char& transB,
           const int& M, const int& N, const int& K,
           const complex<double>& alpha,
           const complex<double>* A, const int& ldA,
           const complex<double>* B, const int& ldB,
           const complex<double>& beta,
           complex<double>* C, const int& ldC);


extern "C"
int zgemv_(const char& trans, const int& M, const int& N,
           const complex<double>& alpha,
           const complex<double>* A, const int& ldA,
           const complex<double>* x, const int& incx,
           const complex<double>& beta,
           complex<double>* y, const int& incy);


extern "C"
int zgesv_(const int& N, const int& nrhs,
           complex<double>* A, const int& ldA, int* ipiv,
           complex<double>* B, const int& ldB, int& info);


extern "C"
int zggev_(const char& jobVL,
           const char& jobVR,
           const int& N,
           const complex<double>* A, const int& ldA,
           const complex<double>* B, const int& ldB,
           const complex<double>* alpha,
           const complex<double>* beta,
           complex<double>* VL, const int& ldVL,
           complex<double>* VR, const int& ldVR,
           complex<double>* work, const int& lwork, double* rwork,
           int& info);


extern "C"
int zgges_(const char& jobVSL, const char& jobVSR,
           const char& sortc, int(*selectg)(complex<double>,complex<double>),
           const int& N, complex<double>* A, const int& ldA,
                         complex<double>* B, const int& ldB,
           int* sdim, complex<double>* alpha, complex<double>* beta,
           complex<double>* VSL, const int& ldVSL,
           complex<double>* VSR, const int& ldVSR,
           complex<double>* work, const int& lwork, double* rwork,
           int* bwork, int& info);


extern "C"
int ztgexc_(int& wantq, int& wantz, const int& N,
           complex<double>* A, const int& ldA,
           complex<double>* B, const int& ldB,
           complex<double>* Q, const int& ldQ,
           complex<double>* Z, const int& ldZ, int& ifst, int& ilst,
           int& info);


extern "C"
int ztgsen_(int& ijob, int& wantq, int& wantz, const int* select,
           const int& N,
           complex<double>* A, const int& ldA,
           complex<double>* B, const int& ldB,
           complex<double>* alpha, complex<double>* beta,
           complex<double>* Q, const int& ldQ,
           complex<double>* Z, const int& ldZ,
           int& m, double* pl, double* pr, double* dif,
           complex<double>* work, const int& lwork,
           int* iwork, const int& liwork, int& info);


// -- Wrapped C function calls

int zgemm(char* trans, int m, int n, int k,
          complex<double> alpha,
          complex<double>* A, int ldA,
          complex<double>* B, int ldB,
          complex<double> beta,
          complex<double>* C, int ldC)
{
    char transA = trans[0];
    char transB = trans[1];

    int info = zgemm_(transA, transB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);

    return info;
}


int zgemv(char* trans, int m, int n,
          complex<double> alpha,
          complex<double>* A, int ldA,
          complex<double>* x, int incx,
          complex<double> beta,
          complex<double>* y, int incy)
{
    int info = zgemv_(*trans, m, n, alpha, A, ldA, x, incx, beta, y, incy);
    return info;
}



int zgesv(int n, int nrhs,
          complex<double>* A, int ldA,
          complex<double>* B, int ldB,
          int* ipiv)
{
    int info = 0;
    zgesv_(n, nrhs, A, ldA, ipiv, B, ldB, info);
    return info;
}


int zggev(int n, complex<double>* A, int ldA,
                 complex<double>* B, int ldB,
                 complex<double>* alpha,
                 complex<double>* beta,
                 complex<double>* VL, int ldVL,
                 complex<double>* VR, int ldVR)
{
    // -- work arrays
    complex<double>* work;
    complex<double>  work1(0,0);
    int             lwork;
    double*         rwork = new double[8*n];
    int              info;

    // -- copy A,B since zggev_ overwrites
    complex<double>* Ac = new complex<double>[n*n];
    complex<double>* Bc = new complex<double>[n*n];
    int ldAc = n;
    int ldBc = n;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            Ac[ldAc*j+i] = A[ldA*j+i];
            Bc[ldBc*j+i] = B[ldB*j+i];
        }

    // -- job options
    char jobVL = (VL) ? 'V' : 'N';
    char jobVR = (VR) ? 'V' : 'N';
    int  ldVLc = (ldVL > 1) ? ldVL : 1;
    int  ldVRc = (ldVR > 1) ? ldVR : 1;

    // -- query work array
    lwork = -1;
    zggev_(jobVL, jobVR, n, Ac, ldAc, Bc, ldBc, alpha, beta,
                           VL, ldVLc, VR, ldVRc,
                           &work1, lwork, rwork, info);

    // -- compute
    lwork = (int)std::real(work1);
    work  = new complex<double>[lwork];
    zggev_(jobVL, jobVR, n, Ac, ldAc, Bc, ldBc, alpha, beta,
                           VL, ldVLc, VR, ldVRc,
                           work, lwork, rwork, info);

    // -- clean up
    delete[] rwork;
    delete[] work;
    delete[] Ac;
    delete[] Bc;

    return info;
}


int zgges(int n, complex<double>* A, int ldA,
                 complex<double>* B, int ldB,
                 complex<double>* alpha,
                 complex<double>* beta,
                 complex<double>* VSL, int ldVSL,
                 complex<double>* VSR, int ldVSR)
{
    // -- work arrays
    complex<double>* work;
    complex<double>  work1(0,0);
    int             lwork;
    double*         rwork = new double[8*n];
    int              info;


    // -- job options
    char sortc  = 'N';
    int  sdim   = 0;
    char jobVSL = (VSL) ? 'V' : 'N';
    char jobVSR = (VSR) ? 'V' : 'N';
    int  ldVSLc = (ldVSL > 1) ? ldVSL : 1;
    int  ldVSRc = (ldVSR > 1) ? ldVSR : 1;

    // -- query work array
    lwork = -1;
    zgges_(jobVSL, jobVSR, sortc, NULL, 
           n, A, ldA, B, ldB, &sdim, alpha, beta,
                           VSL, ldVSLc, VSR, ldVSRc,
                           &work1, lwork, rwork, NULL, info);

    // -- compute
    lwork = (int)std::real(work1);
    work  = new complex<double>[lwork];
    zgges_(jobVSL, jobVSR, sortc, NULL, 
           n, A, ldA, B, ldB, &sdim, alpha, beta,
                           VSL, ldVSLc, VSR, ldVSRc,
                           work, lwork, rwork, NULL, info);

    // -- clean up
    delete[] rwork;
    delete[] work;

    return info;
}


int ztgexc(int n, complex<double>* A, int ldA,
                  complex<double>* B, int ldB,
                  complex<double>* Q, int ldQ,
                  complex<double>* Z, int ldZ,
                  int pos_orig, int pos_new)
{
    int info = 0;

    // -- FORTRAN is 1 based
    pos_orig++;
    pos_new++;

    // -- job options
    int wantq = (Q) ? 1 : 0;
    int wantz = (Z) ? 1 : 0;
    int ldQc  = (ldQ > 1) ? ldQ : 1;
    int ldZc  = (ldZ > 1) ? ldZ : 1;

    // -- compute
    ztgexc_(wantq, wantz, n, A, ldA, B, ldB, Q, ldQc, Z, ldZc, pos_orig, pos_new, info);

    return info;
}


int ztgsen(int n, complex<double>* A, int ldA,
                  complex<double>* B, int ldB,
                  complex<double>* Q, int ldQ,
                  complex<double>* Z, int ldZ,
                  complex<double>* alpha,
                  complex<double>* beta,
                  int* select)
{
    int info;

    // -- work arrays
    complex<double>* work;
    complex<double>  work1(0,0);
    int             lwork;
    int*            iwork;
    int             iwork1;
    int             liwork;
    
    // -- job options
    int ijob  = 0;
    int wantq = (Q) ? 1 : 0;
    int wantz = (Z) ? 1 : 0;
    int ldQc  = (ldQ > 1) ? ldQ : 1;
    int ldZc  = (ldZ > 1) ? ldZ : 1;
    int msize;

    // -- query work array
    lwork = -1;
    liwork=  1;
    ztgsen_(ijob, wantq, wantz, select, n, A, ldA, B, ldB,
            alpha, beta, Q, ldQc, Z, ldZc, msize, NULL, NULL, NULL,
            &work1, lwork, &iwork1, liwork, info);

    // -- compute

    //FIXME: Due to a bug in memory allocation in ztgsen_ in LAPACK v3.0
    //       work must at LEAST be 2 ints. This is fixed if one is linking to
    //       LAPACK v3.1.1. Since ACML may have v3.0, to be safe we add 1,
    //       so that in total lwork>=2.
    lwork = (int)std::real(work1);
    lwork++;

    liwork= iwork1;
    work  = new complex<double>[lwork];
    iwork = new int[liwork];
    ztgsen_(ijob, wantq, wantz, select, n, A, ldA, B, ldB,
            alpha, beta, Q, ldQc, Z, ldZc, msize, NULL, NULL, NULL,
            work, lwork, iwork, liwork, info);

    // -- clean up
    delete[] work;
    delete[] iwork;

    return info;
}
