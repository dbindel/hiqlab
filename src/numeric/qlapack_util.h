#ifndef QLAPACK_UTIL_H
#define QLAPACK_UTIL_H

#include <complex>

/** ComplexDouble versions of LAPACK routines
 */

// -- General matrix multiply
int zgemm(char* trans, int m, int n, int k,
          std::complex<double> alpha,
          std::complex<double>* A, int ldA,
          std::complex<double>* B, int ldB,
          std::complex<double> beta,
          std::complex<double>* C, int ldC);


// -- General matrix vector multiply
int zgemv(char* trans, int m, int n,
          std::complex<double> alpha,
          std::complex<double>* A, int ldA,
          std::complex<double>* x, int incx,
          std::complex<double> beta,
          std::complex<double>* y, int incy);


// -- General solve
int zgesv(int n, int nrhs,
          std::complex<double>* A, int ldA,
          std::complex<double>* B, int ldB,
          int* ipiv);

          
// -- Compute generalized eigenvalue problem
int zggev(int n, std::complex<double>* A, int ldA,
                 std::complex<double>* B, int ldB,
                 std::complex<double>* alpha,
                 std::complex<double>* beta,
                 std::complex<double>* VL, int ldVL,
                 std::complex<double>* VR, int ldVR);


// -- Compute generalized Schur form
int zgges(int n, std::complex<double>* A, int ldA,
                 std::complex<double>* B, int ldB,
                 std::complex<double>* alpha,
                 std::complex<double>* beta,
                 std::complex<double>* VSL, int ldVSL,
                 std::complex<double>* VSR, int ldVSR);


// -- Reorder Schur form to move row at [pos_orig] to [pos_new]
int ztgexc(int n, std::complex<double>* A, int ldA,
                  std::complex<double>* B, int ldB,
                  std::complex<double>* Q, int ldQ,
                  std::complex<double>* Z, int ldZ,
                  int pos_orig, int pos_new);


// -- Reorder Schur form to move selected values flagged in [select]
int ztgsen(int n, std::complex<double>* A, int ldA,
                  std::complex<double>* B, int ldB,
                  std::complex<double>* Q, int ldQ,
                  std::complex<double>* Z, int ldZ,
                  std::complex<double>* alpha,
                  std::complex<double>* beta,
                  int* select);

#endif /* QLAPACK_UTIL_H */
