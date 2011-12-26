#ifndef UMFMATRIX_H
#define UMFMATRIX_H

extern "C" {
#include "umfpack.h"
}

#include "cscmatrix.h"


/** UMFPACK factorization of a compressed sparse column matrix
 */
class UMFMatrix {
public:

    UMFMatrix(CSCMatrix*);
    UMFMatrix(int n, int* jc, int* ir, double* Ax, double* Az);
    UMFMatrix(const UMFMatrix&);
    ~UMFMatrix();

    UMFMatrix& operator=(const UMFMatrix& matrix);

    int get_n()   { return n;         }
    int is_real() { return (Az == 0); }

    /** Factor the matrix */
    int factor();

    /** Solve a linear system A*x = b. */
    int solve(double* xx, double* xz, double* bx, double* bz);
    int solve(dcomplex* x, dcomplex* b);
    int solve(double* x, double* b);

    double& umf_control(int i) { return Control[i]; }
    double& umf_info(int i)    { return Info[i];    }

 private:

    int n;
    int* jc;
    int* ir;
    double* Ax;
    double* Az;

    void* Symbolic; /* UMFPACK structs */
    void* Numeric;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];

    void umfpack_defaults();
    void umfpack_free();
};

#endif /* UMFMATRIX_H */
