/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLSCALAR21_H
#define PMLSCALAR21_H

#include "pmlelement.h"
#include "qmatrix.h"


/** One-dimensional scalar PML element.
 */
class PMLScalar1d : public PMLElement {
 public:

    /** Construct a new isotropic scalar wave element.
     *
     * @param kappa  Constitutive coefficent
     * @param rho    Mass density
     */
    PMLScalar1d(double kappa, double rho);
    ~PMLScalar1d();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

    double* stress(Mesh* mesh, int eltid, double* X, double* stress);

 private:
    double D;
    double rho;
    int nen;

    template<class T> inline void add_K(T* dNdx,   T W, QMatrix<T>& K);
    template<class T> inline void add_M(double* N, T W, QMatrix<T>& M);

    template<class T> inline void add_Ku(T* dNdx,   T* u, T W, QMatrix<T>& F);
    template<class T> inline void add_Ma(double* N, T* a, T W, QMatrix<T>& F);
};


#endif /* PMLSCALAR1D_H */
