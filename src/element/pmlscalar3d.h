/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLSCALAR3D_H
#define PMLSCALAR3D_H

#include "pmlelement.h"
#include "qmatrix.h"


/** Three-dimensional scalar PML element.
 */
class PMLScalar3d : public PMLElement {
 public:

    /** Construct a new isotropic scalar wave element.
     *
     * @param kappa  Constitutive coefficent
     * @param rho    Mass density
     */
    PMLScalar3d(double kappa, double rho);

    /** Construct a new 2D element type
     *
     * @param D    Material property matrix (2-by-2)
     * @param rho  Mass density
     */
    PMLScalar3d(double* D, double rho);

    ~PMLScalar3d();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

    double* stress(Mesh* mesh, int eltid, double* X, double* stress);

 private:
    double D[9];
    double rho;
    int nen;

    template<class T> inline T dotD(T* a, T* b);

    template<class T> inline void add_K(T* dNdx,   T W, QMatrix<T>& K);
    template<class T> inline void add_M(double* N, T W, QMatrix<T>& M);

    template<class T> inline void add_Ku(T* dNdx,   T* u, T W, QMatrix<T>& F);
    template<class T> inline void add_Ma(double* N, T* a, T W, QMatrix<T>& F);
};


#endif /* PMLSCALAR3D_H */
