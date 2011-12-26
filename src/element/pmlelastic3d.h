/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELASTIC3D_H
#define PMLELASTIC3D_H

#include "pmlelement.h"
#include "qmatrix.h"


/** Three-dimensional elastic PML element.
 */
class PMLElastic3d : public PMLElement {
 public:

    /** Construct a new 3D isotropic elasticity element.
     *
     * @param E           Young's modulus
     * @param nu          Poisson ratio
     * @param rho         Mass density
     */
    PMLElastic3d(double E, double nu, double rho);

    /** Construct a new 3D isotropic elasticity element.
     *
     * @param D           Material property matrix (6-by-6)
     * @param rho         Mass density
     */
    PMLElastic3d(double* D, double rho);

    ~PMLElastic3d();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

    double* stress(Mesh* mesh, int eltid, double* X, double* stress);

 private:
    QMatrix1<double,6,6> D;
    double rho;
    int nen;

    template<class Ts, class T> inline
        void compute_stress(Ts* dNdx, T* u, T* sig);

    template<class T> inline void add_K(T* dNdx,   T W, QMatrix<T>& K);
    template<class T> inline void add_M(double* N, T W, QMatrix<T>& M);
    template<class T> inline void add_Ku(T* dNdx,   T* u, T W, QMatrix<T>& F);
    template<class T> inline void add_Ma(double* N, T* a, T W, QMatrix<T>& F);
};


#endif /* PMLELASTIC3D_H */
