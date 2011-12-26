/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELASTICAX_H
#define PMLELASTICAX_H

#include "pmlelement.h"
#include "qmatrix.h"


/** Axisymmetric elastic PML element.
 */
class PMLElasticAxis : public PMLElement {
 public:

    /** Construct a new axisymmetric isotropic elasticity element.
     *
     * @param E    Young's modulus
     * @param nu   Poisson ratio
     * @param rho  Mass density
     */
    PMLElasticAxis(double E, double nu, double rho);

    /** Construct a new axisymmetric elasticity element type.
     *
     * @param D    Material property matrix (4-by-4)
     * @param rho  Mass density
     */
    PMLElasticAxis(double* D, double rho);

    ~PMLElasticAxis();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

    void project_L2(Mesh* mesh, int eltid, int nfields,
                    double* Mdiag, double* xfields,
                    FieldEval& xfunc);

    double* mean_power(Mesh* mesh, int eltid, double* X, double* EX);
    double* stress(Mesh* mesh, int eltid, double* X, double* stress);

 private:
    QMatrix1<double,4,4> D;
    double rho;
    int nen;

    template<class Ts, class T> inline
        void compute_stress(Ts* dNdx, double* N, double r, T* u, T* sig);

    template<class T> inline void add_K(T* dNdx, double* N, double r,
                                        T W, QMatrix<T>& K);
    template<class T> inline void add_M(double* N, T W, QMatrix<T>& M);

    template<class T> inline void add_Ku(T* dNdx, double* N, double r,
                                         T* u, T W, QMatrix<T>& F);
    template<class T> inline void add_Ma(double* N, T* a, T W, QMatrix<T>& F);
};

#endif /* PMLELASTICAX_H */
