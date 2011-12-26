/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELASTICTAX_H
#define PMLELASTICTAX_H

#include "pmlelement.h"
#include "qmatrix.h"


/** Generalized axisymmetric elastic PML element.
 */
class PMLElasticTAxis : public PMLElement {
 public:

    /** Construct a new axisymmetric isotropic elasticity element.
     *
     * @param E    Young's modulus
     * @param nu   Poisson ratio
     * @param rho  Mass density
     * @param l    Radial wave number
     */
    PMLElasticTAxis(double E, double nu, double rho, int l);

    /** Construct a new axisymmetric elasticity element type.
     *
     * @param D    Material property matrix (6-by-6)
     * @param rho  Mass density
     * @param l    Radial wave number
     */
    PMLElasticTAxis(double* D, double rho, int l);

    ~PMLElasticTAxis();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

    void project_L2(Mesh* mesh, int eltid, int nfields,
                    double* Mdiag, double* xfields,
                    FieldEval& xfunc);

    double* mean_power(Mesh* mesh, int eltid, double* X, double* EX);
    double* stress(Mesh* mesh, int eltid, double* X, double* stress);

 private:
    QMatrix1<double,6,6> D;
    double rho;
    int ltheta;
    int nen;

    template<class T> inline void
        compute_DB (T* dNdx, double Ndr,
                    QMatrix1<T,6,3>& DB1, QMatrix1<T,6,3>& DB2);
    template<class T> inline void
        compute_BDB(T* dNdx, double Ndr,
                    QMatrix1<T,6,3>& DB1, QMatrix1<T,6,3>& DB2,
                    QMatrix1<T,3,3>& BDB);

    template<class T> inline void add_BDBu   (T* dNdx, double Ndr,
                                              T* DB1, T* DB2, T Wt, T* BDB);
    template<class T> inline void add_K(T* dNdx, double* N, double r,
                                        T W, QMatrix<T>& K);
    template<class T> inline void add_M(double* N, T W, QMatrix<T>& M);

    template<class T> inline void add_Ku(T* dNdx, double* N, double r,
                                         T* u, T W, QMatrix<T>& F);
    template<class T> inline void add_Ma(double* N, T* a, T W, QMatrix<T>& F);
};

#endif /* PMLELASTICTAX_H */
