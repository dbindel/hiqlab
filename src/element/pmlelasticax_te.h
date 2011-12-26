/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELASTICAX_TE_H
#define PMLELASTICAX_TE_H

#include "pmlelement.h"
#include "qmatrix.h"

/** Axisymmetric thermoelastic PML element.
 */
class PMLElasticAxis_te : public PMLElement {
 public:

    /** Construct a new axisymmetric isotropic thermoelastic element.
     *
     * @param E             Youngs' modulus
     * @param nu            Poisson's ratio
     * @param rho           Mass density
     * @param at            Coefficient of thermal expansion
     * @param cp            Thermal capacity at constant pressure
     * @param kt            Thermal conductivity
     * @param T0            Reference temperature
     */
    PMLElasticAxis_te(double E, double nu, double rho,
                      double at, double cp, double kt, double T0);

    /** Construct a new axisymmetric thermoelastic element.
     *
     * @param D             Material property matrix (5-by-5)
     * @param rho           Mass density
     * @param at            Coefficient of thermal expansion
     * @param cp            Thermal capacity at constant pressure
     * @param kt            Thermal conductivity
     * @param T0            Reference temperature
     */
    PMLElasticAxis_te(double* D, double rho,
                      double at, double cp, double kt, double T0);

    ~PMLElasticAxis_te();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double D[16],vp[4];
    double rho;
    double at,cp,kt,T0;
    int nen;

    void split_Db(double* Db);

    template<class T> inline 
        void add_K(T* dNdx, double* N, double r, T W, QMatrix<T>& K);
    template<class T> inline 
        void add_M(         double* N,           T W, QMatrix<T>& M);
    template<class T> inline 
        void add_C(T* dNdx, double* N, double r, T W, QMatrix<T>& K);

    template<class T> inline 
        void add_Ku(T* dNdx, double* N, double r, T* u, T W, QMatrix<T>& F);
    template<class T> inline
        void add_Ma(         double* N,           T* a, T W, QMatrix<T>& F);
    template<class T> inline 
        void add_Cv(T* dNdx, double* N, double r, T* v, T W, QMatrix<T>& F);
};

#endif /* PMLElasticAX_TE_H */
