/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELASTIC3D_TE_H
#define PMLELASTIC3D_TE_H

#include "pmlelement.h"
#include "qmatrix.h"

/** Three-dimensional thermoelastic PML element.
 */
class PMLElastic3d_te : public PMLElement {
 public:

    /** Construct a new 3D isotropic thermoelastic element.
     *
     * @param E             Young's modulus
     * @param nu            Poisson's ratio
     * @param rho           Mass density
     * @param at            Coefficient of thermal expansion
     * @param cp            Thermal capacity at constant pressure
     * @param kt            Thermal conductivity
     * @param T0            Reference temperature
     */

    PMLElastic3d_te(double E, double nu, double rho,
                     double at, double cp, double kt, double T0);

    /** Construct a new 3D thermoelastic element.
     *
     * @param Db            Material property matrix (7-by-7)
     * @param rho           Mass density
     * @param at            Coefficient of thermal expansion
     * @param cp            Thermal capacity at constant pressure
     * @param kt            Thermal conductivity
     * @param T0            Reference temperature
     */

    PMLElastic3d_te(double* Db, double rho,
                     double at, double cp, double kt, double T0);

    ~PMLElastic3d_te();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double D[36],vp[6];
    double rho;
    double at,cp,kt,T0;
    int nen;

    void split_Db(double* Db);

    template<class T> inline 
        void add_K(double* N, T* dNdx, T W, QMatrix<T>& K);
    template<class T> inline 
        void add_M(double* N,          T W, QMatrix<T>& M);
    template<class T> inline 
        void add_C(double* N, T* dNdx, T W, QMatrix<T>& C);

    template<class T> inline 
        void add_Ku(double* N, T* dNdx, T* u, T W, QMatrix<T>& F);
    template<class T> inline 
        void add_Ma(double* N,          T* a, T W, QMatrix<T>& F);
    template<class T> inline 
        void add_Cv(double* N, T* dNdx, T* v, T W, QMatrix<T>& F);
};

#endif /* PMLELASTIC3D_TE_H */
