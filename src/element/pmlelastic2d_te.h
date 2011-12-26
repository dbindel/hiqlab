/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELASTIC2D_TE_H
#define PMLELASTIC2D_TE_H

#include "pmlelement.h"
#include "qmatrix.h"

/** Two-dimensional thermoelastic PML element.
 */
class PMLElastic2d_te : public PMLElement {
 public:

    /** Construct a new plane stress/plane strain isotropic
     *  thermoelastic element.
     *
     * @param E             Young's modulus
     * @param nu            Poisson's ratio
     * @param rho           Mass density
     * @param at            Coefficient of thermal expansion
     * @param cp            Thermal capacity at constant pressure
     * @param kt            Thermal conductivity
     * @param T0            Reference temperature
     * @param plane_type    [0] for plane strain, [1] for plane stress
     */

    PMLElastic2d_te(double E, double nu, double rho,
                    double at, double cp, double kt, double T0,
                    int plane_type);

    /** Construct a new plane plane strain thermoelastic element.
     * @param Db            Material property matrix (4-by-4)
     * @param rho           Mass density
     * @param at            Coefficient of thermal expansion
     * @param cp            Thermal capacity at constant pressure
     * @param kt            Thermal conductivity
     * @param T0            Reference temperature
     */

    PMLElastic2d_te(double* Db, double rho,
                    double at, double cp, double kt, double T0);

    ~PMLElastic2d_te();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double D[9],vp[3],vCooiv;
    double rho;
    double at,cp,kt,T0;
    int nen;

    void split_Db(double* Db);

    template<class T> inline 
        void add_K(double* N, T* dNdx,   T W, QMatrix<T>& K);
    template<class T> inline 
        void add_M(double* N,            T W, QMatrix<T>& M);
    template<class T> inline 
        void add_C(double* N, T* dNdx,   T W, QMatrix<T>& C);

    template<class T> inline 
        void add_Ku(double* N, T* dNdx,   T* u, T W, QMatrix<T>& F);
    template<class T> inline 
        void add_Ma(double* N,            T* a, T W, QMatrix<T>& F);
    template<class T> inline 
        void add_Cv(double* N, T* dNdx,   T* v, T W, QMatrix<T>& F);
};

#endif /* PMLELASTIC2D_TE_H */
