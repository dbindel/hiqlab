/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELASTIC2D_H
#define PMLELASTIC2D_H

#include "pmlelement.h"
#include "qmatrix.h"


/** Two-dimensional elastic PML element.
 */
class PMLElastic2d : public PMLElement {
 public:

    /** Construct a new plane stress/plane strain isotropic elasticity element.
     *
     * @param E           Young's modulus
     * @param nu          Poisson ratio
     * @param rho         Mass density
     * @param plane_type  Zero for plane strain, one for plane stress
     */
    PMLElastic2d(double E, double nu, double rho, int plane_type);

    /** Construct a new 2D element type
     *
     * @param D    Material property matrix (3-by-3)
     * @param rho  Mass density
     */
    PMLElastic2d(double* D, double rho);

    /** Construct an element in which the properties are given by
     *  a Lua function call.
     *
     * @param L    Lua state
     * @param f    Function handle
     */
    PMLElastic2d(lua_State* L, int f);

    ~PMLElastic2d();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

    void project_L2(Mesh* mesh, int eltid, int nfields,
                    double* Mdiag, double* xfields,
                    FieldEval& xfunc);

    double* mean_power(Mesh* mesh, int eltid, double* X, double* EX);
    double* stress(Mesh* mesh, int eltid, double* X, double* stress);

 private:
    lua_State* L;
    QMatrix1<double,3,3> D;
    double rho;
    int nen;

    void compute_D(Mesh* mesh, int eltid, double* N);

    template<class Ts, class T> inline void
        compute_stress(Ts* dNdx, T* u, T* sig);

    template<class T> inline void add_K(T* dNdx,   T W, QMatrix<T>& K);
    template<class T> inline void add_M(double* N, T W, QMatrix<T>& M);
    template<class T> inline void add_Ku(T* dNdx,   T* u, T W, QMatrix<T>& F);
    template<class T> inline void add_Ma(double* N, T* a, T W, QMatrix<T>& F);
};


#endif /* PMLELASTIC2D_H */
