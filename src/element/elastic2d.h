/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef ELASTIC2D_H
#define ELASTIC2D_H

#include "element.h"
#include "qmatrix.h"


/** Two-dimensional elastic element.
 */
class Elastic2d : public Element {
 public:

    /** Construct a new plane stress/plane strain isotropic elasticity element.
     *
     * @param E           Young's modulus
     * @param nu          Poisson ratio
     * @param rho         Mass density
     * @param plane_type  Zero for plane strain, one for plane stress
     */
    Elastic2d(double E, double nu, double rho, int plane_type);

    /** Construct a new 2D element type
     *
     * @param D    Material property matrix (3-by-3)
     * @param rho  Mass density
     */
    Elastic2d(double* D, double rho);

    ~Elastic2d();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

    void project_L2(Mesh* mesh, int eltid, int nfields,
                    double* Mdiag, double* xfields,
                    FieldEval& xfunc);

    double* mean_power(Mesh* mesh, int eltid, double* X, double* EX);
    double* stress(Mesh* mesh, int eltid, double* X, double* stress);

 private:
    QMatrix1<double,3,3> D;
    double rho;
    int nen;

    void compute_DB (double* dNdx, QMatrix1<double,3,2>& DB);
    void compute_BDB(double* dNdx, QMatrix1<double,3,2>& DB,
                     QMatrix1<double,2,2>& BDB);
    void add_BDBu   (double* dNdx, double* DB, double Wt, double* BDB);
    void add_K(double* dNdx,   double W, QMatrix<double>& K);
    void add_M(double* N, double W, QMatrix<double>& M);
    void add_Ku(double* dNdx,   double* u, double W, QMatrix<double>& F);
    void add_Ma(double* N, double* a, double W, QMatrix<double>& F);
};


#endif /* ELASTIC2D_H */
