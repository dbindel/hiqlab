/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELASTIC2D_PZ_H
#define PMLELASTIC2D_PZ_H

#include "pmlelement.h"
#include "qmatrix.h"

/** Two-dimensional piezoelectric-elastic PML element.
 */
class PMLElastic2d_pz : public PMLElement {
 public:

    /** Construct a new plane strain/stress isotropic piezoelectric-elastic
     *  element.(Assuming mechanical and dielectric isotropy)
     *
     * @param E             Young's modulus
     * @param nu            Poisson's ratio
     * @param rho           Mass density
     * @param pz            Piezoelectric coefficients(6-by-3)
     * @param kds           Dielectric coefficient(constant stress)
     *                                            (constant strain -> kde)
     * @param plane_type    [0] for plane strain, [1] for plane stress
     */

    PMLElastic2d_pz(double E, double nu, double rho,
                    double* pz, double kds, int plane_type);

    /** Construct a new 2D plane piezoelectric elastic element.
     *
     * @param Db            Material property matrix ( 5-by- 5)
     *                      [ sigma        ] = [   D        -D*pzt] [epsilon]
     *                      [ electric disp] = [pz*D  kds-pz*D*pzt] [Efield ]
     * @param rho           Mass density
     */
    PMLElastic2d_pz(double* Db, double rho);

    ~PMLElastic2d_pz();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double D[9], Dpzt[6], kde[4];
    double rho;
    int nen;

    void split_Db(double* Db);

    template<class T> inline void 
        add_K(double* N, T* dNdx,   T W, QMatrix<T>& K);
    template<class T> inline void 
        add_M(double* N,            T W, QMatrix<T>& M);

    template<class T> inline void 
        add_Ku(double* N, T* dNdx,   T* u, T W, QMatrix<T>& F);
    template<class T> inline void 
        add_Ma(double* N,            T* a, T W, QMatrix<T>& F);
};

#endif /* PMLELASTIC2D_PZ_H */
