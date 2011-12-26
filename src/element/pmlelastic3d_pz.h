/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELASTIC3D_PZ_H
#define PMLELASTIC3D_PZ_H

#include "pmlelement.h"
#include "qmatrix.h"

/** Three-dimensional piezoelectric-elastic PML element.
 */
class PMLElastic3d_pz : public PMLElement {
 public:

    /** Construct a new 3D isotropic piezoelectric-elastic element.
     *  (Assuming mechanical and dielectric isotropy)
     *
     * @param E             Young's modulus
     * @param nu            Poisson's ratio
     * @param rho           Mass density
     * @param pz            Piezoelectric coefficients(6-by-3)
     * @param kds           Dielectric coefficient(constant stress)
     *                                            (constant strain -> kde)
     */

    PMLElastic3d_pz(double E, double nu, double rho,
                     double* pz, double kds);

    /** Construct a new 3D piezo-electric element.
     *
     * @param Db            Material property matrix ( 9-by- 9)
     *                      [ sigma        ] = [   D        -D*pzt] [epsilon]
     *                      [ electric disp] = [pz*D  kds-pz*D*pzt] [Efield ]
     * @param pz            Piezoelectric coefficients(6-by-3)
     * @param kds           Dielectric coefficient(constant stress)(3-by-3)
     *                                            (constant strain -> kde)
     * @param rho           Mass density
     */

    PMLElastic3d_pz(double* Db, double rho);

    ~PMLElastic3d_pz();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double D[36],Dpzt[18],kde[9];
    double rho;
    int nen;

    void split_Db(double* Db);

    template<class T> inline void 
        add_K(double* N, T* dNdx, T W, QMatrix<T>& K);
    template<class T> inline void 
        add_M(double* N,          T W, QMatrix<T>& M);

    template<class T> inline void 
        add_Ku(double* N, T* dNdx, T* u, T W, QMatrix<T>& F);
    template<class T> inline void 
        add_Ma(double* N,          T* a, T W, QMatrix<T>& F);
};

#endif /* PMLELASTIC3D_PZ_H */
