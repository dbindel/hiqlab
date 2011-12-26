/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELASTIC2HD_PZ_H
#define PMLELASTIC2HD_PZ_H

#include "pmlelement.h"
#include "qmatrix.h"

/** Two-and-a-half-dimensional piezoelectric-elastic PML element.
 */
class PMLElastic2hd_pz : public PMLElement {
 public:

    /** Construct a new plane stress isotropic piezoelectric-elastic
     *  element.  (Assuming mechanical and dielectric isotropy)
     *
     * @param E             Young's modulus
     * @param nu            Poisson's ratio
     * @param rho           Mass density
     * @param pz            Piezoelectric coefficients(6-by-3)
     * @param kds           Dielectric coefficient(constant stress)
     *                                            (constant strain -> kde)
     */

    PMLElastic2hd_pz(double E, double nu, double rho,
                    double* pz, double kds);

    /** Construct a new 2D plane stress piezoelectric elastic element.
     *
     * @param Db            Material property matrix ( 4-by- 4)
     *                      [ sigma        ] = [   D        -D*pzt] [epsilon]
     *                      [ electric disp] = [pz*D  kds-pz*D*pzt] [Efield ]
     * @param rho           Mass density
     */

    PMLElastic2hd_pz(double* Db, double rho);

    ~PMLElastic2hd_pz();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double D[9],Dpzt[3];
    double rho;
    int nen;

    void split_Db(double* Db);

    template<class T> inline void compute_DB (T* dNdx, T* DB              );
    template<class T> inline void compute_BDB(T* dNdx, T* DB,       T* BDB);
    template<class T> inline void add_BDBu   (T* dNdx, T* DB, T Wt, T* BDB);

    template<class T> inline void compute_BDpztN(T* dNdx, T* Dpztb, T* BDpztb);


    template<class T> inline void add_K_pz(double* N, T* dNdx,   T W, QMatrix<T>& K);
    template<class T> inline void add_M_pz(double* N,            T W, QMatrix<T>& M);

    template<class T> inline void add_Ku_pz(double* N, T* dNdx,   T* u, T W, QMatrix<T>& F);
    template<class T> inline void add_Ma_pz(double* N,            T* a, T W, QMatrix<T>& F);
};

#endif /* PMLELASTIC2HD_PZ_H */
