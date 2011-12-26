/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "qmatrix.h"
#include "gaussquad.h"
#include "qcomplex.h"
#include "shapes.h"
#include "pmlelastic2hd_pz.h"
#include "material_model.h"
#include "mesh.h"

#define ME          PMLElastic2hd_pz
#define MAXNEN      25


inline void Bmatrix() 
{
    /* <generator matexpr>
      function Bmatrix(dN) = [dN(1), 0;
                                  0, dN(2);
                              dN(2), dN(1)];
    */
}


ME::ME(double E, double nu, double rho, double* pz, double kds):
    PMLElement(3), rho(rho)
{
    double lambda = E*nu/(1+nu)/(1-2*nu);
    double mu     = E/2/(1+nu);
    double Db[16];

    // -- Plane Stress Isotropic
    piezo_elasticity_2HD_stress(Db, lambda, mu, pz, kds);
    split_Db(Db);
}


ME::ME(double* Db, double rho):
    PMLElement(3), rho(rho)
{
    split_Db(Db);
}


ME::~ME()
{
}


void ME::split_Db(double* Db)
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            D[j+3*i] = Db[j+4*i];
        Dpzt[i] = -Db[3+4*i];
    }
}


template<class T> inline
void ME::compute_DB(T* dNdx, T* DB)
{
    for (int k = 0; k < 3; ++k) {
        DB[k  ] = D[k  ]*dNdx[0] + D[k+6]*dNdx[1];
        DB[k+3] = D[k+3]*dNdx[1] + D[k+6]*dNdx[0];
    }
}


template<class T> inline
void ME::compute_BDB(T* dNdx, T* DB, T* BDB)
{
    for (int i = 0; i < 2; ++i) {
        BDB[0] = dNdx[0]*DB[0] + dNdx[1]*DB[2];
        BDB[1] = dNdx[1]*DB[1] + dNdx[0]*DB[2];
        DB += 3; BDB += 2;
    }
}


template<class T> inline
void ME::add_BDBu(T* dNdx, T* DB, T Wt, T* BDB)
{
    BDB[0] += (dNdx[0]*DB[0] + dNdx[1]*DB[2])*Wt;
    BDB[1] += (dNdx[1]*DB[1] + dNdx[0]*DB[2])*Wt;
}


template <class T> inline
void ME::compute_BDpztN(T* dNdx, T* DpztN, T* BDpztN)
{
    BDpztN[0] = (dNdx[0]*DpztN[0] + dNdx[1]*DpztN[2]);
    BDpztN[1] = (dNdx[1]*DpztN[1] + dNdx[0]*DpztN[2]);
}


template<class T> inline
void ME::add_K_pz(double* N, T* dNdx, T W, QMatrix<T>& K)
{
    // -- Mechanical part of the stiffness
    // --   Form Kmm += B'*D*B * W
    // -- Coupling term of the stiffness
    // --   Form KmE -= B' * Dpzt * N * W

    T DBj[6], DpztN[3];
    T   A[4],  Ac[2];
    for (int j = 0; j < nen; ++j) {

        // -- Mechanical
        compute_DB(dNdx+2*j, DBj);
        // -- Coupling
        DpztN[0] = Dpzt[0] * N[j];
        DpztN[1] = Dpzt[1] * N[j];
        DpztN[2] = Dpzt[2] * N[j];

        for (int i = 0; i < nen; ++i) {
            // -- Mechanical
            compute_BDB(dNdx+2*i, DBj, A);
            for (int k = 0; k < 2; ++k)
                for (int l = 0; l < 2; ++l)
                    K(i*3+k, j*3+l) += A[k+2*l] * W;
            // -- Coupling
            compute_BDpztN(dNdx+2*i, DpztN, Ac);
            for (int k = 0; k < 2; ++k)
                    K(i*3+k, j*3+2) -= Ac[k] * W;

        }
    }
}


template<class T> inline
void ME::add_M_pz(double* N, T W, QMatrix<T>& K)
{
    // -- Form M += rho*kron(N*N', I2) * W
    for (int j = 0; j < nen; ++j) {
        for (int i = 0; i < nen; ++i) {
            T Mij = rho*N[i]*N[j]*W;
            for (int k = 0; k < 2; ++k)
                K(i*3+k, j*3+k) += Mij;
        }
    }
}


template<class T> inline
void ME::add_Ku_pz(double* N, T* dNdx, T* u, T W, QMatrix<T>& F)
{
    // -- Mechanical part
    // -- Form F += B'*(D*B*u + Dpzt*b*phi) * W
    T stress[3] = {0e0, 0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        T DBj[6];
        compute_DB(dNdx+2*j, DBj);
        for (int k = 0; k < 3; ++k) {
            stress[k] += DBj[k]*u[3*j] + DBj[k+3]*u[3*j+1]
                        -Dpzt[k]*N[j]*u[3*j+2];
        }
    }
    for (int i = 0; i < nen; ++i)
        add_BDBu(dNdx+2*i, stress, W, &F(0,i));

    // -- Electrical part
    // -- Form F += b'*(Dpz*B*u - kde*b*phi) * W
    T edisp_flux[2] = {0e0, 0e0};
    for (int i = 0; i < nen; ++i)
        F(2,i) = (dNdx[2*i+0]*edisp_flux[0] + dNdx[2*i+1]*edisp_flux[1])*W;
}


template<class T> inline
void ME::add_Ma_pz(double* N, T* a, T W, QMatrix<T>& F)
{
    // -- Form F += rho*kron(N*N', I2) * a * dcomplex
    T aa[2] = {0e0, 0e0};
    W *= rho;
    for (int j = 0; j < nen; ++j)
        for (int k = 0; k < 2; ++k)
            aa[k] += N[j]*a[3*j+k] * W;
    for (int i = 0; i < nen; ++i)
        for (int k = 0; k < 2; ++k)
            F(k,i) += N[i]*aa[k];
}


void ME::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx, double cv, double ca)
{
    int id[MAXNEN*3];
    nen = mesh->get_nen(eltid);
    // Are the functions set_local_id, set_local_arrays defined??
    set_local_id(mesh, eltid, id);
    Quad2d shape(mesh, eltid, nen);

    if (has_stretch()) {

       dcomplex nodestretch[MAXNEN*2];
       set_local_stretch(mesh, eltid, nodestretch, nen, 2);
       PMLShape pshape(shape, nodestretch);

       QMatrix<dcomplex> Ke(NULL, 3*nen, 3*nen);
       for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
           pshape.eval(iter.X());
           dcomplex Wt = iter.W()*pshape.Jt();
           if (cx)     add_K_pz(pshape.N(), pshape.dNdxt(), cx * Wt, Ke);
           if (ca)     add_M_pz(pshape.N(),                 ca * Wt, Ke);
       }
       K->add(id, 3*nen, Ke.data);

    } else {

        QMatrix<double> Ke(NULL, 3*nen, 3*nen);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            if (cx)     add_K_pz(shape.N(), shape.dNdx(), cx * Wt, Ke);
            if (ca)     add_M_pz(shape.N(),               ca * Wt, Ke);
        }
        K->add(id, 3*nen, Ke.data);
    }
}


void ME::assemble_R(Mesh* mesh, int eltid)
{
    nen = mesh->get_nen(eltid);
    Quad2d shape(mesh, eltid, nen);
    if (has_stretch()) {

        dcomplex nodestretch[MAXNEN*2];
        set_local_stretch(mesh, eltid, nodestretch, nen, 2);
        PMLShape pshape(shape, nodestretch);

        QMatrix<dcomplex> Fe(NULL, 3, nen);
        dcomplex nodeu[MAXNEN*3];
        dcomplex nodea[MAXNEN*3];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.Jt();
            add_Ku_pz(pshape.N(), pshape.dNdxt(), nodeu, Wt, Fe);
            add_Ma_pz(pshape.N(),                 nodea, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);

    } else {

        QMatrix<double> Fe(NULL, 3, nen);
        double nodeu[MAXNEN*3];
        double nodea[MAXNEN*3];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            add_Ku_pz(shape.N(), shape.dNdx(), nodeu, Wt, Fe);
            add_Ma_pz(shape.N(),               nodea, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);
    }
}
