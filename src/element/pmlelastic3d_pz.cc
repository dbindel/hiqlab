/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstring>
#include <cstdlib>
#include <cmath>

#include "qmatrix.h"
#include "gaussquad.h"
#include "qcomplex.h"
#include "shapes.h"
#include "pmlelastic3d_pz.h"
#include "mesh.h"
#include "material_model.h"

#define ME          PMLElastic3d_pz
#define MAXNEN      125


inline void Bmatrix()
{
    /* <generator matexpr>
      function Bmatrix(dN) = 
          [dN(1),     0,     0;
               0, dN(2),     0;
               0,     0, dN(3);
           dN(2), dN(1),     0;
               0, dN(3), dN(2);
           dN(3),     0, dN(1)];
    */
}


ME::ME(double E, double nu, double rho, double* pz, double kds):
    PMLElement(4), rho(rho)
{
    double lambda = E*nu/(1+nu)/(1-2*nu);
    double mu     = E/2/(1+nu);
    double Db[81];

    // -- 3D Isotropic
    piezo_elasticity_3D(Db,lambda,mu,pz,kds);
    split_Db(Db);
}


ME::ME(double* Db, double rho):
    PMLElement(4), rho(rho)
{
    split_Db(Db);
}


ME::~ME()
{
}


void ME::split_Db(double* Db)
{
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++)
            D[j+6*i] = Db[j+9*i];
        for (int j = 6; j < 9; j++)
            Dpzt[(j-6)+3*i] = -Db[j+9*i];
    }
    for (int i = 6; i < 9; i++) {
        for (int j = 6; j < 9; j++)
            kde[(j-6)+3*(i-6)] = Db[j+9*i];
    }
}


template<class T> inline
void ME::add_K(double* N, T* dNdx, T W, QMatrix<T>& K)
{
    // Mechanical part: Kmm += B'*D*B  * W
    // Electrical part: Kee -= b'*kde*b  * W
    // Coupling term:   Kmt += B'*Dpzt*b * W

    int nen4 = 4*nen;
    for (int j = 0; j < nen; ++j) {
        for (int i = 0; i < nen; ++i) {
            T* dNj = dNdx + 3*j;
            T* dNi = dNdx + 3*i;
            T* Kij = &(K(4*i,4*j));
            /* <generator matexpr complex="T">
              input D symmetric(6), kde symmetric(3), Dpzt(6,3);
              complex input dNi(3), dNj(3), W;
              complex inout Kij[nen4](4,4);
              Bi = Bmatrix(dNi);
              Bj = Bmatrix(dNj);
              Kij += [ Bi'*(D*Bj),     Bi'*(Dpzt*dNj);
                      (dNi'*Dpzt')*Bj, dNi'*kde *dNj] * W;
            */
        }
    }
}


template <class T> inline
void ME::add_M(double* N, T W, QMatrix<T>& K)
{
    // -- Form M += rho*kron(N*N', I3) * W
    for (int j = 0; j < nen; ++j) {
        for (int i = 0; i < nen; ++i) {
            T Mij = rho*N[i]*N[j]*W;
            for (int k = 0; k < 3; ++k)
                K(i*4+k, j*4+k) += Mij;
        }
    }
}


template<class T> inline
void ME::add_Ku(double* N, T* dNdx, T* u, T W, QMatrix<T>& F)
{
    T stress[9] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    T e_disp[3] = {0e0, 0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        T* uj = u+4*j;
        T* dNj = dNdx+3*j;
        /* <generator matexpr complex="T">
          input D symmetric(6), kde symmetric(3), Dpzt(6,3);
          complex input uj(4), dNj(3);
          complex inout stress(6) += [    D*Bmatrix(dNj), Dpzt*dNj]*uj; 
          complex inout e_disp(3) += [Dpzt'*Bmatrix(dNj), -kde*dNj]*uj;
        */
    }
    for (int i = 0; i < nen; ++i) {
        T* dNi = dNdx+3*i;
        T* Fi  = &(F(0,i));
        /* <generator matexpr complex="T">
          complex input dNi(3), stress(6), e_disp(3), W;
          complex inout Fi(4) += [Bmatrix(dNi)'*stress; 
                                  dNi'*e_disp]*W;
        */
    }
}


template<class T> inline
void ME::add_Ma(double* N, T* a, T W, QMatrix<T>& F)
{
    // -- Form F += rho*kron(N*N', I3) * a * dcomplex
    T aa[3] = {0e0, 0e0, 0e0};
    W *= rho;
    for (int j = 0; j < nen; ++j)
        for (int k = 0; k < 3; ++k)
            aa[k] += N[j]*a[4*j+k] * W;
    for (int i = 0; i < nen; ++i)
        for (int k = 0; k < 3; ++k)
            F(k,i) += N[i]*aa[k];
}


void ME::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                               double cx, double cv, double ca)
{
    int id[MAXNEN*4];
    nen = mesh->get_nen(eltid);
    set_local_id(mesh, eltid, id);
    Quad3d shape(mesh, eltid, nen);

    if (has_stretch()) {

        dcomplex nodestretch[MAXNEN*3];
        set_local_stretch(mesh, eltid, nodestretch, nen, 3);
        PMLShape pshape(shape, nodestretch);

        QMatrix<dcomplex> Ke(NULL, 4*nen, 4*nen);
        for (Quad3dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.Jt();
            if (cx)  add_K(pshape.N(), pshape.dNdxt(), cx * Wt, Ke);
            if (ca)  add_M(pshape.N(),                 ca * Wt, Ke);
        }
        K->add(id, 4*nen, Ke.data);

    } else {

        QMatrix<double> Ke(NULL, 4*nen, 4*nen);
        for (Quad3dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            if (cx)  add_K(shape.N(), shape.dNdx(), cx * Wt, Ke);
            if (ca)  add_M(shape.N(),               ca * Wt, Ke);
        }
        K->add(id, 4*nen, Ke.data);

    }
}


void ME::assemble_R(Mesh* mesh, int eltid)
{
    nen = mesh->get_nen(eltid);
    Quad3d shape(mesh, eltid, nen);
    if (has_stretch()) {

        dcomplex nodestretch[MAXNEN*3];
        set_local_stretch(mesh, eltid, nodestretch, nen, 3);
        PMLShape pshape(shape, nodestretch);

        QMatrix<dcomplex> Fe(NULL, 4, nen);
        dcomplex nodeu[MAXNEN*4];
        dcomplex nodea[MAXNEN*4];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        for (Quad3dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.Jt();
            add_Ku(pshape.N(), pshape.dNdxt(), nodeu, Wt, Fe);
            add_Ma(pshape.N(),                 nodea, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);

    } else {

        QMatrix<double> Fe(NULL, 4, nen);
        double nodeu[MAXNEN*4];
        double nodea[MAXNEN*4];
        double nodev[MAXNEN*4];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        set_local_v(mesh, eltid, nodev);
        for (Quad3dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            add_Ku(shape.N(), shape.dNdx(), nodeu, Wt, Fe);
            add_Ma(shape.N(),               nodea, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);

    }
}

