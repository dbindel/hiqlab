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
#include "pmlelastic2d_pz.h"
#include "material_model.h"
#include "mesh.h"

#define ME          PMLElastic2d_pz
#define MAXNEN      25


inline void Bmatrix() 
{
    /* <generator matexpr>
      function Bmatrix(dN) = [dN(1), 0;
                                  0, dN(2);
                              dN(2), dN(1)];
    */
}


ME::ME(double E, double nu, double rho,
       double* pz, double kds, int plane_type):
  PMLElement(3), rho(rho)
{
    double lambda = E*nu/(1+nu)/(1-2*nu);
    double mu     = E/2/(1+nu);
    double Db[25];

    // Plane strain / stress isotopric
    if (plane_type == 0)
      piezo_elasticity_2D_strain(Db, lambda, mu, pz, kds);
    else
      piezo_elasticity_2D_stress(Db, lambda, mu, pz, kds);
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
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++)
            D[j+3*i] = Db[j+5*i];
        for(int j = 3; j < 5; j++)
            Dpzt[(j-3)+2*i] = -Db[j+5*i];
    }
    for(int i = 3; i < 5; i++) {
        for(int j = 3; j < 5; j++)
            kde[(j-3)+2*(i-3)] = Db[j+5*i];
    }
}


template<class T> inline
void ME::add_K(double* N, T* dNdx, T W, QMatrix<T>& K)
{
    // Mechanical part: Kmm += B'*D*B  * W
    // Electrical part: Kee -= b'*kde*b  * W
    // Coupling term:   Kmt += B'*Dpzt*b * W

    int nen3 = 3*nen;
    for (int j = 0; j < nen; ++j) {
        for (int i = 0; i < nen; ++i) {
            T* dNj = dNdx + 2*j;
            T* dNi = dNdx + 2*i;
            T* Kij = &(K(3*i,3*j));
            /* <generator matexpr complex="T">
              input D symmetric(3), kde symmetric(2), Dpzt(3,2);
              complex input dNi(2), dNj(2), W;
              complex inout Kij[nen3](3,3);
              Bi = Bmatrix(dNi);
              Bj = Bmatrix(dNj);
              Kij += [ Bi'*(D*Bj),     Bi'*(Dpzt*dNj);
                      (dNi'*Dpzt')*Bj, dNi'*kde *dNj] * W;
            */
        }
    }
}


template<class T> inline
void ME::add_M(double* N, T W, QMatrix<T>& K)
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
void ME::add_Ku(double* N, T* dNdx, T* u, T W, QMatrix<T>& F)
{
    T stress[3] = {0e0, 0e0, 0e0};
    T e_disp[2] = {0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        T* uj = u+3*j;
        T* dNj = dNdx+2*j;
        /* <generator matexpr complex="T">
          input D symmetric(3), kde symmetric(2), Dpzt(3,2);
          complex input uj(3), dNj(2);
          complex inout stress(3) += [    D*Bmatrix(dNj), Dpzt*dNj]*uj; 
          complex inout e_disp(2) += [Dpzt'*Bmatrix(dNj), -kde*dNj]*uj;
        */
    }
    for (int i = 0; i < nen; ++i) {
        T* dNi = dNdx+2*i;
        T* Fi  = &(F(0,i));
        /* <generator matexpr complex="T">
          complex input dNi(2), stress(3), e_disp(2), W;
          complex inout Fi(3) += [Bmatrix(dNi)'*stress; 
                                  dNi'*e_disp]*W;
        */
    }
}


template<class T> inline
void ME::add_Ma(double* N, T* a, T W, QMatrix<T>& F)
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
           if (cx)     add_K(pshape.N(), pshape.dNdxt(), cx * Wt, Ke);
           if (ca)     add_M(pshape.N(),                 ca * Wt, Ke);
       }
       K->add(id, 3*nen, Ke.data);

    } else {

        QMatrix<double> Ke(NULL, 3*nen, 3*nen);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            if (cx)     add_K(shape.N(), shape.dNdx(), cx * Wt, Ke);
            if (ca)     add_M(shape.N(),               ca * Wt, Ke);
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
            add_Ku(pshape.N(), pshape.dNdxt(), nodeu, Wt, Fe);
            add_Ma(pshape.N(),                 nodea, Wt, Fe);
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
            add_Ku(shape.N(), shape.dNdx(), nodeu, Wt, Fe);
            add_Ma(shape.N(),               nodea, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);
    }
}
