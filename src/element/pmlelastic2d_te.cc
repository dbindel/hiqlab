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
#include "pmlelastic2d_te.h"
#include "material_model.h"
#include "mesh.h"

#define ME          PMLElastic2d_te
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
       double at, double cp, double kt, double T0,
       int plane_type):
    PMLElement(3), rho(rho), at(at), cp(cp), kt(kt), T0(T0)
{

    double lambda = E*nu/(1+nu)/(1-2*nu);
    double mu     = E/2/(1+nu);
    double Db[16];

    // Isotropic plane strain / plane stress
    if (plane_type == 0)
        thermoelasticity_2D_strain(Db,lambda,mu,at);
    else
        thermoelasticity_2D_stress(Db,lambda,mu,at);
    split_Db(Db);
}


ME::ME(double* Db, double rho, double at, double cp, double kt, double T0):
    PMLElement(3), rho(rho), at(at), cp(cp), kt(kt), T0(T0)
{
    split_Db(Db);
}


void ME::split_Db(double* Db)
{
    /* <generator matexpr>
      input Db(4,4);
      Pm = [eye(3); 0, 0, 0];     // Just mechanical components
      Pt = [0; 0; 0; 1];          // Just thermal component
      output D(3,3) = Pm'*Db*Pm;  // Purely mechanical
      output vp(3)  = Pm'*Db*Pt;  // Coupling
      output vCooiv = Pt'*Db*Pt;  // Purely thermal
    */
}


ME::~ME()
{
}


template<class T> inline
void ME::add_K(double* N, T* dNdx, T W, QMatrix<T>& K)
{
    // Mechanical part of the stiffness: Kmm += B'*D*B * W
    // Thermal part of the stiffness:    Ktt += b'*kt(I2)*b * W
    // Coupling term of the stiffness:   Kmt -= alphat * B' * Vp* N * W

    int nen3 = 3*nen;
    for (int j = 0; j < nen; ++j) {
        for (int i = 0; i < nen; ++i) {
            double Nj = N[j];
            T* dNj = dNdx + 2*j;
            T* dNi = dNdx + 2*i;
            T* Kij = &(K(3*i,3*j));
            /* <generator matexpr complex="T">
              input D symmetric(3), vp(3), kt, Nj;
              complex input dNi(2), dNj(2), W;
              complex inout Kij[nen3](3,3);
              Bi = Bmatrix(dNi);
              Bj = Bmatrix(dNj);
              Kij += [Bi'*D*Bj, -Bi'*vp*Nj;
                      0,0,       kt*(dNi'*dNj)] * W;
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
void ME::add_C(double* N, T* dNdx, T W, QMatrix<T>& K)
{
    // Thermal term:   Ctt += const * N*N' * W
    // Coupling term:  Ctm += alphat * theta0 * N' * Vp'* B * W

    int nen3 = 3*nen;
    double Ctherm = (rho*cp+T0*vCooiv);
    for (int j = 0; j < nen; ++j) {
        for (int i = 0; i < nen; ++i) {
            double Nj = N[j];
            double Ni = N[i];
            T* dNj = dNdx + 2*j;
            T* Ktx_ij = &(K(3*i+2,3*j));
            /* <generator matexpr complex="T">
              input D symmetric(3), vp(3), Ni, Nj, Ctherm, T0;
              complex input dNj(2), W;
              complex inout Ktx_ij[nen3](1,3);
              Ktx_ij += [T0*(Ni*vp'*Bmatrix(dNj)), Ctherm*Ni*Nj] * W;
            */
        }
    }
}


template<class T> inline
void ME::add_Ku(double* N, T* dNdx, T* u, T W, QMatrix<T>& F)
{
    T stress[3] = {0e0, 0e0, 0e0};
    T heat_flux[2] = {0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        T* utj     = u+3*j;
        T* dNj     = dNdx+2*j;
        double Nj  = N[j];
        /* <generator matexpr complex="T">
          input D symmetric(3), vp(3), kt, Nj;
          complex input utj(3), dNj(2);
          complex inout stress(3)    += [D*Bmatrix(dNj), -vp*Nj]*utj;
          complex inout heat_flux(2) += kt*dNj*utj(3);
        */
    }
    for (int i = 0; i < nen; ++i) {
        T* dNi = dNdx+2*i;
        T* Fi  = &(F(0,i));
        /* <generator matexpr complex="T">
          complex input dNi(2), stress(3), heat_flux(2), W;
          complex inout Fi(3) += [Bmatrix(dNi)'*stress; 
                                  dNi'*heat_flux] * W;
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


template<class T> inline
void ME::add_Cv(double* N, T* dNdx, T* v, T W, QMatrix<T>& F)
{
    // -- Form F += N'*( theta0*VpB*v(1:2)+ const*v(3) )
    T thetad = 0e0;
    double Ctherm = rho*cp+T0*vCooiv;
    for (int j = 0; j < nen; ++j) {
        double Nj = N[j];
        T* dNj = dNdx+2*j;
        T* utj = v+3*j;
        /* <generator matexpr complex="T">
          input vp(3), T0, Ctherm, Nj;
          complex input utj(3), dNj(2);
          complex inout thetad += [T0*(vp'*Bmatrix(dNj)), Ctherm*Nj]*utj;
        */
    }
    for (int i = 0; i < nen; ++i)
        F(2,i) += N[i] * thetad * W;
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
           if (cv)     add_C(pshape.N(), pshape.dNdxt(), cv * Wt, Ke);
       }
       K->add(id, 3*nen, Ke.data);

    } else {

        QMatrix<double> Ke(NULL, 3*nen, 3*nen);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            if (cx)     add_K(shape.N(), shape.dNdx(), cx * Wt, Ke);
            if (ca)     add_M(shape.N(),               ca * Wt, Ke);
            if (cv)     add_C(shape.N(), shape.dNdx(), cv * Wt, Ke);
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
        dcomplex nodev[MAXNEN*3];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        set_local_v(mesh, eltid, nodev);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.Jt();
            add_Ku(pshape.N(), pshape.dNdxt(), nodeu, Wt, Fe);
            add_Ma(pshape.N(),                 nodea, Wt, Fe);
            add_Cv(pshape.N(), pshape.dNdxt(), nodev, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);

    } else {

        QMatrix<double> Fe(NULL, 3, nen);
        double nodeu[MAXNEN*3];
        double nodea[MAXNEN*3];
        double nodev[MAXNEN*3];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        set_local_v(mesh, eltid, nodev);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            add_Ku(shape.N(), shape.dNdx(), nodeu, Wt, Fe);
            add_Ma(shape.N(),               nodea, Wt, Fe);
            add_Cv(shape.N(), shape.dNdx(), nodev, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);
    }
}
