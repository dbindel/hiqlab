/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "qmatrix.h"
#include "gaussquad.h"
#include "qcomplex.h"
#include "shapes.h"
#include "pmlelastic3d.h"
#include "mesh.h"
#include "material_model.h"

#define ME           PMLElastic3d
#define MAXNEN       125


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


ME::ME(double E, double nu, double rho) :
    PMLElement(3), rho(rho)
{
    double lambda = E*nu/(1+nu)/(1-2*nu);
    double mu     = E/2/(1+nu);

    elasticity_3D(D.data,lambda,mu);
}


ME::ME(double* Ddata, double rho) :
    PMLElement(3), rho(rho)
{
    D = Ddata;
}


ME::~ME()
{
}


template<class Ts, class T> inline
void ME::compute_stress(Ts* dNdx, T* u, T* sig)
{
    T eps[6] = {0e0, 0e0, 0e0, 0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        Ts* dNj = dNdx+3*j;
        T*  uj  = u+3*j;
        /* <generator matexpr complex="T">
          complex input dNj(3), uj(3);
          complex inout eps(6) += Bmatrix(dNj)*uj;
        */
    }
    /* <generator matexpr complex="T">
      input D symmetric(6);
      complex input eps(6); 
      output sig(6) = D*eps;
    */
}


template<class T> inline
void ME::add_K(T* dNdx, T W, QMatrix<T>& K)
{
    // -- Form K += B'*D*B * W
    int nen3 = 3*nen;
    for (int j = 0; j < nen; ++j) {
        for (int i = 0; i < nen; ++i) {
            T* dNj = dNdx + 3*j;
            T* dNi = dNdx + 3*i;
            T* Kij = &(K(3*i,3*j));
            /* <generator matexpr complex="T">
              input D symmetric(6);
              complex input dNi(3), dNj(3), W;
              complex inout Kij[nen3](3,3);
              Kij += (Bmatrix(dNi)'*D*Bmatrix(dNj)) * W;
            */
        }
    }
}


template<class T> inline
void ME::add_M(double* N, T W, QMatrix<T>& K)
{
    // -- Form M += rho*kron(N*N', I3) * W
    for (int j = 0; j < nen; ++j) {
        for (int i = 0; i < nen; ++i) {
            T Mij = rho*N[i]*N[j]*W;
            for (int k = 0; k < 3; ++k)
                K(i*3+k, j*3+k) += Mij;
        }
    }
}


template<class T> inline
void ME::add_Ku(T* dNdx, T* u, T W, QMatrix<T>& F)
{
    // -- Form F += B'*(D*B*u) * W
    T sig[6];
    compute_stress(dNdx, u, sig);
    for (int i = 0; i < nen; ++i) {
        T* dNi = dNdx+3*i;
        T* Fi  = &(F(0,i));
        /* <generator matexpr complex="T">
          complex input dNi(3), sig(6), W;
          complex inout Fi(3) += (Bmatrix(dNi)'*sig) * W;
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
            aa[k] += N[j]*a[3*j+k] * W;
    for (int i = 0; i < nen; ++i)
        for (int k = 0; k < 3; ++k)
            F(k,i) += N[i]*aa[k];
}


void ME::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx, double cv, double ca)
{
    int id[MAXNEN*3];
    nen = mesh->get_nen(eltid);
    set_local_id(mesh, eltid, id);
    Quad3d shape(mesh, eltid, nen);

    if (has_stretch()) {

        dcomplex nodestretch[MAXNEN*3];
        set_local_stretch(mesh, eltid, nodestretch, nen, 3);
        PMLShape pshape(shape, nodestretch);

        QMatrix<dcomplex> Ke(NULL, 3*nen, 3*nen);
        for (Quad3dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.Jt();
            if (cx)        add_K(pshape.dNdxt(), cx * Wt, Ke);
            if (ca)        add_M(pshape.N(),     ca * Wt, Ke);
        }
        K->add(id, 3*nen, Ke.data);

    } else {

        QMatrix<double> Ke(NULL, 3*nen, 3*nen);
        for (Quad3dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            if (cx)        add_K(shape.dNdx(), cx * Wt, Ke);
            if (ca)        add_M(shape.N(),    ca * Wt, Ke);
        }
        K->add(id, 3*nen, Ke.data);

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

        QMatrix<dcomplex> Fe(NULL, 3, nen);
        dcomplex nodeu[MAXNEN*3];
        dcomplex nodea[MAXNEN*3];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        for (Quad3dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.Jt();
            add_Ku(pshape.dNdxt(), nodeu, Wt, Fe);
            add_Ma(pshape.N(),     nodea, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);

    } else {

        QMatrix<double> Fe(NULL, 3, nen);
        double nodeu[MAXNEN*3];
        double nodea[MAXNEN*3];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        for (Quad3dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            add_Ku(shape.dNdx(), nodeu, Wt, Fe);
            add_Ma(shape.N(),    nodea, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);

    }
}


double* ME::stress(Mesh* mesh, int eltid, double* X, double* stress)
{
    nen = mesh->get_nen(eltid);
    double nodeu[MAXNEN*3];
    set_local_u(mesh, eltid, nodeu);
    Quad2d shape(mesh, eltid, nen);
    shape.eval(X);
    compute_stress(shape.dNdx(), nodeu, stress);
    return stress+6;
}

