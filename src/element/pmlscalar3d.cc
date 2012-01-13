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
#include "pmlscalar3d.h"
#include "mesh.h"

#define ME           PMLScalar3d
#define MAXNEN       125


ME::ME(double kappa, double rho) :
    PMLElement(1), rho(rho)
{
    D[0] = kappa;  D[3] = 0;      D[6] = 0;
    D[1] = 0;      D[4] = kappa;  D[7] = 0;
    D[2] = 0;      D[5] = 0;      D[8] = kappa;
}


ME::ME(double* D, double rho) :
    PMLElement(3), rho(rho)
{
    std::copy(D, D+9, this->D);
}


ME::~ME()
{
}


template<class T> inline
T ME::dotD(T* a, T* b)
{
    return (a[0]*D[0]*b[0] + a[1]*D[1]*b[0] + a[2]*D[2]*b[0] +
            a[0]*D[1]*b[1] + a[1]*D[4]*b[1] + a[2]*D[5]*b[1] +
            a[0]*D[2]*b[2] + a[1]*D[5]*b[2] + a[2]*D[8]*b[2]);
}


template<class T> inline
void ME::add_K(T* dNdx, T W, QMatrix<T>& K)
{
    // -- Form K += dN*D*dN' * W
    for (int j = 0; j < nen; ++j)
        for (int i = 0; i < nen; ++i)
            K(i,j) += dotD(dNdx+3*i, dNdx+3*j) * W;
}


template<class T> inline
void ME::add_M(double* N, T W, QMatrix<T>& K)
{
    // -- Form M += rho*N*N' * W
    for (int j = 0; j < nen; ++j)
        for (int i = 0; i < nen; ++i)
            K(i,j) += rho*N[i]*N[j] * W;
}


template<class T> inline
void ME::add_Ku(T* dNdx, T* u, T W, QMatrix<T>& F)
{
    // -- Form K += dN*D*dN'*u * W
    T dudx[3] = {0e0, 0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        dudx[0] += dNdx[3*j+0] * u[j];
        dudx[1] += dNdx[3*j+1] * u[j];
    }
    for (int i = 0; i < nen; ++i)
        F(i) += dotD(dNdx+3*i, dudx) * W;
}


template<class T> inline
void ME::add_Ma(double* N, T* a, T W, QMatrix<T>& F)
{
    // -- Form F += rho*N*N' * a
    T aa = 0e0;
    W *= rho;
    for (int j = 0; j < nen; ++j)
        aa += N[j]*a[j];
    for (int i = 0; i < nen; ++i)
        F(i) += N[i]*aa * W;
}


void ME::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx, double cv, double ca)
{
    nen = mesh->get_nen(eltid);
    int id[MAXNEN];
    set_local_id(mesh, eltid, id);
    Quad3d shape(mesh, eltid, nen);

    if (has_stretch()) {

        dcomplex nodestretch[MAXNEN*3];
        set_local_stretch(mesh, eltid, nodestretch, nen, 3);
        PMLShape pshape(shape, nodestretch);

        QMatrix<dcomplex> Ke(NULL, nen, nen);
        for (Quad3dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.Jt();
            if (cx)        add_K(pshape.dNdxt(), cx * Wt, Ke);
            if (ca)        add_M(pshape.N(),     ca * Wt, Ke);
        }
        K->add(id, nen, Ke.data);

    } else {

        QMatrix<double> Ke(NULL, nen, nen);
        for (Quad3dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            if (cx)        add_K(shape.dNdx(), cx * Wt, Ke);
            if (ca)        add_M(shape.N(),    ca * Wt, Ke);
        }
        K->add(id, nen, Ke.data);

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

        QMatrix<dcomplex> Fe(NULL, nen, 1);
        dcomplex nodeu[MAXNEN];
        dcomplex nodea[MAXNEN];
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

        QMatrix<double> Fe(NULL, 1, nen);
        double nodeu[MAXNEN];
        double nodea[MAXNEN];
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
    double u[MAXNEN];
    set_local_u(mesh, eltid, u);
    Quad3d shape(mesh, eltid, nen);
    shape.eval(X);

    double dudx[3] = {0e0, 0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        dudx[0] += shape.dNdx(0,j) * u[j];
        dudx[1] += shape.dNdx(1,j) * u[j];
        dudx[2] += shape.dNdx(2,j) * u[j];
    }
    stress[0] = D[0]*dudx[0] + D[1]*dudx[1] + D[2]*dudx[2];
    stress[1] = D[1]*dudx[0] + D[4]*dudx[1] + D[5]*dudx[2];
    stress[2] = D[2]*dudx[0] + D[5]*dudx[1] + D[8]*dudx[2];

    return stress+3;
}
