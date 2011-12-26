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
#include "pmlscalar1d.h"
#include "mesh.h"

#define ME           PMLScalar1d
#define MAXNEN       5


ME::ME(double kappa, double rho) :
    PMLElement(1), rho(rho), nen(nen)
{
    D = kappa;
}


ME::~ME()
{
}


template<class T> inline
void ME::add_K(T* dNdx, T W, QMatrix<T>& K)
{
    // -- Form K += dN*D*dN' * W
    for (int j = 0; j < nen; ++j)
        for (int i = 0; i < nen; ++i)
            K(i,j) += dNdx[i]*D*dNdx[j] * W;
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
    T dudx = 0e0;
    for (int j = 0; j < nen; ++j)
        dudx += dNdx[j]*u[j];
    for (int i = 0; i < nen; ++i)
        F(i) += dNdx[i]*D*dudx * W;
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
    Quad1d shape(mesh, eltid, nen);

    if (has_stretch()) {

        dcomplex nodestretch[MAXNEN];
        set_local_stretch(mesh, eltid, nodestretch, nen, 1);
        PMLShape pshape(shape, nodestretch);

        QMatrix<dcomplex> Ke(NULL, nen, nen);
        for (Quad1dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.Jt();
            if (cx)        add_K(pshape.dNdxt(), cx * Wt, Ke);
            if (ca)        add_M(pshape.N(),     ca * Wt, Ke);
        }
        K->add(id, nen, Ke.data);

    } else {

        QMatrix<double> Ke(NULL, nen, nen);
        for (Quad1dInt iter(nen); !iter.end(); iter.next()) {
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
    Quad1d shape(mesh, eltid, nen);
    if (has_stretch()) {

        dcomplex nodestretch[MAXNEN];
        set_local_stretch(mesh, eltid, nodestretch, nen, 1);
        PMLShape pshape(shape, nodestretch);

        QMatrix<dcomplex> Fe(NULL, nen, 1);
        dcomplex nodeu[MAXNEN];
        dcomplex nodea[MAXNEN];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        for (Quad1dInt iter(nen); !iter.end(); iter.next()) {
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
        for (Quad1dInt iter(nen); !iter.end(); iter.next()) {
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
    Quad1d shape(mesh, eltid, nen);
    shape.eval(X);

    double dudx = 0e0;
    for (int j = 0; j < nen; ++j)
        dudx += shape.dNdx(0,j) * u[j];
    stress[0] = D*dudx;

    return stress+1;
}

