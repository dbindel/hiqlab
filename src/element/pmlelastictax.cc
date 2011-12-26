/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

#include "qmatrix.h"
#include "gaussquad.h"
#include "qcomplex.h"
#include "shapes.h"
#include "pmlelastictax.h"
#include "mesh.h"

using std::vector;

#define ME           PMLElasticTAxis
#define MAXNEN       25


/*
 * Some notes on formulation:
 *
 * The algebra for higher-order axisymmetric simulations is similar to
 * the algebra for ordinary axisymmetry... but much messier.  The full
 * B matrix is given in equation 9.27 of Z&T, vol 2.  However, the B
 * matrix can be written in two parts as
 *
 *    B = B1 cos(l theta) + B2 sin(l theta)
 *
 * which is useful, since we don't really want to form D*B and B'*D*B;
 * we want to form integrals of D*B and B'*D*B with respect to theta.
 * Integrating out theta from B'*D*B, for example, gives something
 * like
 *
 *    2*pi*B1'*D*B1             for l == 0
 *    pi*(B1'*D*B1 + B2'*D*B2)  for l != 0.
 *
 * The details are not spelled out in Z&T, but they can be worked
 * out.  The key thing to remember in the code below is that, where in
 * the 2D or 3D code I might have some product or the other, here I
 * have the *integral* of the product with respect to the polar angle
 * -- and the integrand is not constant with respect to that variable.
 *
 * I drop factors of 2*pi (for l == 0) and pi (for l != 0) in the
 * interest of convenience.  It might be smart to reconsider this later.
 *
 * The other thing to remember is that there are three degrees of
 * freedom per node, even though the mesh should be 2D.  The local
 * dofs are in r, theta, z order.
 */


ME::ME(double E, double nu, double rho, int l) :
    PMLElement(3), rho(rho), ltheta(l)
{
    double a = E/(1+nu)/(1-2*nu);
    D.clear();

    D(0,0)=(1-nu)*a;  D(0,1)=   nu *a;  D(0,2)=   nu *a;
    D(1,0)=   nu *a;  D(1,1)=(1-nu)*a;  D(1,2)=   nu *a;
    D(2,0)=   nu *a;  D(2,1)=   nu *a;  D(2,2)=(1-nu)*a;

    D(3,3)=(1-2*nu)/2*a;
    D(4,4)=(1-2*nu)/2*a;
    D(5,5)=(1-2*nu)/2*a;
}


ME::ME(double* Ddata, double rho, int l) :
    PMLElement(3), rho(rho), ltheta(l)
{
    D = Ddata;
}


ME::~ME()
{
}


template<class T> inline
void ME::compute_DB(T* dNdx, double Ndr,
                    QMatrix1<T,6,3>& DB1,
                    QMatrix1<T,6,3>& DB2)
{
    for (int k = 0; k < 6; ++k) {

        DB1(k,0) = D(k,0)*dNdx[0] + D(k,2)*Ndr + D(k,3)*dNdx[1];
        DB1(k,1) = D(k,2)*(ltheta*Ndr);
        DB1(k,2) = D(k,1)*dNdx[1] + D(k,3)*dNdx[0];

        DB2(k,0) = D(k,5)*(-ltheta*Ndr);
        DB2(k,1) = D(k,4)*dNdx[1] + D(k,5)*(dNdx[0]-Ndr);
        DB2(k,2) = D(k,4)*(-ltheta*Ndr);

    }
}


template<class T> inline
void ME::compute_BDB(T* dNdx, double Ndr,
                     QMatrix1<T,6,3>& DB1,
                     QMatrix1<T,6,3>& DB2,
                     QMatrix1<T,3,3>& BDB)
{
    for (int i = 0; i < 3; ++i) {

        BDB(0,i) =  dNdx[0]*DB1(0,i) + Ndr*DB1(2,i) + dNdx[1]*DB1(3,i);
        BDB(1,i) =  ltheta*Ndr*DB1(2,i);
        BDB(2,i) =  dNdx[1]*DB1(1,i) + dNdx[0]*DB1(3,i);

        BDB(0,i) += -ltheta*Ndr*DB2(5,i);
        BDB(1,i) += dNdx[1]*DB2(4,i) + (dNdx[0]-Ndr)*DB2(5,i);
        BDB(2,i) += -ltheta*Ndr*DB2(4,i);

    }
}


template<class T> inline
void ME::add_BDBu(T* dNdx, double Ndr, T* DB1, T* DB2,
                               T Wt, T* BDB)
{
    BDB[0] += (dNdx[0]*DB1[0] + Ndr*DB1[2] + dNdx[1]*DB1[3]) * Wt;
    BDB[1] += (ltheta*Ndr*DB1[2]                           ) * Wt;
    BDB[2] += (dNdx[1]*DB1[1]              + dNdx[0]*DB1[3]) * Wt;

    BDB[0] += (-ltheta*Ndr*DB2[5]                          ) * Wt;
    BDB[1] += (dNdx[1]*DB2[4] + (dNdx[0]-Ndr)*DB2[5]       ) * Wt;
    BDB[2] += (-ltheta*Ndr*DB2[4]                          ) * Wt;
}


template<class T> inline
void ME::add_K(T* dNdx, double* N, double r, T W, QMatrix<T>& K)
{
    // -- Form K += B'*D*B * W
    QMatrix1<T,6,3> DBj1;
    QMatrix1<T,6,3> DBj2;
    QMatrix1<T,3,3> A;
    for (int j = 0; j < nen; ++j) {
        compute_DB(dNdx+2*j, N[j]/r, DBj1, DBj2);
        for (int i = 0; i < nen; ++i) {
            compute_BDB(dNdx+2*i, N[i]/r, DBj1, DBj2, A);
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    K(i*3+k, j*3+l) += A(k,l) * W;
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
void ME::add_Ku(T* dNdx, double* N, double r, T* u, T W,
                            QMatrix<T>& F)
{
    // -- Form F += B'*(D*B*u) * W
    T stress1[6] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    T stress2[6] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        QMatrix1<T,6,3> DBj1;
        QMatrix1<T,6,3> DBj2;
        compute_DB(dNdx+2*j, N[j]/r, DBj1, DBj2);
        for (int k = 0; k < 6; ++k) {

            stress1[k] += DBj1(k,0)*u[3*j+0] +
                          DBj1(k,1)*u[3*j+1] +
                          DBj1(k,2)*u[3*j+2];

            stress2[k] += DBj2(k,0)*u[3*j+0] +
                          DBj2(k,1)*u[3*j+1] +
                          DBj2(k,2)*u[3*j+2];
        }
    }
    for (int i = 0; i < nen; ++i)
        add_BDBu(dNdx+2*i, N[i]/r, stress1, stress2, W, &F(0,i));
}


template<class T> inline
void ME::add_Ma(double* N, T* a, T W, QMatrix<T>& F)
{
    // -- Form F += rho*kron(N*N', I2) * a * dcomplex
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
    nen = mesh->get_nen(eltid);
    int id[MAXNEN*3];
    set_local_id(mesh, eltid, id);
    Quad2d shape(mesh, eltid, nen);
    if (has_stretch()) {

        dcomplex nodestretch[MAXNEN*2];
        set_local_stretch(mesh, eltid, nodestretch, nen, 2);
        PMLShape pshape(shape, nodestretch);

        QMatrix<dcomplex> Ke(NULL, 3*nen, 3*nen);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.x(0)*pshape.Jt();
            if (cx)  add_K(pshape.dNdxt(), pshape.N(), pshape.x(0), cx*Wt, Ke);
            if (ca)  add_M(pshape.N(),                              ca*Wt, Ke);
        }
        K->add(id, 3*nen, Ke.data);

    } else {

        QMatrix<double> Ke(NULL, 3*nen, 3*nen);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.x(0)*shape.J();
            if (cx)  add_K(shape.dNdx(), shape.N(), shape.x(0), cx*Wt, Ke);
            if (ca)  add_M(shape.N(),                           ca*Wt, Ke);
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
            dcomplex Wt = iter.W()*pshape.x(0)*pshape.Jt();
            add_Ku(pshape.dNdxt(), pshape.N(), pshape.x(0), nodeu, Wt, Fe);
            add_Ma(pshape.N(),                              nodea, Wt, Fe);
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
            double Wt = iter.W()*shape.x(0)*shape.J();
            add_Ku(shape.dNdx(), shape.N(), shape.x(0), nodeu, Wt, Fe);
            add_Ma(shape.N(),                           nodea, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);

    }
}


void ME::project_L2(Mesh* mesh, int eltid, int nfields,
                    double* Mdiag, double* xfields,
                    FieldEval& xfunc)
{
    nen = mesh->get_nen(eltid);
    Quad2d shape(mesh, eltid, nen);
    vector<double> Xfields(nfields);
    for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
        shape.eval(iter.X());
        double Wt = iter.W()*shape.J();
        double* N = shape.N();
        for (int i = 0; i < nen; ++i) {
            int ii = mesh->ix(i,eltid);
            Mdiag[ii] += N[i]*Wt;
            xfunc(mesh, eltid, iter.X(), &(Xfields[0]));
            for (int j = 0; j < nfields; ++j)
                xfields[ii*nfields+j] += Xfields[j]*N[i]*Wt;
        }
    }
}


double* ME::mean_power(Mesh* mesh, int eltid, double* X, double* EX)
{
    nen = mesh->get_nen(eltid);
    dcomplex nodeu[MAXNEN*3];
    dcomplex nodev[MAXNEN*3];
    QMatrix<dcomplex> u(nodeu, 3,nen);
    QMatrix<dcomplex> v(nodev, 3,nen);
    set_local_u(mesh, eltid, nodeu);
    set_local_v(mesh, eltid, nodev);

    Quad2d shape(mesh, eltid, nen);
    shape.eval(X);
    dcomplex s1[6] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    dcomplex s2[6] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    dcomplex vx[3]      = {0e0, 0e0, 0e0};

    for (int j = 0; j < nen; ++j) {
        QMatrix1<double,6,3> DB1;
        QMatrix1<double,6,3> DB2;
        compute_DB(shape.dNdx()+2*j, shape.N(j)/shape.x(0), DB1, DB2);
        for (int k = 0; k < 6; ++k) {
            s1[k] += DB1(k,0)*u(0,j) + DB1(k,1)*u(1,j) + DB1(k,2)*u(2,j);
            s2[k] += DB2(k,0)*u(0,j) + DB2(k,1)*u(1,j) + DB2(k,2)*u(2,j);
        }
        for (int k = 0; k < 3; ++k)
            vx[k] += shape.N(j) * v(k,j);
    }

    /*
     * This is a little strange.  The displacements are
     *   u = (u_r, u_theta, u_z)
     * and the stress is ordered
     *   (e_rr, e_zz, e_tt, g_rz, g_zt, g_tr)
     */

    // First the part that gets multiplied by cos(l theta)*cos(l theta)

    EX[0] = -real( s1[0]*conj(vx[0])+s1[3]*conj(vx[2])+s1[5]*conj(vx[1]) )/2;
    EX[1] = -real( s1[2]*conj(vx[1])+s1[4]*conj(vx[2])+s1[5]*conj(vx[0]) )/2;
    EX[2] = -real( s1[1]*conj(vx[2])+s1[3]*conj(vx[0])+s1[4]*conj(vx[1]) )/2;

    // Then the part that gets multiplied by cos(l theta)*sin(l theta)

    EX[3] = -real( s2[0]*conj(vx[0])+s2[3]*conj(vx[2])+s2[5]*conj(vx[1]) )/2;
    EX[4] = -real( s2[2]*conj(vx[1])+s2[4]*conj(vx[2])+s2[5]*conj(vx[0]) )/2;
    EX[5] = -real( s2[1]*conj(vx[2])+s2[3]*conj(vx[0])+s2[4]*conj(vx[1]) )/2;

    return EX+6;
}


double* ME::stress(Mesh* mesh, int eltid, double* X, double* stress)
{
    nen = mesh->get_nen(eltid);
    double nodeu[MAXNEN*3];
    QMatrix<double> u(nodeu, 3,nen);
    set_local_u(mesh, eltid, nodeu);
    Quad2d shape(mesh, eltid, nen);
    shape.eval(X);

    std::fill(stress, stress+12, 0);
    double* stress1 = stress;
    double* stress2 = stress + 6;
    for (int j = 0; j < nen; ++j) {
        QMatrix1<double,6,3> DB1;
        QMatrix1<double,6,3> DB2;
        compute_DB(shape.dNdx()+2*j, shape.N(j)/shape.x(0), DB1, DB2);
        for (int k = 0; k < 6; ++k) {
            stress1[k] += DB1(k,0)*u(0,j) + DB1(k,1)*u(1,j) + DB1(k,2)*u(2,j);
            stress2[k] += DB1(k,0)*u(0,j) + DB1(k,1)*u(1,j) + DB1(k,2)*u(2,j);
        }
    }

    return stress+12;
}

