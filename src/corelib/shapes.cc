/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "qmatrix.h"
#include "gaussquad.h"
#include "shapes.h"
#include "mesh.h"

using std::copy;
using std::swap;

#define MAXNEN       125


/* ==== Shape base class ==== */


Shape::Shape(int nen, int ndm)
{
    owner = 1;
    nen_ = nen;
    ndm_ = ndm;
    N_    = new double[nen];
    dNdx_ = new double[nen*ndm];
}


Shape::Shape()
{
    owner = 0;
}


Shape::~Shape()
{
    if (owner) {
        delete[] N_;
        delete[] dNdx_;
    }
}


void Shape::interp(double* f, double* fx, int nfields)
{
    for (int j = 0; j < nfields; ++j) {
        fx[j] = 0e0;
        for (int i = 0; i < nen_; ++i)
            fx[j] += N_[i]*f[nfields*i+j];
    }
}


void Shape::interp(dcomplex* f, dcomplex* fx, int nfields)
{
    for (int j = 0; j < nfields; ++j) {
        fx[j] = 0e0;
        for (int i = 0; i < nen_; ++i)
            fx[j] += N_[i]*f[nfields*i+j];
    }
}


double Shape::interp(double* f)
{
    double fx = 0;
    for (int i = 0; i < nen_; ++i)
        fx += N_[i]*f[i];
    return fx;
}


dcomplex Shape::interp(dcomplex* f)
{
    dcomplex fx = 0e0;
    for (int i = 0; i < nen_; ++i)
        fx += N_[i]*f[i];
    return fx;
}


ShapeInt::~ShapeInt()
{
}


/* ==== PML shape wrapper ==== */


PMLShape::PMLShape(Shape& base, dcomplex* stretch) :
    base_(base)
{
    owner = 0;
    nen_  = base.nen();
    ndm_  = base.ndm();
    N_    = base.N();
    dNdx_ = base.dNdx();

    if (stretch) {
        stretch_.resize(nen()*ndm());
        dNdxt_.resize(nen()*ndm());
        copy(stretch, stretch+(nen_*ndm_), stretch_.begin());
//        copy(stretch_.begin(), stretch_.end(), stretch); // usage of copy not correct
    }
}


void PMLShape::eval(double* X)
{
    dcomplex stretchX[3];
    base_.eval(X);
    base_.interp(&stretch_[0], stretchX, base_.ndm());

    for (int j = 0; j < nen_; ++j)
        for (int i = 0; i < ndm_; ++i)
            dNdxt(i,j) = base_.dNdx(i,j) / stretchX[i];

    Jt_ = base_.J();
    for (int i = 0; i < ndm_; ++i)        Jt_ *= stretchX[i];
    for (int i = 0; i < ndm_; ++i)        x_[i] = base_.x(i);
}


/* ==== 1D case ==== */


Shape1d::Shape1d(int nen) :
    Shape(nen,1), nodex(nen)
{
}


/* Compute the Jacobian matrix and determinant for the
 * isoparametric map.
 */
void Shape1d::compute_jacobian(double* dNdX, double* J, double& Jdet)
{
    *J = 0;
    for (int k = 0; k < nen_; ++k)
        *J += nx(k) * dNdX[k];
    Jdet = *J;
}


/* Compute spatial derivatives of shape functions
 */
void Shape1d::remap_gradients(double* J, double* dNdX)
{
    for (int i = 0; i < nen_; ++i)
        dNdx(0,i) = dNdX[i] / (*J);
}


/* ==== 2D case ==== */


Shape2d::Shape2d(int nen) :
    Shape(nen,2), nodex(2*nen)
{
}


/* Compute the Jacobian matrix and determinant for the
 * isoparametric map.
 */
void Shape2d::compute_jacobian(double* dNdX1, double* J1, double& Jdet)
{
    QMatrix<double> J    (J1,    2, 2  );
    QMatrix<double> dNdX (dNdX1, 2, nen_);

    J.clear();
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < nen_; ++k)
                J(i,j) += nx(i,k) * dNdX(j,k);

    Jdet = J(0,0)*J(1,1) - J(0,1)*J(1,0);
}


/* Compute spatial derivatives of shape functions
 */
void Shape2d::remap_gradients(double* J1, double* dNdX1)
{
    QMatrix<double> J   (J1,    2, 2  );
    QMatrix<double> dNdX(dNdX1, 2, nen_);

    double Jdet = J(0,0)*J(1,1) - J(1,0)*J(0,1);

    for (int i = 0; i < nen_; ++i) {
        dNdx(0,i) = ( dNdX(0,i)*J(1,1) - dNdX(1,i)*J(1,0) )/Jdet;
        dNdx(1,i) = (-dNdX(0,i)*J(0,1) + dNdX(1,i)*J(0,0) )/Jdet;
    }
}


/* ==== 3D case ==== */


Shape3d::Shape3d(int nen) :
    Shape(nen,3), nodex(3*nen)
{
}


/* Compute the Jacobian matrix and determinant for the
 * isoparametric map.
 */
void Shape3d::compute_jacobian(double* dNdX1, double* J1, double& Jdet)
{
    QMatrix<double> J    (J1,    3, 3   );
    QMatrix<double> dNdX (dNdX1, 3, nen_);

    J.clear();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < nen_; ++k)
                J(i,j) += nx(i,k) * dNdX(j,k);

    Jdet =
        J(0,0) * (J(1,1)*J(2,2)-J(2,1)*J(1,2)) -
        J(1,0) * (J(0,1)*J(2,2)-J(2,1)*J(0,2)) +
        J(2,0) * (J(0,1)*J(1,2)-J(1,1)*J(0,2));
}


/* Compute spatial derivatives of shape functions
 */
void Shape3d::remap_gradients(double* J1, double* dNdX1)
{
    double LU1[9];
    QMatrix<double> J   (J1,    3, 3   );
    QMatrix<double> dNdX(dNdX1, 3, nen_);
    QMatrix<double> LU  (LU1,   3, 3   );

    LU   = J;
    copy(dNdX1, dNdX1+3*nen_, dNdx_);
//    copy(dNdx_, dNdx_+3*nen_, dNdX1);

    for (int i = 0; i < 2; ++i) {

        // Choose pivot
        int piv = i;
        for (int k = i+1; k < 3; ++k)
            if (fabs(LU(i,k)) > fabs(LU(i,piv)))
                piv = k;

        // Apply pivot
        if (piv != i) {
            for (int j = 0; j < 3;    ++j) swap( LU(j,i),   LU(j,piv)   );
            for (int j = 0; j < nen_; ++j) swap( dNdx(i,j), dNdx(piv,j) );
        }

        // Eliminate
        for (int j = i+1; j < 3; ++j)
            LU(i,j) /= LU(i,i);
        for (int j = i+1; j < 3; ++j)
            for (int k = i+1; k < 3; ++k)
                LU(k,j) -= LU(i,j)*LU(k,i);
    }

    for (int j = 0; j < nen_; ++j) {

        // Forward substitute with L
        dNdx(1,j) -=  dNdx(0,j)*LU(0,1);
        dNdx(2,j) -= (dNdx(0,j)*LU(0,2) + dNdx(1,j)*LU(1,2));

        // Back substitute with U
        dNdx(2,j) = (                            dNdx(2,j)        ) / LU(2,2);
        dNdx(1,j) = (          dNdx(1,j)        -dNdx(2,j)*LU(2,1)) / LU(1,1);
        dNdx(0,j) = (dNdx(0,j)-dNdx(1,j)*LU(1,0)-dNdx(2,j)*LU(2,0)) / LU(0,0);
    }
}


/* ==== Lagrangian shape functions ==== */


/* Compute shape functions and their gradients in the parent
 * coordinate system.
 */
void shape1d(int order, double* N, double* dN, double xx)
{
    if (order == 1) {

        N[0] = (1-xx)/2;
        N[1] = (1+xx)/2;

        dN[0] = -1./2;
        dN[1] =  1./2;

    } else if (order == 2) {

        N[0] =    -xx*(1-xx)/2;
        N[1] = (1-xx)*(1+xx);
        N[2] =     xx*(1+xx)/2;

        dN[0] = (2*xx-1)/2;
        dN[1] = -2*xx;
        dN[2] = (2*xx+1)/2;

    } else if (order == 3) {

        N[0]  =         (1+3*xx)*(1-3*xx)*(1-xx) * -1./16;
        N[1]  =  (1+xx)*         (1-3*xx)*(1-xx) *  9./16;
        N[2]  =  (1+xx)*(1+3*xx)*         (1-xx) *  9./16;
        N[3]  =  (1+xx)*(1+3*xx)*(1-3*xx)        * -1./16;

        dN[0] = (-1 - 18*xx + 27*xx*xx) * -1./16;
        dN[1] = (-3 -  2*xx +  9*xx*xx) *  9./16;
        dN[2] = ( 3 -  2*xx -  9*xx*xx) *  9./16;
        dN[3] = ( 1 - 18*xx - 27*xx*xx) * -1./16;

    } else if (order == 4) {

        N[0] =        (1+2*xx)*xx*(1-2*xx)*(1-xx) /  6.   ;
        N[1] = (1+xx)*         xx*(1-2*xx)*(1-xx) * -4./3 ;
        N[2] = (1+xx)*(1+2*xx)*   (1-2*xx)*(1-xx);
        N[3] = (1+xx)*(1+2*xx)*xx         *(1-xx) *  4./3 ;
        N[4] = (1+xx)*(1+2*xx)*xx*(1-2*xx)        / -6.   ;

        dN[0] = ( -8*xx*xx*(1-xx) + (1-2*xx)*(1+2*xx)*(1-2*xx)      ) / 6. ;
        dN[1] = ( -2*xx*xx*(1-2*xx) + (1-xx)*(1+xx)*(1-4*xx) ) * -4./3;
        dN[2] = xx*( -8*(1-xx)*(1+xx) - 2*(1-2*xx)*(1+2*xx)  );
        dN[3] = ( -2*xx*xx*(1+2*xx) + (1-xx)*(1+xx)*(1+4*xx) ) *  4./3;
        dN[4] = ( -8*xx*xx*(1+xx) + (1+2*xx)*(1+2*xx)*(1-2*xx)      ) / -6. ;

    }
}


/* Compute shape functions and their gradients in the parent
 * coordinate system.
 */
void shape2d(int nen, double* N, double* dN1, double* x)
{
    QMatrix<double> dN(dN1, 2,nen);

    double Nx[5], dNx[5];
    double Ny[5], dNy[5];

    int order = order2d(nen);
    shape1d(order, Nx, dNx, x[0]);
    shape1d(order, Ny, dNy, x[1]);

    int i = 0;
    for (int ix = 0; ix < order+1; ++ix) {
        for (int iy = 0; iy < order+1; ++iy) {
            N [i]   =  Nx[ix]* Ny[iy];
            dN(0,i) = dNx[ix]* Ny[iy];
            dN(1,i) =  Nx[ix]*dNy[iy];
            ++i;
        }
    }
}


/* Compute shape functions and their gradients in the parent
 * coordinate system.
 */
void shape3d(int nen, double* N, double* dN1, double* x)
{
    QMatrix<double> dN(dN1, 3,nen);

    double Nx[3], dNx[3];
    double Ny[3], dNy[3];
    double Nz[3], dNz[3];

    int order = order3d(nen);
    shape1d(order, Nx, dNx, x[0]);
    shape1d(order, Ny, dNy, x[1]);
    shape1d(order, Nz, dNz, x[2]);

    int i = 0;
    for (int ix = 0; ix < order+1; ++ix) {
        for (int iy = 0; iy < order+1; ++iy) {
            for (int iz = 0; iz < order+1; ++iz) {
                N [i]   =  Nx[ix]* Ny[iy]* Nz[iz];
                dN(0,i) = dNx[ix]* Ny[iy]* Nz[iz];
                dN(1,i) =  Nx[ix]*dNy[iy]* Nz[iz];
                dN(2,i) =  Nx[ix]* Ny[iy]*dNz[iz];
                ++i;
            }
        }
    }
}


/* ==== 1D Lagrangian case ==== */


Quad1d::Quad1d(Mesh* mesh, int eltid, int nen) :
    Shape1d(nen)
{
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        nx(j) = mesh->x(0,nodeid);
    }
}


void Quad1d::eval(double* X)
{
    double dNdX[MAXNEN];
    double Jmat;

    shape1d(nen_, N_, dNdX, *X);
//    shape1d(nen_-1, N_, dNdX, *X);
    compute_jacobian(dNdX, &Jmat, J_);
    remap_gradients(&Jmat, dNdX);
    interp(&nodex[0], &x_[0], 1);
}


void Quad1d::inv(double* x, double* X, int maxiter, double tol)
{
    double dNdX[MAXNEN];
    double Jmat;
    double resid;
    double dX;
    int    iter = 0;

    shape1d(nen_, N_, dNdX, *X);
    compute_jacobian(dNdX, &Jmat, J_);
    interp(&nodex[0], &x_[0], 1);
    resid = x_[0] - x[0];

    while (fabs(resid) > tol && iter < maxiter) {
        iter++;
        apply_inv(&Jmat, &resid, &dX);
        *X = *X - dX;
        shape1d(nen_, N_, dNdX, *X);
        compute_jacobian(dNdX, &Jmat, J_);
        interp(&nodex[0], &x_[0], 1);
        resid = x_[0] - x[0];
    }
}


void Quad1d::apply_inv(double* Jmat, double* x, double* xi)
{
    xi[0] = x[0] / (*Jmat);
}


/* ==== 2D Lagrangian case ==== */


Quad2d::Quad2d(Mesh* mesh, int eltid, int nen) :
    Shape2d(nen)
{
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        for (int i = 0; i < 2; ++i)
            nx(i,j) = mesh->x(i,nodeid);
    }
}


void Quad2d::eval(double* X)
{
    double dNdX[MAXNEN*2];
    double Jmat[4];

    shape2d(nen_, N_, dNdX, X);
    compute_jacobian(dNdX, Jmat, J_);
    remap_gradients(Jmat, dNdX);
    interp(&nodex[0], &x_[0], 2);
}


void Quad2d::inv(double* x, double* X, int maxiter, double tol)
{
    double dNdX[MAXNEN*2];
    double Jmat[4];
    double resid[2];
    double nresid;
    double dX[2];
    int    iter = 0;

    shape2d(nen_, N_, dNdX, X);
    compute_jacobian(dNdX, Jmat, J_);
    interp(&nodex[0], &x_[0], 2);
    nresid = 0;
    for (int i = 0; i < 2; ++i)  {
        resid[i] = x_[i] - x[i];
        nresid   += resid[i]*resid[i];
    }
    nresid = sqrt(nresid);

    while ( nresid > tol && iter < maxiter) {
        iter++;
        apply_inv(Jmat, resid, dX);
        for (int i = 0; i < 2; ++i) 
            X[i] -= dX[i];
        shape2d(nen_, N_, dNdX, X);
        compute_jacobian(dNdX, Jmat, J_);
        interp(&nodex[0], &x_[0], 2);
        nresid = 0;
        for (int i = 0; i < 2; ++i)  {
            resid[i] = x_[i] - x[i];
            nresid   += resid[i]*resid[i];
        }
        nresid = sqrt(nresid);
    }
}


void Quad2d::apply_inv(double* J1, double* x, double* xi)
{
    QMatrix<double> J   (J1,    2, 2  );

    double Jdet = J(0,0)*J(1,1) - J(1,0)*J(0,1);

    xi[0] = ( x[0]*J(1,1) - x[1]*J(0,1) )/Jdet;
    xi[1] = (-x[0]*J(1,0) + x[1]*J(0,0) )/Jdet;
}


/* ==== 3D case ==== */


Quad3d::Quad3d(Mesh* mesh, int eltid, int nen) :
    Shape3d(nen)
{
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        for (int i = 0; i < 3; ++i)
            nx(i,j) = mesh->x(i,nodeid);
    }
}


void Quad3d::eval(double* X)
{
    double dNdX[MAXNEN*3];
    double Jmat[9];

    shape3d(nen_, N_, dNdX, X);
    compute_jacobian(dNdX, Jmat, J_);
    remap_gradients(Jmat, dNdX);
    interp(&nodex[0], &x_[0], 3);
}


void Quad3d::inv(double* x, double* X, int maxiter, double tol)
{
    double dNdX[MAXNEN*3];
    double Jmat[9];
    double resid[3];
    double nresid;
    double dX[3];
    int    iter = 0;

    shape3d(nen_, N_, dNdX, X);
    compute_jacobian(dNdX, Jmat, J_);
    interp(&nodex[0], &x_[0], 3);
    nresid = 0;
    for (int i = 0; i < 3; ++i)  {
        resid[i] = x_[i] - x[i];
        nresid   += resid[i]*resid[i];
    }
    nresid = sqrt(nresid);

    while ( nresid > tol && iter < maxiter) {
        iter++;
        apply_inv(Jmat, resid, dX);
        for (int i = 0; i < 3; ++i) 
            X[i]-=dX[i];
        shape3d(nen_, N_, dNdX, X);
        compute_jacobian(dNdX, Jmat, J_);
        interp(&nodex[0], &x_[0], 3);
        nresid = 0;
        for (int i = 0; i < 3; ++i)  {
            resid[i] = x_[i] - x[i];
            nresid   += resid[i]*resid[i];
        }
        nresid = sqrt(nresid);
    }
}


void Quad3d::apply_inv(double* J1, double* x, double* xi)
{
    double LU1[9];
    QMatrix<double> J   (J1,    3, 3   );
    QMatrix<double> LU  (LU1,   3, 3   );

    LU   = J;
    copy(x, x+3, xi);
    for (int i = 0; i < 2; ++i) {

        // Choose pivot
        int piv = i;
        for (int k = i+1; k < 3; ++k)
            if (fabs(LU(k,i)) > fabs(LU(piv,i)))
                piv = k;

        // Apply pivot
        if (piv != i) {
            for (int j = 0; j < 3;    ++j) swap( LU(i,j),   LU(piv,j)   );
                                           swap(    xi[i],  xi[piv]     );
        }

        // Eliminate
        for (int j = i+1; j < 3; ++j)
            LU(j,i) /= LU(i,i);
        for (int j = i+1; j < 3; ++j)
            for (int k = i+1; k < 3; ++k)
                LU(k,j) -= LU(i,j)*LU(k,i);
    }


    // Forward substitute with L
    xi[1] -=  xi[0]*LU(1,0);
    xi[2] -= (xi[0]*LU(2,0) + xi[1]*LU(2,1));

    // Back substitute with U
    xi[2] = (                        xi[2]        ) / LU(2,2);
    xi[1] = (          xi[1]        -xi[2]*LU(1,2)) / LU(1,1);
    xi[0] = (    xi[0]-xi[1]*LU(0,1)-xi[2]*LU(0,2)) / LU(0,0);
}


/* ==== 2D integration iterator ==== */


Quad1dInt::Quad1dInt(int nen)
{
    ix = 0;
    ngauss = nen;
    eval();
}


int Quad1dInt::end()
{
    return (ix >= ngauss);
}


int Quad1dInt::next()
{
    ++ix;
    if (ix >= ngauss) return 0;
    eval();
    return 1;
}


void Quad1dInt::eval()
{
    W_ = gauss_weight(ix, ngauss);
    X_[0] = gauss_point(ix, ngauss);
}


/* ==== 2D integration iterator ==== */


Quad2dInt::Quad2dInt(int nen)
{
    ix = 0;
    iy = 0;
    ngauss = order2d(nen)+1;
    eval();
}


int Quad2dInt::end()
{
    return (ix >= ngauss);
}


int Quad2dInt::next()
{
    ++iy;
    if (iy >= ngauss) { iy = 0; ++ix; }
    if (ix >= ngauss) return 0;
    eval();
    return 1;
}


void Quad2dInt::eval()
{
    W_ = gauss_weight(ix, ngauss)*gauss_weight(iy, ngauss);
    X_[0] = gauss_point(ix, ngauss);
    X_[1] = gauss_point(iy, ngauss);
}


/* ==== 3D integration iterator ==== */


Quad3dInt::Quad3dInt(int nen)
{
    ix = 0;
    iy = 0;
    iz = 0;
    ngauss = order3d(nen)+1;
    eval();
}


int Quad3dInt::end()
{
    return (ix >= ngauss);
}


int Quad3dInt::next()
{
    ++iz;
    if (iz >= ngauss) { iz = 0; ++iy; }
    if (iy >= ngauss) { iy = 0; ++ix; }
    if (ix >= ngauss) return 0;
    eval();
    return 1;
}


void Quad3dInt::eval()
{
    W_ = gauss_weight(ix, ngauss)*
        gauss_weight(iy, ngauss)*
        gauss_weight(iz, ngauss);
    X_[0] = gauss_point(ix, ngauss);
    X_[1] = gauss_point(iy, ngauss);
    X_[2] = gauss_point(iz, ngauss);
}

