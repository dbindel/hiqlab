/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

#include "qmatrix.h"
#include "gaussquad.h"
#include "qcomplex.h"
#include "shapes.h"
#include "elastic2d.h"
#include "mesh.h"
#include "material_model.h"

using std::vector;

#define ME           Elastic2d
#define MAXNEN       25


ME::ME(double E, double nu, double rho, int plane_type) :
    Element(2), rho(rho)
{
    double lambda = E*nu/(1+nu)/(1-2*nu);
    double mu     = E/2/(1+nu);

    if (plane_type == 0) {

        // -- Plane strain
        elasticity_2D_strain(D.data,lambda,mu);

    } else {

        // -- Plane stress
        elasticity_2D_stress(D.data,lambda,mu);

    }
}


ME::ME(double* Ddata, double rho) :
    Element(2), rho(rho)
{
    D = Ddata;
}


ME::~ME()
{
}


void ME::compute_DB(double* dNdx, QMatrix1<double,3,2>& DB)
{
    for (int k = 0; k < 3; ++k) {
        DB(k,0) = D(k,0)*dNdx[0] + D(k,2)*dNdx[1];
        DB(k,1) = D(k,1)*dNdx[1] + D(k,2)*dNdx[0];
    }
}


void ME::compute_BDB(double* dNdx, QMatrix1<double,3,2>& DB,
                     QMatrix1<double,2,2>& BDB)
{
    for (int i = 0; i < 2; ++i) {
        BDB(0,i) = dNdx[0]*DB(0,i) + dNdx[1]*DB(2,i);
        BDB(1,i) = dNdx[1]*DB(1,i) + dNdx[0]*DB(2,i);
    }
}


void ME::add_BDBu(double* dNdx, double* DB, double Wt, double* BDB)
{
    BDB[0] += (dNdx[0]*DB[0] + dNdx[1]*DB[2])*Wt;
    BDB[1] += (dNdx[1]*DB[1] + dNdx[0]*DB[2])*Wt;
}


void ME::add_K(double* dNdx, double W, QMatrix<double>& K)
{
    // -- Form K += B'*D*B * W
    QMatrix1<double,3,2> DBj;
    QMatrix1<double,2,2> A;
    for (int j = 0; j < nen; ++j) {
        compute_DB(dNdx+2*j, DBj);
        for (int i = 0; i < nen; ++i) {
            compute_BDB(dNdx+2*i, DBj, A);
            for (int k = 0; k < 2; ++k)
                for (int l = 0; l < 2; ++l)
                    K(i*2+k, j*2+l) += A(k,l) * W;
        }
    }
}


void ME::add_M(double* N, double W, QMatrix<double>& K)
{
    // -- Form M += rho*kron(N*N', I2) * W
    for (int j = 0; j < nen; ++j) {
        for (int i = 0; i < nen; ++i) {
            double Mij = rho*N[i]*N[j]*W;
            for (int k = 0; k < 2; ++k)
                K(i*2+k, j*2+k) += Mij;
        }
    }
}


void ME::add_Ku(double* dNdx, double* u, double W, QMatrix<double>& F)
{
    // -- Form F += B'*(D*B*u) * W
    double stress[3] = {0e0, 0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        QMatrix1<double,3,2> DBj;
        compute_DB(dNdx+2*j, DBj);
        for (int k = 0; k < 3; ++k)
            stress[k] += DBj(k,0)*u[2*j] + DBj(k,1)*u[2*j+1];
    }
    for (int i = 0; i < nen; ++i)
        add_BDBu(dNdx+2*i, stress, W, &F(0,i));
}


void ME::add_Ma(double* N, double* a, double W, QMatrix<double>& F)
{
    // -- Form F += rho*kron(N*N', I2) * a * dcomplex
    double aa[2] = {0e0, 0e0};
    W *= rho;
    for (int j = 0; j < nen; ++j)
        for (int k = 0; k < 2; ++k)
            aa[k] += N[j]*a[2*j+k] * W;
    for (int i = 0; i < nen; ++i)
        for (int k = 0; k < 2; ++k)
            F(k,i) += N[i]*aa[k];
}


void ME::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx, double cv, double ca)
{
    int id[MAXNEN*2];
    nen = mesh->get_nen(eltid);
    set_local_id(mesh, eltid, id);
    Quad2d shape(mesh, eltid, nen);

    QMatrix<double> Ke(NULL, 2*nen, 2*nen);
    for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
        shape.eval(iter.X());
        double Wt = iter.W()*shape.J();
        if (cx)        add_K(shape.dNdx(), cx * Wt, Ke);
        if (ca)        add_M(shape.N(),    ca * Wt, Ke);
    }
    K->add(id, 2*nen, Ke.data);
}


void ME::assemble_R(Mesh* mesh, int eltid)
{
    nen = mesh->get_nen(eltid);
    Quad2d shape(mesh, eltid, nen);

    QMatrix<double> Fe(NULL, 2, nen);
    double nodeu[MAXNEN*2];
    double nodea[MAXNEN*2];
    set_local_u(mesh, eltid, nodeu);
    set_local_a(mesh, eltid, nodea);
    for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
        shape.eval(iter.X());
        double Wt = iter.W()*shape.J();
        add_Ku(shape.dNdx(), nodeu, Wt, Fe);
        add_Ma(shape.N(),    nodea, Wt, Fe);
    }
    add_local_f(mesh, eltid, Fe.data);
}


void ME::project_L2(Mesh* mesh, int eltid, int nfields,
                    double* Mdiag, double* xfields,
                    FieldEval& xfunc)
{
    nen = mesh->get_nen(eltid);
    vector<double> Xfields(nfields);
    Quad2d shape(mesh, eltid, nen);
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
    dcomplex nodeu[MAXNEN*2];
    dcomplex nodev[MAXNEN*2];
    QMatrix<dcomplex> u(nodeu, 2,nen);
    QMatrix<dcomplex> v(nodev, 2,nen);
    set_local_u(mesh, eltid, nodeu);
    set_local_v(mesh, eltid, nodev);

    Quad2d shape(mesh, eltid, nen);
    shape.eval(X);
    dcomplex stress[3] = {0e0, 0e0, 0e0};
    dcomplex vx[2]     = {0e0, 0e0};

    for (int j = 0; j < nen; ++j) {
        QMatrix1<double,3,2> DB;
        compute_DB(shape.dNdx()+2*j, DB);
        for (int k = 0; k < 3; ++k)
            stress[k] += DB(k,0)*u(0,j) + DB(k,1)*u(1,j);
        for (int k = 0; k < 2; ++k)
            vx[k] += shape.N(j) * v(k,j);
    }

    EX[0] = -real( stress[0]*conj(vx[0]) + stress[2]*conj(vx[1]) )/2;
    EX[1] = -real( stress[1]*conj(vx[1]) + stress[2]*conj(vx[0]) )/2;
    return EX+2;
}


double* ME::stress(Mesh* mesh, int eltid, double* X, double* stress)
{
    nen = mesh->get_nen(eltid);
    double nodeu[MAXNEN*2];
    QMatrix<double> u(nodeu, 2,nen);
    set_local_u(mesh, eltid, nodeu);
    Quad2d shape(mesh, eltid, nen);
    shape.eval(X);

    std::fill(stress, stress+3, 0);
    for (int j = 0; j < nen; ++j) {
        QMatrix1<double,3,2> DB;
        compute_DB(shape.dNdx()+2*j, DB);
        for (int k = 0; k < 3; ++k)
            stress[k] += DB(k,0)*u(0,j) + DB(k,1)*u(1,j);
    }
    return stress+3;
}

