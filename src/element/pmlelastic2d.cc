/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "luasupport.h"
#include "qmatrix.h"
#include "gaussquad.h"
#include "qcomplex.h"
#include "shapes.h"
#include "pmlelastic2d.h"
#include "mesh.h"
#include "material_model.h"

using std::vector;

#define ME           PMLElastic2d
#define MAXNEN       25


inline void Bmatrix() 
{
    /* <generator matexpr>
      function Bmatrix(dN) = [dN(1), 0;
                                  0, dN(2);
                              dN(2), dN(1)];
    */
}


ME::ME(double E, double nu, double rho, int plane_type) :
    PMLElement(2), L(NULL), rho(rho)
{
    double lambda = E*nu/(1+nu)/(1-2*nu);
    double mu     = E/2/(1+nu);

    if (plane_type == 0)   // -- Plane strain
        elasticity_2D_strain(D.data, lambda, mu);
    else                   // -- Plane stress
        elasticity_2D_stress(D.data, lambda, mu);
}


ME::ME(double* Ddata, double rho) :
    PMLElement(2), L(NULL), rho(rho)
{
    D = Ddata;
}


ME::ME(lua_State* L, int func) :
    PMLElement(2)
{
    if (!lua_isfunction(L,func)) {
        this->L = NULL;
        D.clear();
        rho = 0;
    } else {
        this->L = L;
        lua_pushstring(L, "material");
        lua_pushvalue(L,func);
        lua_objsettable(L, this);
    }
}


ME::~ME()
{
}


void ME::compute_D(Mesh* mesh, int eltid, double* N)
{
    if (!L)
        return;

    QMatrix<double> localx(nen,2);
    set_local_x(mesh, eltid, localx.data, 2);

    double px[2] = {0e0, 0e0};
    for (int i = 0; i < nen; ++i) {
        px[0] += N[i]*localx(i,0);
        px[1] += N[i]*localx(i,1);
    }

    lua_pushstring(L, "material");
    lua_objgettable(L, this);

    for (int i = 0; i < 2; ++i)
        lua_pushnumber(L, px[i]);
    lua_pushnumber(L, eltid);

    if (lua_pcall2(L, 3, 4) == 0) {
        if (lua_istable(L,-4)) { // Returns D, rho

            for (int i = 0; i < 9; ++i) {
                lua_rawgeti(L,-4,i+1);
                D.data[i] = lua_tonumber(L,-1);
                lua_pop(L,1);
            }
            rho = lua_tonumber(L,-3);

        } else { // Returns E, nu, rho, plane_type

            double E          = lua_tonumber(L,-4);
            double nu         = lua_tonumber(L,-3);
            rho               = lua_tonumber(L,-2);
            int plane_type    = (int) lua_tonumber(L,-1);

            double lambda = E*nu/(1+nu)/(1-2*nu);
            double mu     = E/2/(1+nu);
            if (plane_type == 0) {
                elasticity_2D_strain(D.data,lambda,mu);
            } else {
                elasticity_2D_stress(D.data,lambda,mu);
            }
        }
        lua_pop(L, 4);
    }
}


template<class Ts, class T> inline
void ME::compute_stress(Ts* dNdx, T* u, T* sig)
{
    T eps[3] = {0e0, 0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        Ts* dNj = dNdx+2*j;
        T*  uj  = u+2*j;
        /* <generator matexpr complex="T">
          complex input dNj(2), uj(2);
          complex inout eps(3) += Bmatrix(dNj)*uj;
        */
    }
    /* <generator matexpr complex="T">
      input D symmetric(3);
      complex input eps(3); 
      output sig(3) = D*eps;
    */
}


template<class T> inline
void ME::add_K(T* dNdx, T W, QMatrix<T>& K)
{
    // -- Form K += B'*D*B * W
    int nen2 = 2*nen;
    for (int j = 0; j < nen; ++j) {
        for (int i = 0; i < nen; ++i) {
            T* dNj = dNdx + 2*j;
            T* dNi = dNdx + 2*i;
            T* Kij = &(K(2*i,2*j));
            /* <generator matexpr complex="T">
              input D symmetric(3);
              complex input dNi(2), dNj(2), W;
              complex inout Kij[nen2](2,2);
              Kij += (Bmatrix(dNi)'*D*Bmatrix(dNj)) * W;
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
                K(i*2+k, j*2+k) += Mij;
        }
    }
}


template<class T> inline
void ME::add_Ku(T* dNdx, T* u, T W, QMatrix<T>& F)
{
    // -- Form F += B'*(D*B*u) * W
    T sig[3];
    compute_stress(dNdx, u, sig);
    for (int i = 0; i < nen; ++i) {
        T* dNi = dNdx+2*i;
        T* Fi  = &(F(0,i));
        /* <generator matexpr complex="T">
          complex input dNi(2), sig(3), W;
          complex inout Fi(2) += (Bmatrix(dNi)'*sig) * W;
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

    if (has_stretch()) {

        dcomplex nodestretch[MAXNEN*2];
        set_local_stretch(mesh, eltid, nodestretch, nen, 2);
        PMLShape pshape(shape, nodestretch);

        QMatrix<dcomplex> Ke(NULL, 2*nen, 2*nen);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.Jt();
            compute_D(mesh, eltid, pshape.N());
            if (cx)        add_K(pshape.dNdxt(), cx * Wt, Ke);
            if (ca)        add_M(pshape.N(),     ca * Wt, Ke);
        }
        K->add(id, 2*nen, Ke.data);

    } else {

        QMatrix<double> Ke(NULL, 2*nen, 2*nen);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            compute_D(mesh, eltid, shape.N());
            if (cx)        add_K(shape.dNdx(), cx * Wt, Ke);
            if (ca)        add_M(shape.N(),    ca * Wt, Ke);
        }
        K->add(id, 2*nen, Ke.data);

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

        QMatrix<dcomplex> Fe(NULL, 2, nen);
        dcomplex nodeu[MAXNEN*2];
        dcomplex nodea[MAXNEN*2];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            pshape.eval(iter.X());
            dcomplex Wt = iter.W()*pshape.Jt();
            compute_D(mesh, eltid, pshape.N());
            add_Ku(pshape.dNdxt(), nodeu, Wt, Fe);
            add_Ma(pshape.N(),     nodea, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);

    } else {

        QMatrix<double> Fe(NULL, 2, nen);
        double nodeu[MAXNEN*2];
        double nodea[MAXNEN*2];
        set_local_u(mesh, eltid, nodeu);
        set_local_a(mesh, eltid, nodea);
        for (Quad2dInt iter(nen); !iter.end(); iter.next()) {
            shape.eval(iter.X());
            double Wt = iter.W()*shape.J();
            compute_D(mesh, eltid, shape.N());
            add_Ku(shape.dNdx(), nodeu, Wt, Fe);
            add_Ma(shape.N(),    nodea, Wt, Fe);
        }
        add_local_f(mesh, eltid, Fe.data);

    }
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
    set_local_u(mesh, eltid, nodeu);
    set_local_v(mesh, eltid, nodev);
    Quad2d shape(mesh, eltid, nen);
    shape.eval(X);

    dcomplex vx[2] = {0e0, 0e0};
    for (int j = 0; j < nen; ++j) {
        double Nj = shape.N(j);
        vx[0] += nodev[2*j+0] * Nj;
        vx[1] += nodev[2*j+1] * Nj;
    }

    dcomplex stress[3];
    compute_D(mesh, eltid, shape.N());
    compute_stress(shape.dNdx(), nodeu, stress);
    /* <generator matexpr complex="dcomplex">
      complex input stress(3);
      complex input vx(2);
      sigt = [stress(1), stress(3); stress(3), stress(2)];
      output EX(2) = -real( sigt*conj(vx) )/2;
    */
    return EX+2;
}


double* ME::stress(Mesh* mesh, int eltid, double* X, double* stress)
{
    nen = mesh->get_nen(eltid);
    double nodeu[MAXNEN*2];
    set_local_u(mesh, eltid, nodeu);
    Quad2d shape(mesh, eltid, nen);
    shape.eval(X);
    compute_D(mesh, eltid, shape.N());
    compute_stress(shape.dNdx(), nodeu, stress);
    return stress+3;
}

