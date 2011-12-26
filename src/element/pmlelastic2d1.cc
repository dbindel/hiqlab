#line 1 "pmlelastic2d.cc"
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
        complex input dNj(2);
        complex input uj(2);
        complex inout eps(3);

        Bj = [dNj(1), 0;
                   0, dNj(2);
              dNj(2), dNj(1)];

        eps += Bj*uj;
        */
        /* <generated matexpr> */ {
#line 125 "pmlelastic2d.cc"
        T tmp2_ = dNj[0*2+0];
        T tmp3_ = dNj[0*2+1];
#line 126 "pmlelastic2d.cc"
        T tmp5_ = uj[0*2+0];
        T tmp6_ = uj[0*2+1];
#line 127 "pmlelastic2d.cc"
        T tmp8_ = eps[0*3+0];
        T tmp9_ = eps[0*3+1];
        T tmp10_ = eps[0*3+2];
#line 131 "pmlelastic2d.cc"
#line 133 "pmlelastic2d.cc"
        T tmp13_ = tmp2_ * tmp5_;
        T tmp14_ = tmp3_ * tmp6_;
        T tmp15_ = tmp3_ * tmp5_;
        T tmp16_ = tmp2_ * tmp6_;
        T tmp17_ = tmp15_ + tmp16_;
        T tmp18_ = tmp8_ + tmp13_;
        T tmp19_ = tmp9_ + tmp14_;
        T tmp20_ = tmp10_ + tmp17_;
#line 127 "pmlelastic2d.cc"
        eps[0*3+0] = tmp18_;
        eps[0*3+1] = tmp19_;
        eps[0*3+2] = tmp20_;
        } /* </generated matexpr> */
#line 135 "pmlelastic2d.cc"
    }

    /* <generator matexpr complex="T">
    input D symmetric(3);
    complex input eps(3); 
    output sig(3);
    sig = D*eps;
    */
    /* <generated matexpr> */ {
#line 138 "pmlelastic2d.cc"
    double tmp2_ = D[0*3+0];
    double tmp3_ = D[1*3+0];
    double tmp4_ = D[1*3+1];
    double tmp5_ = D[2*3+0];
    double tmp6_ = D[2*3+1];
    double tmp7_ = D[2*3+2];
#line 139 "pmlelastic2d.cc"
    T tmp9_ = eps[0*3+0];
    T tmp10_ = eps[0*3+1];
    T tmp11_ = eps[0*3+2];
#line 141 "pmlelastic2d.cc"
    T tmp13_ = tmp2_ * tmp9_;
    T tmp14_ = tmp3_ * tmp10_;
    T tmp15_ = tmp13_ + tmp14_;
    T tmp16_ = tmp5_ * tmp11_;
    T tmp17_ = tmp15_ + tmp16_;
    T tmp18_ = tmp3_ * tmp9_;
    T tmp19_ = tmp4_ * tmp10_;
    T tmp20_ = tmp18_ + tmp19_;
    T tmp21_ = tmp6_ * tmp11_;
    T tmp22_ = tmp20_ + tmp21_;
    T tmp23_ = tmp5_ * tmp9_;
    T tmp24_ = tmp6_ * tmp10_;
    T tmp25_ = tmp23_ + tmp24_;
    T tmp26_ = tmp7_ * tmp11_;
    T tmp27_ = tmp25_ + tmp26_;
#line 140 "pmlelastic2d.cc"
    sig[0*3+0] = tmp17_;
    sig[0*3+1] = tmp22_;
    sig[0*3+2] = tmp27_;
    } /* </generated matexpr> */
#line 143 "pmlelastic2d.cc"
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
            complex input W;
            complex input dNi(2);
            complex input dNj(2);
            complex inout Kij[nen2](2,2);

            Bi = [dNi(1), 0;
                       0, dNi(2);
                  dNi(2), dNi(1)];

            Bj = [dNj(1), 0;
                       0, dNj(2);
                  dNj(2), dNj(1)];

            Kij += (Bi'*D*Bj) * W;
            */
            /* <generated matexpr> */ {
#line 158 "pmlelastic2d.cc"
            double tmp2_ = D[0*3+0];
            double tmp3_ = D[1*3+0];
            double tmp4_ = D[1*3+1];
            double tmp5_ = D[2*3+0];
            double tmp6_ = D[2*3+1];
            double tmp7_ = D[2*3+2];
#line 159 "pmlelastic2d.cc"
            T tmp9_ = W;
#line 160 "pmlelastic2d.cc"
            T tmp11_ = dNi[0*2+0];
            T tmp12_ = dNi[0*2+1];
#line 161 "pmlelastic2d.cc"
            T tmp14_ = dNj[0*2+0];
            T tmp15_ = dNj[0*2+1];
#line 162 "pmlelastic2d.cc"
            T tmp17_ = Kij[0*nen2+0];
            T tmp18_ = Kij[0*nen2+1];
            T tmp19_ = Kij[1*nen2+0];
            T tmp20_ = Kij[1*nen2+1];
#line 166 "pmlelastic2d.cc"
#line 170 "pmlelastic2d.cc"
#line 172 "pmlelastic2d.cc"
            T tmp24_ = tmp11_ * tmp2_;
            T tmp25_ = tmp12_ * tmp5_;
            T tmp26_ = tmp24_ + tmp25_;
            T tmp27_ = tmp12_ * tmp3_;
            T tmp28_ = tmp11_ * tmp5_;
            T tmp29_ = tmp27_ + tmp28_;
            T tmp30_ = tmp11_ * tmp3_;
            T tmp31_ = tmp12_ * tmp6_;
            T tmp32_ = tmp30_ + tmp31_;
            T tmp33_ = tmp12_ * tmp4_;
            T tmp34_ = tmp11_ * tmp6_;
            T tmp35_ = tmp33_ + tmp34_;
            T tmp36_ = tmp12_ * tmp7_;
            T tmp37_ = tmp28_ + tmp36_;
            T tmp38_ = tmp11_ * tmp7_;
            T tmp39_ = tmp31_ + tmp38_;
            T tmp40_ = tmp26_ * tmp14_;
            T tmp41_ = tmp37_ * tmp15_;
            T tmp42_ = tmp40_ + tmp41_;
            T tmp43_ = tmp29_ * tmp14_;
            T tmp44_ = tmp39_ * tmp15_;
            T tmp45_ = tmp43_ + tmp44_;
            T tmp46_ = tmp32_ * tmp15_;
            T tmp47_ = tmp37_ * tmp14_;
            T tmp48_ = tmp46_ + tmp47_;
            T tmp49_ = tmp35_ * tmp15_;
            T tmp50_ = tmp39_ * tmp14_;
            T tmp51_ = tmp49_ + tmp50_;
            T tmp52_ = tmp42_ * tmp9_;
            T tmp53_ = tmp45_ * tmp9_;
            T tmp54_ = tmp48_ * tmp9_;
            T tmp55_ = tmp51_ * tmp9_;
            T tmp56_ = tmp17_ + tmp52_;
            T tmp57_ = tmp18_ + tmp53_;
            T tmp58_ = tmp19_ + tmp54_;
            T tmp59_ = tmp20_ + tmp55_;
#line 162 "pmlelastic2d.cc"
            Kij[0*nen2+0] = tmp56_;
            Kij[0*nen2+1] = tmp57_;
            Kij[1*nen2+0] = tmp58_;
            Kij[1*nen2+1] = tmp59_;
            } /* </generated matexpr> */
#line 174 "pmlelastic2d.cc"
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
        complex inout Fi(2);

        Bi = [dNi(1), 0;
                   0, dNi(2);
              dNi(2), dNi(1)];

        Fi += (Bi'*sig) * W;
        */
        /* <generated matexpr> */ {
#line 204 "pmlelastic2d.cc"
        T tmp2_ = dNi[0*2+0];
        T tmp3_ = dNi[0*2+1];
#line 204 "pmlelastic2d.cc"
        T tmp5_ = sig[0*3+0];
        T tmp6_ = sig[0*3+1];
        T tmp7_ = sig[0*3+2];
#line 204 "pmlelastic2d.cc"
        T tmp9_ = W;
#line 205 "pmlelastic2d.cc"
        T tmp11_ = Fi[0*2+0];
        T tmp12_ = Fi[0*2+1];
#line 209 "pmlelastic2d.cc"
#line 211 "pmlelastic2d.cc"
        T tmp15_ = tmp2_ * tmp5_;
        T tmp16_ = tmp3_ * tmp7_;
        T tmp17_ = tmp15_ + tmp16_;
        T tmp18_ = tmp3_ * tmp6_;
        T tmp19_ = tmp2_ * tmp7_;
        T tmp20_ = tmp18_ + tmp19_;
        T tmp21_ = tmp17_ * tmp9_;
        T tmp22_ = tmp20_ * tmp9_;
        T tmp23_ = tmp11_ + tmp21_;
        T tmp24_ = tmp12_ + tmp22_;
#line 205 "pmlelastic2d.cc"
        Fi[0*2+0] = tmp23_;
        Fi[0*2+1] = tmp24_;
        } /* </generated matexpr> */
#line 213 "pmlelastic2d.cc"
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
    output EX(2);

    sigt = [stress(1), stress(3); stress(3), stress(1)];
    EX = -real( sigt*conj(vx) )/2;
    */
    /* <generated matexpr> */ {
#line 360 "pmlelastic2d.cc"
    dcomplex tmp2_ = stress[0*3+0];
    dcomplex tmp3_ = stress[0*3+2];
#line 361 "pmlelastic2d.cc"
    dcomplex tmp5_ = vx[0*2+0];
    dcomplex tmp6_ = vx[0*2+1];
#line 364 "pmlelastic2d.cc"
#line 365 "pmlelastic2d.cc"
    dcomplex tmp9_ = conj(tmp5_);
    dcomplex tmp10_ = conj(tmp6_);
    dcomplex tmp11_ = tmp2_ * tmp9_;
    dcomplex tmp12_ = tmp3_ * tmp10_;
    dcomplex tmp13_ = tmp11_ + tmp12_;
    dcomplex tmp14_ = tmp3_ * tmp9_;
    dcomplex tmp15_ = tmp2_ * tmp10_;
    dcomplex tmp16_ = tmp14_ + tmp15_;
    double tmp17_ = real(tmp13_);
    double tmp18_ = real(tmp16_);
    double tmp19_ = -tmp17_;
    double tmp20_ = -tmp18_;
    double tmp22_ = tmp19_ / 2.0;
    double tmp23_ = tmp20_ / 2.0;
#line 362 "pmlelastic2d.cc"
    EX[0*2+0] = tmp22_;
    EX[0*2+1] = tmp23_;
    } /* </generated matexpr> */
#line 367 "pmlelastic2d.cc"
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

