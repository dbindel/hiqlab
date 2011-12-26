/* HiQLab
 * Copyright (c): Regents of the University of California
 */

extern "C" {
  #include <lua.h>
  #include <tolua++.h>
}

#include <cstdio>
#include <cstring>
#include <cmath>

#include "qcomplex.h"
#include "leigs.h"
#include "luasupport.h"


// ---- Symmetric version ----


ArpackDSLua::ArpackDSLua()
{
}


ArpackDSLua::~ArpackDSLua()
{
}


void ArpackDSLua::compute_eigs(lua_State* L_, int op1_, int op2_, int opM_,
                               QArray* d, QArray* v)
{
    L = L_;
    op1 = op1_;
    op2 = op2_;
    opM = opM_;
    ArpackDS::compute_eigs(d->data_r(), (v == NULL ? NULL : v->data_r()));
}


void ArpackDSLua::times_OP1(double* x, double* opx)
{
    QArray y1(n,1, 0,1, x);
    QArray y2(n,1, 0,1, opx);
    lua_pushvalue(L, op1);
    tolua_pushusertype(L, (void*) &y1, "QArray");
    tolua_pushusertype(L, (void*) &y2, "QArray");
    lua_pcall2(L, 2, 0);
}


void ArpackDSLua::times_OP2(double* x, double* opx, double* Mx)
{
    QArray y1(n,1, 0,1, x);
    QArray y2(n,1, 0,1, opx);
    QArray y3(n,1, 0,1, Mx);
    lua_pushvalue(L, op2);
    tolua_pushusertype(L, (void*) &y1, "QArray");
    tolua_pushusertype(L, (void*) &y2, "QArray");
    tolua_pushusertype(L, (void*) &y3, "QArray");
    lua_pcall2(L, 3, 0);
}


void ArpackDSLua::times_M(double* x, double* Mx)
{
    QArray y1(n,1, 0,1, x);
    QArray y2(n,1, 0,1, Mx);
    lua_pushvalue(L, opM);
    tolua_pushusertype(L, (void*) &y1, "QArray");
    tolua_pushusertype(L, (void*) &y2, "QArray");
    lua_pcall2(L, 2, 0);
}


// ---- Nonsymmetric version ----


ArpackDNLua::ArpackDNLua()
{
}


ArpackDNLua::~ArpackDNLua()
{
}


void ArpackDNLua::compute_eigs(lua_State* L_, int op1_, int op2_, int opM_,
                               QArray* d, QArray* v)
{
    L = L_;
    op1 = op1_;
    op2 = op2_;
    opM = opM_;
    ArpackDN::compute_eigs(d->data_r(),
                           d->data_i(),
                           (v == NULL ? NULL : v->data_r()));
}


void ArpackDNLua::times_OP1(double* x, double* opx)
{
    QArray y1(n,1, 0,1, x);
    QArray y2(n,1, 0,1, opx);
    lua_pushvalue(L, op1);
    tolua_pushusertype(L, (void*) &y1, "QArray");
    tolua_pushusertype(L, (void*) &y2, "QArray");
    lua_pcall2(L, 2, 0);
}


void ArpackDNLua::times_OP2(double* x, double* opx, double* Mx)
{
    QArray y1(n,1, 0,1, x);
    QArray y2(n,1, 0,1, opx);
    QArray y3(n,1, 0,1, Mx);
    lua_pushvalue(L, op2);
    tolua_pushusertype(L, (void*) &y1, "QArray");
    tolua_pushusertype(L, (void*) &y2, "QArray");
    tolua_pushusertype(L, (void*) &y3, "QArray");
    lua_pcall2(L, 3, 0);
}


void ArpackDNLua::times_M(double* x, double* Mx)
{
    QArray y1(n,1, 0,1, x);
    QArray y2(n,1, 0,1, Mx);
    lua_pushvalue(L, opM);
    tolua_pushusertype(L, (void*) &y1, "QArray");
    tolua_pushusertype(L, (void*) &y2, "QArray");
    lua_pcall2(L, 2, 0);
}


// ---- Nonsymmetric complex version ----


ArpackZNLua::ArpackZNLua()
{
}


ArpackZNLua::~ArpackZNLua()
{
}


void ArpackZNLua::compute_eigs(lua_State* L_, int op1_, int op2_, int opM_,
                               QArray* d, QArray* v)
{
    L = L_;
    op1 = op1_;
    op2 = op2_;
    opM = opM_;
    ArpackZN::compute_eigs((dcomplex*) d->data_r(),
                           (dcomplex*) (v == NULL ? NULL : v->data_r()));
}


void ArpackZNLua::times_OP1(dcomplex* x, dcomplex* opx)
{
    QArray y1(n,1, 1,1, (double*) x);
    QArray y2(n,1, 1,1, (double*) opx);
    lua_pushvalue(L, op1);
    tolua_pushusertype(L, (void*) &y1, "QArray");
    tolua_pushusertype(L, (void*) &y2, "QArray");
    lua_pcall2(L, 2, 0);
}


void ArpackZNLua::times_OP2(dcomplex* x, dcomplex* opx, dcomplex* Mx)
{
    QArray y1(n,1, 1,1, (double*) x);
    QArray y2(n,1, 1,1, (double*) opx);
    QArray y3(n,1, 1,1, (double*) Mx);
    lua_pushvalue(L, op2);
    tolua_pushusertype(L, (void*) &y1, "QArray");
    tolua_pushusertype(L, (void*) &y2, "QArray");
    tolua_pushusertype(L, (void*) &y3, "QArray");
    lua_pcall2(L, 3, 0);
}


void ArpackZNLua::times_M(dcomplex* x, dcomplex* Mx)
{
    QArray y1(n,1, 1,1, (double*) x);
    QArray y2(n,1, 1,1, (double*) Mx);
    lua_pushvalue(L, opM);
    tolua_pushusertype(L, (void*) &y1, "QArray");
    tolua_pushusertype(L, (void*) &y2, "QArray");
    lua_pcall2(L, 2, 0);
}
