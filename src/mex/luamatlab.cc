/* HiQLab
 * Copyright (c): Regents of the University of California
 */

extern "C" {
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

#include <tolua++.h>
#include <mex.h>
#include <algorithm>

#include "luamex.h"
#include "qhelperslua.h"
#include "qassemblylua.h"
#include "qarraylua.h"
#include "cscmatrixlua.h"
#include "meshlua.h"
#include "elementlua.h"
#include "shapeslua.h"
#include "material_modellua.h"
#include "tedlinearlua.h"
#include "pzlinearlua.h"

#include "qarray.h"
#include "qassembly.h"
#include "cscmatrix.h"
#include "element.h"
#include "mesh.h"
#include "luasupport.h"

#include "luamatlab.h"


lua_State* Lua_open()
{
    lua_State* state = lua_open();
    luaopen_base(state);
    luaopen_io(state);
    luaopen_math(state);
    luaopen_table(state);
    luaopen_string(state);
    luaopen_debug(state);
    tolua_open(state);
    tolua_qhelpers_open(state);
    tolua_qassembly_open(state);
    tolua_qarray_open(state);
    tolua_cscmatrix_open(state);
    tolua_mesh_open(state);
    tolua_element_open(state);
    tolua_shapes_open(state);
    tolua_luamex_open(state);
    tolua_material_model_open(state);
    tolua_tedlinear_open(state);
    tolua_pzlinear_open(state);
    return state;
}


int Lua_dofunc(int status, lua_State* L, const char* name)
{
    if (status != 0) {
        const char* msg = lua_tostring(L, -1);
        if (msg == NULL)
            msg = "(error with no message)";
        const char* type = "Unknown error";
        if (status == LUA_ERRSYNTAX) type = "Syntax error";
        if (status == LUA_ERRMEM)    type = "Memory error";
        mexPrintf("%s loading %s: %s\n", type, name, msg);
        lua_pop(L,1);
        mexErrMsgTxt("Error in Lua execution");
    }
    return lua_pcall2(L,0,0);
}


int Lua_dofile(lua_State* L, const char* fname)
{
    return Lua_dofunc(luaL_loadfile(L, fname), L, fname);
}


int Lua_dostring(lua_State* L, const char* s)
{
    return Lua_dofunc(luaL_loadbuffer(L, s, strlen(s), "[User]"), L, "[User]");
}


void Lua_set_mesh(lua_State* L, const char* name, Mesh* mesh)
{
    tolua_pushusertype(L, mesh, "Mesh");
    lua_setglobal(L, name);
}


Mesh* Lua_get_mesh(lua_State* L, const char* name)
{
    tolua_Error tolua_err;
    Mesh* mesh = NULL;
    lua_getglobal(L, name);
    if (tolua_isusertype(L,-1,"Mesh",0,&tolua_err))
        mesh = (Mesh*) tolua_tousertype(L,-1,0);
    lua_pop(L,1);
    return mesh;
}


void Lua_set_string(lua_State* L, const char* name, const char* s)
{
    lua_pushstring(L, s);
    lua_setglobal(L, name);
}


mxArray* Lua_get_string(lua_State* L, const char* name)
{
    lua_getglobal(L, name);
    if (!lua_isstring(L,-1)) {
        lua_pop(L, 1);
        mexErrMsgTxt("Variable does not exist");
    }
    mxArray* s = mxCreateString(lua_tostring(L, -1));
    lua_pop(L,1);
    return s;
}


void Lua_set_double(lua_State* L, const char* name, double x)
{
    lua_pushnumber(L, x);
    lua_setglobal(L, name);
}


double Lua_get_double(lua_State* L, const char* name)
{
    lua_getglobal(L, name);
    if (!lua_isnumber(L,-1)) {
        lua_pop(L, 1);
        mexErrMsgTxt("Variable does not exist");
    }
    double x = lua_tonumber(L, -1);
    lua_pop(L, 1);
    return x;
}


void Lua_set_cscmatrix(lua_State* L, const mxArray* array)
{
    if (mxGetM(array) != mxGetN(array))
        mexErrMsgTxt("Non-square sparse arrays currently not supported");
    int n = mxGetN(array);
    CSCMatrix* mat = new CSCMatrix(mxGetJc(array), mxGetIr(array),
                                   mxGetPr(array), mxGetPi(array), n, n);
    tolua_pushusertype_and_takeownership(L, mat, "CSCMatrix");
}


void Lua_set_qarrayr(lua_State* L, const mxArray* array, int a_type)
{
    int size = mxGetM(array)*mxGetN(array);
    QArray* m = new QArray(mxGetM(array), mxGetN(array), a_type);
    std::copy(mxGetPr(array), mxGetPr(array)+size, m->data_r());
    tolua_pushusertype_and_takeownership(L, m, "QArray");
}


void Lua_set_qarrayz(lua_State* L, const mxArray* array)
{
    int size = mxGetM(array)*mxGetN(array);
    QArray* m = new QArray(mxGetM(array), mxGetN(array), 2);
    std::copy(mxGetPr(array), mxGetPr(array)+size, m->data_r());
    std::copy(mxGetPi(array), mxGetPi(array)+size, m->data_i());
    tolua_pushusertype_and_takeownership(L, m, "QArray");
}


void Lua_set_array(lua_State* L, const char* name, const mxArray* array)
{
    if (mxIsSparse(array))
        Lua_set_cscmatrix(L, array);
    else if (mxIsDouble(array)) {
        if (mxGetPi(array)) 
            Lua_set_qarrayz(L, array);
        else
            Lua_set_qarrayr(L, array, 0);
    }
    lua_setglobal(L, name);
}


void Lua_set_arrayz(lua_State* L, const char* name, const mxArray* array)
{
    if (mxIsSparse(array))
        Lua_set_cscmatrix(L, array);
    else if (mxIsDouble(array)) {
        if (mxGetPi(array)) 
            Lua_set_qarrayz(L, array);
        else
            Lua_set_qarrayr(L, array, 2);
    }
    lua_setglobal(L, name);
}


mxArray* Lua_get_cscmatrix(lua_State* L)
{
    using std::copy;
    CSCMatrix* matrix = (CSCMatrix*) tolua_tousertype(L,-2,0);
    int n   = matrix->get_n();
    int nnz = matrix->get_nnz();
    int is_real = (matrix->get_Az() == 0);
    mxArray* A = mxCreateSparse(n, n, nnz, is_real ? mxREAL : mxCOMPLEX);
    copy(matrix->get_jc(), matrix->get_jc()+n+1, mxGetJc(A));
    copy(matrix->get_ir(), matrix->get_ir()+nnz, mxGetIr(A));
    copy(matrix->get_Ax(), matrix->get_Ax()+nnz, mxGetPr(A));
    if (!is_real)
        copy(matrix->get_Ax(), matrix->get_Ax()+nnz, mxGetPi(A));
    return A;
}


mxArray* Lua_get_qarray(lua_State* L)
{
    QArray* matrix = (QArray*) tolua_tousertype(L,-2,0);
    int m = matrix->m();
    int n = matrix->n();
    int b = matrix->base();

    mxArray* A;
    if (matrix->type() == 0) {
        A = mxCreateDoubleMatrix(m, n, mxREAL);
        double* pr = mxGetPr(A);
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i)
                *pr++ = matrix->get(i+b,j+b);
    } else {
        A = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
        double* pr = mxGetPr(A);
        double* pi = mxGetPi(A);
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i) {
                *pr++ = matrix->get (i+b,j+b);
                *pi++ = matrix->geti(i+b,j+b);
            }
    }
    return A;
}


mxArray* Lua_get_array(lua_State* L, const char* name)
{
    lua_getglobal(L, name);
    const char* t = tolua_typename(L,-1);
    mxArray* A = NULL;
    if (strcmp(t,"CSCMatrix") == 0)
        A = Lua_get_cscmatrix(L);
    else if (strcmp(t,"QArray") == 0)
        A = Lua_get_qarray(L);
    else {
        lua_pop(L,1);
        mexErrMsgTxt("Variable does not exist");
    }
    lua_pop(L,1);
    return A;
}


void luamex_print(const char* s)
{
    mexPrintf("%s", s);
}
