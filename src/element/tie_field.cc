/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "luasupport.h"
#include "qmatrix.h"
#include "tie_field.h"
#include "mesh.h"

#define ME TieFieldElement


ME::ME(lua_State* L, int func) :
    Element(0), L(L)
{
    lua_pushlightuserdata(L, this);
    if (!lua_isfunction(L,func)) {
        this->L = NULL;
        lua_pushnil(L);
    } else {
        this->L = L;
        lua_pushvalue(L,func);
    }
    lua_settable(L, LUA_REGISTRYINDEX);
}


ME::~ME()
{
    // Clean up the Lua table
    lua_pushlightuserdata(L, this);
    lua_pushnil(L);
    lua_settable(L, LUA_REGISTRYINDEX);
}


void ME::initialize(Mesh* mesh, int eltid)
{
    call_func(mesh, eltid, 0, 0, NULL);
    mesh->nbranch(eltid) = nbranch+1;
}


void ME::assign_ids(Mesh* mesh, int eltid)
{
    call_func(mesh, eltid, 1, 0, NULL);
    for (int i = 0; i <= nbranch; ++i)
        mesh->branchid(i,eltid) = 1;
}


void ME::assemble_R(Mesh* mesh, int eltid)
{
    call_func(mesh, eltid, 2, 0, NULL);
}


void ME::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                                  double cx, double cv, double ca)
{
    call_func(mesh, eltid, 3, cx, K);
}


void ME::call_func(Mesh* mesh, int eltid, int mode,
                   double cx, QAssembler* K)
{
    // Value returned from function should not be nondimensionalized??
    int N      = mesh->numnp();
    int maxndf = mesh->get_ndf();
    int ndm    = mesh->get_ndm();

    nbranch = 0;

    lua_pushlightuserdata(L, this);
    lua_gettable(L, LUA_REGISTRYINDEX);
    int func = lua_gettop(L);

    if (!lua_isfunction(L,func)) {
        lua_pop(L,1);
        return;
    }

    for (int j = 0; j < N; ++j) {
        int t = lua_gettop(L)+1;

        // Make the call
        lua_pushvalue(L, func);
        for (int i = 0; i < ndm; ++i)
            lua_pushnumber(L, mesh->x(i,j));

        if (lua_pcall2(L, ndm, maxndf+1) == 0) {

            // And process the return values
            const char* s = lua_tostring(L, t++);
            if (s) {
                for (unsigned i = 0; i < strlen(s); ++i) {
                    if (s[i] == 'u') {

                        double v;
                        if (mode != 0) {
                            v = lua_tonumber(L, t++);
//                            v = lua_tonumber(L, t++);
                        }
                        if (mode == 0) { // Find number of potential dofs
                            ++nbranch;
                        } else if (mode == 1 && mesh->id(i,j)>0 ) { // Assignids
                            ++nbranch;
                        } else if (mode == 2 && mesh->id(i,j)>=0 ) { // Get R
                            ++nbranch;
                            double s = mesh->branchu(0,      eltid);
                            double u = mesh->u(i,j);
                            double l = mesh->branchu(nbranch,eltid);

                            mesh->branchf(0,eltid)       +=  l*v;
                            mesh->f(i,j)                 += -l;
                            mesh->branchf(nbranch,eltid) +=  s*v-u;

                        } else if (mode == 3 && mesh->id(i,j) >=0 ) { // Get K
                            ++nbranch;
                            int id[3];
                            id[0] = mesh->ibranch(0,      eltid);
                            id[1] = mesh->inode(i,j);
                            id[2] = mesh->ibranch(nbranch,eltid);

                            double Ke[9] = {
                                     0,      0,    v*cx,
                                     0,      0,     -cx,
                                   v*cx,    -cx,     0 };

                            K->add(id, 3, Ke);

                        }
                    }
                }
            }

            // Finally, clean up
            lua_pop(L, maxndf+1);

        }
    }

    lua_pop(L,1);
}
