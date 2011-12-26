/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef LUAMATLAB_H
#define LUAMATLAB_H

extern "C" {
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

#include <mex.h>

#include "qarray.h"
#include "qassembly.h"
#include "element.h"
#include "mesh.h"
#include "material_model.h"

lua_State* Lua_open();
int        Lua_dofile(lua_State* L, const char* fname);
int        Lua_dostring(lua_State* L, const char* s);

void     Lua_set_mesh   (lua_State* L, const char* name, Mesh* mesh);
Mesh*    Lua_get_mesh   (lua_State* L, const char* name);
void     Lua_set_string (lua_State* L, const char* name, const char* s);
mxArray* Lua_get_string (lua_State* L, const char* name);
void     Lua_set_double (lua_State* L, const char* name, double x);
double   Lua_get_double (lua_State* L, const char* name);

void     Lua_set_array (lua_State* L, const char* name, const mxArray* array);
void     Lua_set_arrayz(lua_State* L, const char* name, const mxArray* array);
mxArray* Lua_get_array (lua_State* L, const char* name);

void luamex_print(const char* s);

#endif /* LUAMATLAB_H */
