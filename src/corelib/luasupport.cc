/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include "luasupport.h"
#include <stdio.h>

/*@T
 * \section{Lua support routines}
 *
 * We provide a very small number of extensions to the support library
 * in Lua for general use.  Because these extensions are used by the
 * mesh object and by some of the elements, we put them here rather than
 * in the Lua subdirectory.
 *
 *@q*/


/*@T
 * \subsection{Protected calls}
 *
 * The [[lua_pcall]] routine makes a protected call -- that is, it
 * either completes successfully (status code zero) or pushes a message
 * on the stack and returns an error code.  It does not throw an exception.
 * The [[lua_pcall2]] routine is like [[lua_pcall]], except that it
 * makes an effort to report a problem in the case of an error.
 *@c*/

int lua_pcall2(lua_State* L, int nargin, int nargout)
{
    int base = lua_gettop(L) - nargin;
    lua_pushliteral(L, "_TRACEBACK");
    lua_rawget(L, LUA_GLOBALSINDEX);
    lua_insert(L, base);
    int status = lua_pcall(L, nargin, nargout, base);
    int had_traceback = lua_isfunction(L,base);
    lua_remove(L, base);
    if (status == 0)
        return 0;

    // If there was an error, try to handle it with _ALERT
    lua_getglobal(L, "_ALERT");
    if (lua_isfunction(L,-1)) {
        lua_pushvalue(L,-2);
        if (lua_pcall(L,1,0,0) == 0) {
            lua_pop(L,1);
            return status;
        } else {
            lua_pop(L,1);
        }
    } else {
        lua_pop(L,1);
    }

    // If these doesn't work, just printf
    const char* msg = lua_tostring(L, -1);
    fprintf(stderr, "%s\n", msg);
    return status;
}


/*@T
 * \subsection{Object data association}
 * 
 * The Lua registry can be used by the C code to store various data
 * (particularly things like Lua function handles, which can't be
 * stored directly in a C type).  The [[lua_objtable]] routine
 * retrieves a table in the registry (or creates a new one) keyed by
 * an object's [[this]] pointer.  The [[lua_getobjtable]] and
 * [[lua_setobjtable]] routines get and set entries in the table
 * associated with an object in this way.
 *
 * Note 1: Right now there's a possibility of a memory leak if the
 * destructor of an object using this mechanism does not explicitly
 * set its registry entry to NULL on destruction!
 *
 * Note 2: This is delicate when different levels in a multiple
 * inheritance hierarchy want to access the same data.  Don't do it.
 *
 *@c*/


void lua_objtable(lua_State* L, void* p)
{
    lua_pushlightuserdata(L, p);
    lua_gettable(L, LUA_REGISTRYINDEX);
    if (lua_isnil(L,-1)) {
        lua_pop(L,1);
        lua_newtable(L);
        lua_pushlightuserdata(L, p);
        lua_pushvalue(L,-2);
        lua_settable(L, LUA_REGISTRYINDEX);
    }
}


void lua_objgettable(lua_State* L, void* p)
{
    lua_objtable(L, p);
    lua_insert(L,-2);
    lua_gettable(L,-2);
    lua_remove(L,-2);
}


void lua_objsettable(lua_State* L, void* p)
{
    lua_objtable(L, p);
    lua_insert(L,-3);
    lua_settable(L,-3);
    lua_pop(L,1);
}
