/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef LUASUPPORT_H
#define LUASUPPORT_H

extern "C" {
  #include <lua.h>
}

int lua_pcall2(lua_State* L, int nargin, int nargout);
void lua_objgettable(lua_State* L, void* p);
void lua_objsettable(lua_State* L, void* p);

#endif /* LUASUPPORT_H */
