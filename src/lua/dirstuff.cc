#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include "dirstuff.h"


static int lua_cd(lua_State* L)
{
    int n = lua_gettop(L);
    if (n == 0) {
        char buf[512];
        lua_pushstring(L, getcwd(buf,sizeof(buf)));
        return 1;
    } else if (n == 1) {
        if (!lua_isstring(L,1)) {
            lua_pushstring(L, "Path name must be a string in cd()");
            lua_error(L);
        }
        int status = chdir(lua_tostring(L,1));
        if (status < 0) {
            lua_pushfstring(L, "Could not cd to '%s'", lua_tostring(L,1));
            lua_error(L);
        }
        return 0;
    } else {
        lua_pushstring(L, "cd takes zero or one argument");
        lua_error(L);
        return 0;
    }
}


static int lua_dir(lua_State* L)
{
    int n = lua_gettop(L);
    const char* dirname = ".";
    if (n > 2) {
        lua_pushstring(L, "Path name must be a string in cd()");
        lua_error(L);
    }
    if (n == 1) {
        if (!lua_isstring(L,1)) {
            lua_pushstring(L, "Path name must be a string in cd()");
            lua_error(L);
        }
        dirname = lua_tostring(L,1);
    }
    DIR* dir = opendir(dirname);
    if (!dir) {
        lua_pushfstring(L, "Could not open '%s'", dirname);
        lua_error(L);
    }
    lua_newtable(L);
    int t = lua_gettop(L);
    struct dirent* entry;
    for (int i = 1; entry = readdir(dir); ++i) {
        lua_pushstring(L, entry->d_name);
        lua_rawseti(L, t, i);
    }
    closedir(dir);
    return 1;
}


int lua_dirstuff_open(lua_State* L)
{
    const char *lscode =
        "function ls(name, pattern) "
        " table.foreach(dir(name or '.'), "
        "   function(k,s) "
        "     if not pattern or (pattern and string.find(s, pattern)) then "
        "       print(s) "
        "     end"
        "   end) "
        "end";
    lua_register(L, "cd",  lua_cd);
    lua_register(L, "dir", lua_dir);
    lua_dobuffer(L, lscode, strlen(lscode), "dir: embedded Lua code");
    return 1;
}
