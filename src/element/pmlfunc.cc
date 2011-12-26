extern "C" {
#include "lua.h"
#include "lauxlib.h"
}

#include "luasupport.h"
#include "pmlfunc.h"
#include <cmath>


PMLFunc::~PMLFunc()
{
}


LuaPMLFunc::LuaPMLFunc(lua_State* L, int func) :
    L(L)
{
    lua_pushvalue(L, func);
    ref_ = luaL_ref(L, LUA_REGISTRYINDEX);
}


LuaPMLFunc::~LuaPMLFunc()
{
    luaL_unref(L, LUA_REGISTRYINDEX, ref_);
}


void LuaPMLFunc::operator()(double* xx, int ndm_, int npts)
{
    for (int j = 0; j < npts; ++j) {
        lua_rawgeti(L, LUA_REGISTRYINDEX, ref_);
        for (int i = 0; i < ndm_; ++i)
            lua_pushnumber(L, xx[i]);
        if (lua_pcall2(L, ndm_, ndm_) == 0) {
            for (int i = 0; i < ndm_; ++i)
                xx[i] = lua_tonumber(L, i-ndm_);
            lua_pop(L, ndm_);
        }
        xx += ndm_;
    }
}


template<class T>
inline T max(T x, T y)
{
    return x > y ? x : y;
}

void BoxPMLFunc::operator()(double* xx, int ndm, int npts)
{
    using namespace std;
    for (int j = 0; j < npts; ++j) {
        xx[0] = max( (abs(xx[0])-xpml_)/(xrad_-xpml_) * f0_, 0.0 );
        xx[1] = max( (abs(xx[1])-ypml_)/(yrad_-ypml_) * f0_, 0.0 );
        xx += ndm;
    }
}
