#ifndef PMLFUNC_H
#define PMLFUNC_H

extern "C" {
#include "lua.h"
}

// FIXME: ndm should be specified at function construct time!


class PMLFunc {
public:
    virtual ~PMLFunc();
    virtual void operator()(double* xx, int ndm, int npts = 1) = 0;
};


class LuaPMLFunc : public PMLFunc {
public:
    LuaPMLFunc(lua_State* L, int func);
    ~LuaPMLFunc();
    void operator()(double* xx, int ndm, int npts = 1);

private:
    lua_State* L;
    int ref_;
    // int ndm_; FIXME
};


class BoxPMLFunc : public PMLFunc {
public:
    BoxPMLFunc(double xrad, double yrad,
               double xpml, double ypml,
               double f0) :
        xrad_(xrad), yrad_(yrad),
        xpml_(xpml), ypml_(ypml),
        f0_(f0) {}

    void operator()(double* xx, int ndm, int npts = 1);

private:
    double xrad_, yrad_;
    double xpml_, ypml_;
    double f0_;
};


#endif /* PMLFUNC_H */
