#ifndef LEIGS_H
#define LEIGS_H

#include "qarray.h"
#include "areigs.h"


class ArpackDSLua : public ArpackDS {
public:
    ArpackDSLua();
    ~ArpackDSLua();
    void compute_eigs(lua_State* L, int op1, int op2, int opM,
                      QArray* d, QArray* v);
private:
    lua_State* L;
    int op1, op2, opM;
    void times_OP1(double* x, double* opx);
    void times_OP2(double* x, double* opx, double* Mx);
    void times_M  (double* x, double* Mx);
};


class ArpackDNLua : public ArpackDN {
public:
    ArpackDNLua();
    ~ArpackDNLua();
    void compute_eigs(lua_State* L, int op1, int op2, int opM,
                      QArray* d, QArray* v);
private:
    lua_State* L;
    int op1, op2, opM;
    void times_OP1(double* x, double* opx);
    void times_OP2(double* x, double* opx, double* Mx);
    void times_M  (double* x, double* Mx);
};


class ArpackZNLua : public ArpackZN {
public:
    ArpackZNLua();
    ~ArpackZNLua();
    void compute_eigs(lua_State* L, int op1, int op2, int opM,
                      QArray* d, QArray* v);
private:
    lua_State* L;
    int op1, op2, opM;
    void times_OP1(dcomplex* x, dcomplex* opx);
    void times_OP2(dcomplex* x, dcomplex* opx, dcomplex* Mx);
    void times_M  (dcomplex* x, dcomplex* Mx);
};


#endif /* LEIGS_H */
