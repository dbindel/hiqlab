/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstdio>
#include <algorithm>

#include "luasupport.h"
#include "qmatrix.h"
#include "pmlelement.h"
#include "mesh.h"

#define ME PMLElement


ME::~ME()
{
    if (func_)
        delete func_;
    if (stretch_size)
        delete[] stretch;
}


void ME::set_stretch(PMLFunc* func)
{
    if (func_)
        delete func_;
    func_ = func;
}


void ME::set_stretch(lua_State* L, int func)
{
    set_stretch( new LuaPMLFunc(L, func) );
}


void ME::set_stretch(double* stretch, int ndm, int nen)
{
    stretch_size = ndm*nen;
    if (stretch_size == 0)
        this->stretch = stretch;
    else {
        this->stretch = new double[stretch_size];
        std::copy(stretch, stretch+stretch_size, this->stretch);
    }
}


void ME::compute_stretch(double* xx, int ndm)
{
    if (!func_) {
        for (int i = 0; i < ndm; ++i)
            xx[i] = 0;
        return;
    }

    if (func_)
        (*func_)(xx, ndm, 1);
}


void ME::compute_stretch(Mesh* mesh, int nodeid, double* xx, int ndm)
{
    if (stretch) {
        QMatrix<double> stretch(this->stretch, ndm,0);
        for (int j = 0; j < ndm; ++j)
            xx[j] = stretch(j,nodeid);

    } else {
        for (int j = 0; j < ndm; ++j)
            xx[j] = mesh->x(j,nodeid);
        compute_stretch(xx, ndm);
    }
}


void ME::set_local_stretch(Mesh* mesh, int eltid,
                           dcomplex* nodestretch1,
                           int nen, int ndm)
{
    QMatrix<dcomplex> nodestretch(nodestretch1, ndm,nen);

    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        double xx[3];
        compute_stretch(mesh, nodeid, xx, ndm);
        for (int i = 0; i < ndm; ++i)
            nodestretch(i,j) =  dcomplex(1.0, -xx[i]);
    }
}
