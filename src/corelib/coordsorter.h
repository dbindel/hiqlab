/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef COORDSORTER_H
#define COORDSORTER_H

#include "mesh.h"

class CoordSorter {
public:
    CoordSorter(int ndm, Mesh& mesh, double tol) :
        ndm(ndm), mesh(mesh), tol(tol) {}

    int operator()(int i, int j)
    {
        for (int k = 0; k < ndm; ++k) {
            if (mesh.x(k,i) < mesh.x(k,j)-tol)   return  1;
            if (mesh.x(k,i) > mesh.x(k,j)+tol)   return  0;
        }
        return 0;
    }

    int compare(int i, int j)
    {
        for (int k = 0; k < ndm; ++k) {
            if (mesh.x(k,i) < mesh.x(k,j)-tol)   return  -1;
            if (mesh.x(k,i) > mesh.x(k,j)+tol)   return   1;
        }
        return 0;
    }

private:
    int ndm;
    Mesh& mesh;
    double tol;
};

#endif /* COORDSORTER_H */
