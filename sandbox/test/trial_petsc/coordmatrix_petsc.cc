/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 * $Id: coordmatrix_petsc.cc,v 1.3 2006/06/18 04:07:23 tkoyama Exp $
 */
 
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>

#include "coordmatrix_petsc.h"

#define MAXNEN 64

using namespace std;

CoordMatrix_Petsc::CoordMatrix_Petsc(int M, int N) :
    CoordMatrix(M,N)
{
}

CoordMatrix_Petsc::CoordMatrix_Petsc(int N) :
    CoordMatrix(N)
{
}

CoordMatrix_Petsc::~CoordMatrix_Petsc()
{
}

void CoordMatrix_Petsc::set_mapping_functions(int (* f_r)(int), int (* f_c)(int) )
{
    f_row = f_r;
    f_col = f_c;
}

void CoordMatrix_Petsc::add(int* eltid, int n, dcomplex* Ke)
{
    int eltid_m[n];
    for (int j = 0; j < n; ++j)
        if ( f_r(eltid[j]) >= min_r && f_r(eltid[j]) <= max_r ) {
            eltid_m[j] = eltid[j];
        } else
           eltid_m[j] = -3;
    CoordMatrix::add(&eltid_m, n, Ke);
}

void CoordMatrix_Petsc::add(int* eltid, int n, double* Ke)
{
    int eltid_m[n];
    for (int j = 0; j < n; ++j)
        if ( f_r(eltid[j]) >= min_r && f_r(eltid[j]) <= max_r ) {
            eltid_m[j] = eltid[j];
        } else
           eltid_m[j] = -3;
    CoordMatrix::add(&eltid_m, n, Ke);
}
