/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstdio>

#include "fieldeval.h"
#include "mesh.h"


FieldEval::~FieldEval()
{
}


FieldEvalPower::~FieldEvalPower()
{
}


void FieldEvalPower::operator()(Mesh* mesh, int eltid, double* X, double* fX)
{
    mesh->etype(eltid)->mean_power(mesh, eltid, X, fX);
}
