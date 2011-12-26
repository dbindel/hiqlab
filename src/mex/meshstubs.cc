/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <mex.h>
#include <iostream>

#include "mexutil.h"
#include "mesh.h"
#include "cscindexer.h"
#include "mesh_csc_dR.h"
#include "qassembly.h"
#include "cscassembly.h"
#include "pzlinear.h"
#include "tolua++.h"

#include <vector>
#include <memory>

using std::vector;
using std::auto_ptr;
using std::copy;
using std::fill;


mxArray* Mesh_assemble_struct1(Mesh* m, int reduced)
{
    vector<int> jc;
    vector<int> ir;
    int n = reduced ? m->get_numid() : m->numnp() * m->get_ndf();
    build_csc_matrix(m, n, jc, ir, reduced);
    int nnz = jc[n];
    mxArray* K = mxCreateSparse(n, n, nnz, mxREAL);
    copy(jc.begin(), jc.end(), mxGetJc(K));
    copy(ir.begin(), ir.end(), mxGetIr(K));
    fill(mxGetPr(K), mxGetPr(K)+nnz, 1);
    return K;
}


mxArray* Mesh_assemble_dR1(Mesh* m, const mxArray* Kmat, 
                           double cx, double cv, double ca, int reduced)
{
    auto_ptr<CSCMatrix> K( mx_to_CSCMatrix(Kmat) );
    assemble_dR(m, K.get(), cx, cv, ca, reduced);
    return CSCMatrix_to_mx(K.get());
}


mxArray* Mesh_assemble_dR1(Mesh* m, 
                           double cx, double cv, double ca, int reduced)
{
    int N = (reduced ? m->get_numid() : m->numnp()*m->get_ndf());
    auto_ptr<CoordMatrix> assembler(new CoordMatrix(N));
    m->assemble_dR(assembler.get(), cx, cv, ca, reduced);

    assembler->pack();
    int n   = assembler->get_N();
    int nnz = assembler->get_ncoord();
    mxArray* A = mxCreateSparse(n, n, nnz, mxCOMPLEX);
    assembler->to_sparse(mxGetJc(A), mxGetIr(A), mxGetPr(A), mxGetPi(A));

    return A;
}


mxArray* Mesh_element_dR1(Mesh* m, int eltid, double cx, double cv, double ca)
{
    auto_ptr<CSCMatrix> K( element_dR(m, eltid, cx, cv, ca) );
    return CSCMatrix_to_mx(K.get());
}


void Mesh_get_e1(Mesh* mesh, int* e)
{
    int nen = mesh->get_nen();
    int numelt = mesh->numelt();
    for (int j = 0; j < numelt; ++j) {
        int jnen = mesh->get_nen(j);
        for (int i = 0; i < jnen; ++i)
            e[j*nen+i] = mesh->ix(i,j);
        for (int i = jnen; i < nen; ++i)
            e[j*nen+i] = -1;
    }
}


void Mesh_get_bc1(Mesh* mesh, double* bc)
{
    int ndf = mesh->get_ndf();
    int numnp = mesh->numnp();
    for (int j = 0; j < numnp; ++j)
        for (int i = 0; i < ndf; ++i) {
            char bcij = mesh->bc(i,j);
            if (bcij == 'u')
                *bc++ = 1;
            else if (bcij == 'f')
                *bc++ = 2;
            else
                *bc++ = 0;
        }
}


void Mesh_get_u1(Mesh* mesh, double* u, double* ui)
{
    int n = mesh->numnp();
    int m = mesh->get_ndf();
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            *u++  = mesh->u(i,j);
            *ui++ = mesh->ui(i,j);
        }
    }
}


void Mesh_get_v1(Mesh* mesh, double* u)
{
    int n = mesh->numnp();
    int m = mesh->get_ndf();
    for (int j = 0; j < n; ++j) 
        for (int i = 0; i < m; ++i) 
            *u++  = mesh->v(i,j);
}


void Mesh_get_a1(Mesh* mesh, double* u)
{
    int n = mesh->numnp();
    int m = mesh->get_ndf();
    for (int j = 0; j < n; ++j) 
        for (int i = 0; i < m; ++i) 
            *u++  = mesh->a(i,j);
}


void Mesh_get_f1(Mesh* mesh, double* f, double* fi)
{
    int m = mesh->get_ndf();
    int n = mesh->numnp();
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; ++i) {
            *f++  = mesh->f(i,j);
            *fi++ = mesh->fi(i,j);
        }
}


void Mesh_get_nen_elt(Mesh* m, int* result)
{
    for (int i = 0; i < m->numelt(); ++i)
        result[i] = m->get_nen(i);
}
