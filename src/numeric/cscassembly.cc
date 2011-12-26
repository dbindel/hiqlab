#include <cassert>
#include <algorithm>

#include "qcomplex.h"
#include "mesh.h"
#include "qassembly.h"
#include "cscmatrix.h"
#include "cscassembly.h"

using std::fill;

CSCAssembler::CSCAssembler(CSCMatrix* matrix)
{
    mm = matrix->get_m();
    nn = matrix->get_n();
    jc = matrix->get_jc();
    ir = matrix->get_ir();
    pr = matrix->get_Ax();
    pi = matrix->get_Az();
    int nnz = jc[nn];
    fill(pr, pr+nnz, 0);
    if (pi)
        fill(pi, pi+nnz, 0);
}


CSCAssembler::CSCAssembler(int* jc, int* ir, double* pr, double* pi) :
    jc(jc), ir(ir), pr(pr), pi(pi)
{
    int nnz = jc[nn];
    fill(pr, pr+nnz, 0);
    if (pi)
        fill(pi, pi+nnz, 0);
}


void CSCAssembler::add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke)
{
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; ++i)
            add(eltidm[i], eltidn[j], *Ke++);
}


void CSCAssembler::add(int* eltid, int n, dcomplex* Ke)
{
    add(eltid, n, eltid, n, Ke);
}


void CSCAssembler::add(int* eltidm, int m, int* eltidn, int n, double* Ke)
{
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; ++i)
            add(eltidm[i], eltidn[j], *Ke++);
}


void CSCAssembler::add(int* eltid, int n, double* Ke)
{
    add(eltid, n, eltid, n, Ke);
}


void CSCAssembler::add(int i, int j, dcomplex Ke)
{
    if (i < 0 || j < 0 || i >= mm || j >= nn)
        return;
    for (int ii = jc[j]; ii < jc[j+1]; ++ii) {
        if (ir[ii] == i) {
            pr[ii] += real(Ke);
            if (pi)
                pi[ii] += imag(Ke);
            return;
        }
    }
    assert(0);
}


void CSCAssembler::add(int i, int j, double Ke)
{
    if (i < 0 || j < 0 || i >= mm || j >= nn)
        return;
    for (int ii = jc[j]; ii < jc[j+1]; ++ii) {
        if (ir[ii] == i) {
            pr[ii] += Ke;
            return;
        }
    }
    assert(0);
}


void assemble_dR(Mesh* mesh, CSCMatrix* K, 
                 double cx, double cv, double ca,
                 int reduced)
{
    CSCAssembler assembler(K);
    mesh->assemble_dR(&assembler, cx, cv, ca, reduced);
}
