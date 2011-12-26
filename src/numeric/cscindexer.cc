#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>

#include "qassembly.h"
#include "mesh.h"
#include "cscindexer.h"

using std::vector;
using std::fill;

CSCSparsityCounter::CSCSparsityCounter(int m, int n) : 
    m(m), n(n), col(n)
{
    pool = new EntryPool(0);
    fill(col.begin(), col.end(), static_cast<Entry*>(0));
}

CSCSparsityCounter::CSCSparsityCounter(int n) : 
    m(n), n(n), col(n)
{
    pool = new EntryPool(0);
    fill(col.begin(), col.end(), static_cast<Entry*>(0));
}


CSCSparsityCounter::~CSCSparsityCounter()
{
    while (pool) {
        EntryPool* next = pool->next;
        delete pool;
        pool = next;
    }
}


CSCSparsityCounter::Entry* CSCSparsityCounter::new_entry()
{
    if (pool->used == ENTRY_POOL_LEN)
        pool = new EntryPool(pool);
    return &(pool->data[pool->used++]);
}


void CSCSparsityCounter::add(int i, int j)
{
    if (i < 0 || j < 0 || i >= m || j >= n)
        return;
    Entry** p;
    for (p = &col[j]; *p && (*p)->row < i; p = &((*p)->next));
    if (*p == NULL || (*p)->row != i) {
        Entry* e = new_entry();
        e->row = i;
        e->next = *p;
        *p = e;
    }
}


void CSCSparsityCounter::get_counts(int* counts)
{
    for (int j = 0; j < n; ++j) {
        int count = 0;
        for (Entry* p = col[j]; p; p = p->next, ++count);
        counts[j] = count;
    }
}


void csc_sparsity_count(Mesh* mesh, int n, int* counts, int reduced)
{
    CSCSparsityCounter counter(n);
    mesh->assemble_struct(&counter, reduced);
    counter.get_counts(counts);
    int nnz = 0;
    for (int j = 0; j < n; ++j) {
        int nnzj = counts[j];
        counts[j] = nnz;
        nnz += nnzj;
    }
    counts[n] = nnz;
}


CSCIndexBuilder::CSCIndexBuilder(int m, int n, int* jc, int* ir) :
    m(m), n(n), jc(jc), ir(ir)
{
    fill(ir, ir+jc[n], -1);
}


CSCIndexBuilder::CSCIndexBuilder(int n, int* jc, int* ir) :
    m(n), n(n), jc(jc), ir(ir)
{
    fill(ir, ir+jc[n], -1);
}


void CSCIndexBuilder::add(int i, int j)
{
    if (i < 0 || j < 0 || i >= m || j >= n)
        return;

    // Find the location where entry i should go in sorted order
    int ii;
    for (ii = jc[j]; ir[ii] != -1 && ir[ii] < i; ++ii);
    assert(ii < jc[j+1]);
    if (ir[ii] == i)
        return;

    // Shift all subsequent entries to make room for i, then insert
    for (int kk = jc[j+1]-1; kk > ii; --kk)
        ir[kk] = ir[kk-1];
    ir[ii] = i;
}


void csc_build_index(Mesh* mesh, int n, int* jc, int* ir, int reduced)
{
    CSCIndexBuilder indexer(n, jc, ir);
    mesh->assemble_struct(&indexer, reduced);
}


/*@T
 * \subsubsection{Putting it all together}
 *
 *@c*/
void build_csc_matrix(Mesh* mesh, int n, 
                      vector<int>& jc, vector<int>& ir,
                      int reduced)
{
    jc.resize(n+1);
    csc_sparsity_count(mesh, n, &jc[0], reduced);
    int nnz = jc[n];
    ir.resize(nnz);
    csc_build_index(mesh, n, &jc[0], &ir[0], reduced);
}


CSCMatrix* build_csc_matrix(Mesh* mesh, int reduced)
{
    int n = mesh->get_numid();
    CSCMatrix* result = new CSCMatrix(n, n, 0);
    int* jc = result->get_jc();
    csc_sparsity_count(mesh, n, jc, reduced);
    int nnz = jc[n];
    result->set_nnz(nnz);
    int* ir = result->get_ir();
    csc_build_index(mesh, n, jc, ir, reduced);
    return result;
}

