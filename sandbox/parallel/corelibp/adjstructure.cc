#include <vector>
#include <iostream>

#include "adjstructure.h"
#include "cscindexer.h"
#include "metispart.h"

using std::vector;

/*@t --------
 * \subsection{Sparsity counter for adjacency computation}
 *
 *@c*/

void adj_sparsity_count(Mesh* mesh, int* counts)
{
    int n = mesh->numnp();
    CSCSparsityCounter counter(n);

    int numelt = mesh->numelt();
    for (int ie = 0; ie < numelt; ++ie) {

        int nen = mesh->get_nen(ie);

        // -- Add upper half
        for (int j = 0; j < nen; ++j)
            for (int i = j+1; i < nen; ++i) {
                counter.add(mesh->ix(i,ie), mesh->ix(j,ie));
                counter.add(mesh->ix(j,ie), mesh->ix(i,ie));
            }
    }

    counter.get_counts(counts);
    int nnz = 0;
    for (int j = 0; j < n; ++j) {
        int nnzj = counts[j];
        counts[j] = nnz;
        nnz += nnzj;
    }
    counts[n] = nnz;
}


void adj_build_index(Mesh* mesh, int* jc, int* ir)
{
    int n = mesh->numnp();
    CSCIndexBuilder indexer(n, jc, ir);

    int numelt = mesh->numelt();
    for (int ie = 0; ie < numelt; ++ie) {

        int nen = mesh->get_nen(ie);

        // -- Add upper half
        for (int j = 0; j < nen; ++j)
            for (int i = j+1; i < nen; ++i) {
                indexer.add(mesh->ix(i,ie), mesh->ix(j,ie));
                indexer.add(mesh->ix(j,ie), mesh->ix(i,ie));
            }
    }
}


void build_adj(Mesh* mesh, std::vector<int>& jc, std::vector<int>& ir)
{
    int n = mesh->numnp();
    jc.resize(n+1);
    adj_sparsity_count(mesh, &jc[0]);
    int nnz = jc[n];
    ir.resize(nnz);
    adj_build_index(mesh, &jc[0], &ir[0]);
}


void build_adj(Mesh* mesh, int* jc, int* ir)
{
    vector<int> jcv;
    vector<int> irv;

    build_adj(mesh,jcv,irv);
    int n   = mesh->numnp();
    int nnz = irv.size();

    for (int i = 0; i < n+1; ++i)
        jc[i] = jcv[i];
    for (int i = 0; i < nnz; ++i)
        ir[i] = irv[i];
}


int build_adj_nnz(Mesh* mesh)
{
    vector<int> jcv;
    vector<int> irv;

    build_adj(mesh,jcv,irv);
    int nnz = irv.size();

    return nnz;
}
