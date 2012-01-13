#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>

#include "qassembly.h"
#include "mesh.h"
#include "csrindexer.h"

using std::vector;

void csr_sparsity_count(Mesh* mesh, int n, int* counts, int reduced)
{
    CSRSparsityCounter counter(n);
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


void csr_build_index(Mesh* mesh, int n, int* ir, int* jc, int reduced)
{
    CSRIndexBuilder indexer(n, ir, jc);
    mesh->assemble_struct(&indexer, reduced);
}


/*@T
 * \subsubsection{Putting it all together}
 *
 *@c*/
void build_csr_matrix(Mesh* mesh, int n, 
                      vector<int>& ir, vector<int>& jc,
                      int reduced)
{
    ir.resize(n+1);
    csr_sparsity_count(mesh, n, &ir[0], reduced);
    int nnz = ir[n];
    jc.resize(nnz);
    csr_build_index(mesh, n, &ir[0], &jc[0], reduced);
}
