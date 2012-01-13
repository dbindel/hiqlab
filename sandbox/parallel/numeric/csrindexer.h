#ifndef CSRINDEXER_H
#define CSRINDEXER_H

#include "mesh.h"
#include "cscmatrix.h"
#include "cscindexer.h"

#include <vector>

/*@T
 * \subsection{Accumulating row counts}
 *
 * In the first pass through the matrix, we build the row index
 * structure for a sparse matrix using linked lists (instead of using
 * static arrays as in the usual compressed sparse forms).  We
 * allocate the link in chunks in order to amortize both the time and
 * space overhead of the calls to [[new]].  The linked list
 * representation should take roughly as much space as the full
 * compressed sparse representation -- the space taken by the pointers
 * is about the same as the space that would be taken by the matrix
 * entries.
 *
 * The time and space requirements can probably be improved by storing
 * the finite element matrix in blocked forms, with bitmaps representing
 * small subblocks.  But this is clearly not worthwhile unless the current
 * code proves to be a bottleneck.
 *
 * The [[CSRSparsityCounter]] object keeps all the state for this first
 * pass through the matrix.
 *
 * Essentially this is a transposed version of cscindexer
 *
 *@c*/
#define ENTRY_POOL_LEN 1024


class CSRSparsityCounter : public CSCSparsityCounter {
public:
    CSRSparsityCounter(int m, int n) : CSCSparsityCounter(n,m) {}
    CSRSparsityCounter(int n) : CSCSparsityCounter(n) {}
    ~CSRSparsityCounter() {}

    void add(int i, int j) {CSCSparsityCounter::add(j,i);}
};

/*@T
 * \subsection{Building the index structure}
 *
 * After accumulating the row counts, we want to actually build the
 * index structure.  We do this with a second structural assembly pass,
 * because (a) it keeps the memory overhead down, and (b) it's easy.
 * If this phase ever becomes the bottleneck, there are probably other
 * things we could do that would avoid the O(row_nnz^2) search/insert
 * overhead.
 *
 * The [[CSRIndexBuilder]] object is used to track this second stage.
 * I assume the [[ir]] array is already set up when the object is
 * constructed.
 * 
 *@c*/
class CSRIndexBuilder : public CSCIndexBuilder {
public:
    CSRIndexBuilder(int m, int n, int* ir, int* jc) : CSCIndexBuilder(n,m,ir,jc) {}
    CSRIndexBuilder(int n, int* ir, int* jc) : CSCIndexBuilder(n,ir,jc) {}
    void add(int i, int j) {CSCIndexBuilder::add(j,i);}
};

void build_csr_matrix(Mesh* mesh, int n, 
                      std::vector<int>& ir, std::vector<int>& jc, 
                      int reduced);

#endif /* CSRINDEXER_H */
