#ifndef CSCINDEXER_H
#define CSCINDEXER_H

#include "mesh.h"
#include "cscmatrix.h"

#include <vector>

/*@T
 * \subsection{Accumulating column counts}
 *
 * In the first pass through the matrix, we build the column index
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
 * The [[CSCSparsityCounter]] object keeps all the state for this first
 * pass through the matrix.
 *
 *@c*/
#define ENTRY_POOL_LEN 1024

class CSCSparsityCounter : public QStructAssembler {
public:
    CSCSparsityCounter(int m, int n);
    CSCSparsityCounter(int n);
    ~CSCSparsityCounter();

    void add(int i, int j);
    void get_counts(int* counts);

private:
    struct Entry { // Entry
        int row;
        Entry* next;
    };

    struct EntryPool {
        EntryPool(EntryPool* next) : used(0), next(next) {}
        unsigned int used;           // Number of entries used
        EntryPool* next;             // Pointer for stack
        Entry data[ENTRY_POOL_LEN];  // Pool of entry slots
    };

    int m;                   // Number of rows
    int n;                   // Number of columns
    std::vector<Entry*> col; // Column block lists
    EntryPool* pool;         // Current pool block

    CSCSparsityCounter(const CSCSparsityCounter&);
    CSCSparsityCounter& operator=(const CSCSparsityCounter&);

    Entry* new_entry();
};

/*@T
 * \subsection{Building the index structure}
 *
 * After accumulating the column counts, we want to actually build the
 * index structure.  We do this with a second structural assembly pass,
 * because (a) it keeps the memory overhead down, and (b) it's easy.
 * If this phase ever becomes the bottleneck, there are probably other
 * things we could do that would avoid the O(col_nnz^2) search/insert
 * overhead.
 *
 * The [[CSCIndexBuilder]] object is used to track this second stage.
 * I assume the [[jc]] array is already set up when the object is
 * constructed.
 * 
 *@c*/
class CSCIndexBuilder : public QStructAssembler {
public:
    CSCIndexBuilder(int m, int n, int* jc, int* ir);
    CSCIndexBuilder(int n, int* jc, int* ir);
    void add(int i, int j);
private:
    int m;    // Number of rows
    int n;    // Number of columns
    int* jc;  // Column pointers into ir
    int* ir;  // Row entries
};

void build_csc_matrix(Mesh* mesh, int n, 
                      std::vector<int>& jc, std::vector<int>& ir, 
                      int reduced);
CSCMatrix* build_csc_matrix(Mesh* mesh, int reduced);

#endif /* CSCINDEXER_H */
