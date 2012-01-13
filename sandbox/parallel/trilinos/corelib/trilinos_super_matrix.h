#ifndef _TRILINOS_SUPER_MATRIX_H
#define _TRILINOS_SUPER_MATRIX_H

/** @file trilinos_super_matrix.h
 *
 *  Classes and functions related to Epetra
 */

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "trilinos_epetra_vector.h"
#include "trilinos_epetra_operator.h"
#include "trilinos_indexmap.h"

#include <vector>

/** Data structure for Epetra_Super_CrsMatrix
 */
class Epetra_Super_CrsMatrix {
 public:

    /** Create a new matrix
     */
    Epetra_Super_CrsMatrix(Epetra_Map& RowMap);

    virtual ~Epetra_Super_CrsMatrix();

    /** Add an Epetra_Crsmatrix */
    void InsertMatrix(Epetra_CrsMatrix* A, IndexMap* im);
    void InsertMatrix(Epetra_CrsMatrix* A, int* rind_mapping, int* cind_mapping);
    void InsertMatrix(Epetra_CrsMatrix* A, int* rind_mapping);

    /** Construct matrix */
    void ConstructMatrix();

    /** Get CrsMatrix */
    Epetra_CrsMatrix* GetCrsMatrix() { return C;};

 private:
    int nsm;
    Epetra_CrsMatrix* C;
    Epetra_Map* FullMap;
    std::vector<Epetra_CrsMatrix *> sub_mat;
    std::vector<IndexMap *> indexmap;

    /** Create and Fill matrix */
    void CreateMatrix();
    void FillMatrix();

    /** Set NNZ structure for super matrix */
    void setMyNNZ(int* MyNNZ);

    /** Compute nnz structure of newly distributed matrix from old matrix*/
    void NumMyEntriesPerRow(const Epetra_CrsMatrix& A,
                      const Epetra_Map& TargetMap,
                      int* MyNNZ_B);

    /** Construct mapping based on the function given */
    Epetra_Map s2f_rowmap(int ism);
};

#endif /* _TRILINOS_SUPER_MATRIX_H*/
