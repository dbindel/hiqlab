#ifndef _TRILINOS_SUPER_VECTOR_H
#define _TRILINOS_SUPER_VECTOR_H

/** @file trilinos_super_vector.h
 *
 *  Classes and functions related to Epetra
 */

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "trilinos_epetra_vector.h"
#include "trilinos_epetra_operator.h"
#include "trilinos_indexmap.h"

#include <vector>

/** Data structure for Epetra_Super_MultiVector
 */
class Epetra_Super_MultiVector {
 public:
	
    /** Create a new MultiVector
     */
    Epetra_Super_MultiVector(Epetra_Map& RowMap, int NumVectors);

    virtual ~Epetra_Super_MultiVector();

    /** Add an Epetra_MultiVector */
    void InsertVector(Epetra_MultiVector* A, IndexMap* im);
    void InsertVector(Epetra_MultiVector* A, int* rind_mapping, int* cind_mapping);
    void InsertVector(Epetra_MultiVector* A, int* rind_mapping);

    /** Construct matrix */
    void ConstructVector();

    /** Get MultiVector */
    Epetra_MultiVector* GetMultiVector() { return C;};


 private:
    int nsv;
    Epetra_MultiVector* C;
    Epetra_Map* FullMap;
    std::vector<Epetra_MultiVector *> sub_vec;
    std::vector<IndexMap *> indexmap;

    /** Create and Fill matrix */
    void FillVector();

    /** Construct mapping based on the function given */
    Epetra_Map s2f_rowmap(int ism);
};

/** Extract MultiVector */
void ExtractMultiVector(Epetra_MultiVector* Source, Epetra_MultiVector* B, IndexMap* im);

#endif /* _TRILINOS_SUPER_VECTOR_H*/
