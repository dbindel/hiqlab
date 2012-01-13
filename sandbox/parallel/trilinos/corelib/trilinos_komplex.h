#ifndef _TRILINOS_KOMPLEX_H
#define _TRILINOS_KOMPLEX_H

#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_vector.h"
#include "trilinos_super_matrix.h"
#include "trilinos_super_vector.h"
#include "trilinos_indexmap.h"

/** @file trilinos_komplex.h
 *
 *  Classes and functions related to Komplex form
 */

Epetra_Super_CrsMatrix* Complex2SuperCrsMatrix(Epetra_CrsMatrix_Complex* A,
                                               int form=0);
Epetra_Super_MultiVector* Complex2SuperMultiVector(Epetra_MultiVector_Complex* A,                                               int form=0);
Epetra_Super_MultiVector* Complex2SuperMultiVector(Epetra_Vector_Complex* A,                                               int form=0);
void Multi2ComplexMultiVector(Epetra_MultiVector* Source,
                              Epetra_MultiVector_Complex* A,
                                               int form=0);

#endif /* _TRILINOS_KOMPLEX_H*/
