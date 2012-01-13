#ifndef _QTRILINOS_EPETRA_H
#define _QTRILINOS_EPETRA_H

#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "trilinos_epetra_vector.h"

// -- Display distributed objects in MATLAB format
int CrsMatrix2MATLAB(Epetra_CrsMatrix* Ap );
int Vector2MATLAB(Epetra_Vector* vp);
int MultiVector2MATLAB(Epetra_MultiVector* vp, int index);

// -- Helper functions for distributed matrices
Epetra_CrsMatrix* qRow2CrsMatrix(Epetra_RowMatrix* ERM);
Epetra_RowMatrix* qCrs2RowMatrix(Epetra_CrsMatrix* ECM);

// -- Helper functions for distributed vectors
Epetra_Vector* qMultiVector2Vector(Epetra_MultiVector* EMV);
Epetra_MultiVector* qVector2MultiVector(Epetra_Vector* EV);

// -- Constructor for Epetra_Vectors
Epetra_Vector* qVectorCreate(int n);
Epetra_MultiVector* qMultiVectorCreate(int m, int n);
Epetra_MultiVector_Complex* qMultiVectorComplexCreate(int m, int n);
Epetra_Vector_Complex* qVectorComplexCreate(int n);

// -- Transformation between Multi and Single
Epetra_Vector_Complex*      qMulti2SingleComplexVector(Epetra_MultiVector_Complex* v);
Epetra_MultiVector_Complex* qSingle2MultiComplexVector(Epetra_Vector_Complex* v);

#endif /* _QTRILINOS_EPETRA_H */
