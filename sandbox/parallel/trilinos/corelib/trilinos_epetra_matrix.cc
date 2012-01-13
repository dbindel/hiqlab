/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "trilinos_epetra_matrix.h"

Epetra_CrsMatrix_Complex::Epetra_CrsMatrix_Complex(Epetra_DataAccess CV,
                     const Epetra_Map& RowMap,
                     const int* NumEntriesPerRow, int is_real, bool StaticProfile)
 : Epetra_CrsMatrix(CV, RowMap, NumEntriesPerRow, StaticProfile),
   is_real(is_real), Az(NULL)
{
    if (!is_real)
        Az = new Epetra_CrsMatrix(CV, RowMap, NumEntriesPerRow, StaticProfile);
}

Epetra_CrsMatrix_Complex::Epetra_CrsMatrix_Complex(
                     const Epetra_CrsMatrix_Complex& Matrix)
 : Epetra_CrsMatrix(Matrix), is_real(Matrix.is_real_matrix()), Az(NULL)
{
    if (!is_real)
        Az = new Epetra_CrsMatrix( Matrix.view_Az() );
}

Epetra_CrsMatrix_Complex::Epetra_CrsMatrix_Complex(const Epetra_CrsMatrix& xAx,
                                                   const Epetra_CrsMatrix& xAz)
 : Epetra_CrsMatrix(xAx), is_real(0)
{
    Az = new Epetra_CrsMatrix(xAz);
}

Epetra_CrsMatrix_Complex::~Epetra_CrsMatrix_Complex()
{
    if (Az)
         delete Az;
}

int Epetra_CrsMatrix_Complex::Multiply(bool TransA,
                          const Epetra_MultiVector& X,
                          const Epetra_MultiVector& Xz,
                                Epetra_MultiVector& Y,
                                Epetra_MultiVector& Yz) const
{
    int ierr = 0;
    Epetra_MultiVector Vt(X.Map(),X.NumVectors());

    // -- Compute real part of Y
    ierr += Epetra_CrsMatrix::Multiply(TransA, X, Y);
    ierr += Az->Multiply(TransA, Xz, Vt);
    Y.Update(-1.0, Vt, 1.0);

    // -- Compute imag part of Y
    ierr += Az->Multiply(TransA, X, Yz);
    ierr += Epetra_CrsMatrix::Multiply(TransA, Xz, Vt);
    Yz.Update( 1.0, Vt, 1.0);

    ierr = (ierr==0) ? 0 : 1;

    return ierr;
}

int Epetra_CrsMatrix_Complex::Multiply(bool TransA,
                          const Epetra_MultiVector_Complex& X,
                                Epetra_MultiVector_Complex& Y) const
{
    return Multiply(TransA, X, X.view_Vz(), Y, *(Y.get_Vz()) );
}

int Epetra_CrsMatrix_Complex::Multiply(bool TransA,
                          const Epetra_Vector& X, const Epetra_Vector& Xz,
                                Epetra_Vector& Y,       Epetra_Vector& Yz) const
{
    int ierr = 0;
    Epetra_Vector Vt(X.Map());

    // -- Compute real part of Y
    ierr += Epetra_CrsMatrix::Multiply(TransA, X, Y);
    ierr += Az->Multiply(TransA, Xz, Vt);
    Y.Update(-1.0, Vt, 1.0);

    // -- Compute imag part of Y
    ierr += Az->Multiply(TransA, X, Yz);
    ierr += Epetra_CrsMatrix::Multiply(TransA, Xz, Vt);
    Yz.Update( 1.0, Vt, 1.0);

    ierr = (ierr==0) ? 0 : 1;

    return ierr;
}

int Epetra_CrsMatrix_Complex::Multiply(bool TransA,
                          const Epetra_Vector_Complex& X,
                                Epetra_Vector_Complex& Y) const
{
    return Multiply(TransA, X, X.view_Vz(), Y, *(Y.get_Vz()) );
}
