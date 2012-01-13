#ifndef _TRILINOS_EPETRA_MATRIX_H
#define _TRILINOS_EPETRA_MATRIX_H

/** @file trilinos_epetra_matrix.h
 *
 *  Classes and functions related to Epetra
 */

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "trilinos_epetra_vector.h"
#include "trilinos_epetra_operator.h"

/** Complex data structure for Epetra_CrsMatrix
 */
class Epetra_CrsMatrix_Complex : public Epetra_CrsMatrix, public virtual Epetra_Operator_Complex{
 public:

    /** Create a new matrix
     */
    Epetra_CrsMatrix_Complex(Epetra_DataAccess CV, const Epetra_Map& RowMap,
       const int* NumEntriesPerRow, int is_real, bool StaticProfile = false);

    /** Copy constructors
     */
    Epetra_CrsMatrix_Complex(const Epetra_CrsMatrix_Complex& Matrix);
    Epetra_CrsMatrix_Complex(const Epetra_CrsMatrix& Ax,
                             const Epetra_CrsMatrix& Az);
    virtual ~Epetra_CrsMatrix_Complex();

    /** View or get imag part of matrix */
    const Epetra_CrsMatrix& view_Az() const { return *Az; };
    Epetra_CrsMatrix* get_Az()              { return Az;  };

    /** Is real matrix ? */
    int is_real_matrix() const { return is_real; };


    /** Multiply by multivector */
    int Multiply(bool TransA, const Epetra_Vector_Complex& X,
                                    Epetra_Vector_Complex& Y) const;
    int Multiply(bool TransA, const Epetra_Vector& Xx, const Epetra_Vector& Xz,
                                    Epetra_Vector& Yx,       Epetra_Vector& Yz) const;

    int Multiply(bool TransA, const Epetra_MultiVector_Complex& X,
                                    Epetra_MultiVector_Complex& Y) const;
    int Multiply(bool TransA, const Epetra_MultiVector& Xx,
                              const Epetra_MultiVector& Xz,
                                    Epetra_MultiVector& Yx,
                                    Epetra_MultiVector& Yz) const;

    /** Apply to multivector */
    int Apply(const Epetra_MultiVector_Complex& X,
                    Epetra_MultiVector_Complex& Y) const
       {return(Epetra_CrsMatrix_Complex::Multiply(Epetra_CrsMatrix::UseTranspose(), X, Y));}


    /** Methods required by inheritance of Epetra_Operator_Complex
     *  (Inheritance done by 'virtual' comment since there are classes inherited
     *   twice.)
     */

    /** Multiply by the transpose
     * (This method has no effect and returns -1 as error code.)
     */
    int SetUseTranspose(bool UseTranspose){return(-1);};

    /** Apply inverse to a Epetra_MultiVector_Complex X
     *  (This method has no effect and returns -1 as error code.)
     */
    int ApplyInverse(const Epetra_MultiVector_Complex& X,
                           Epetra_MultiVector_Complex& Y) const{return 0;};

    /** Returns the infinity norm of the global matrix.
     *warning This method must not be called unless HasNormInf() returns true.
     */
    double NormInf() const {return(0.0);};

    /** Returns a character string describing the operator */
    const char * Label() const {return(Epetra_CrsMatrix::Label());};

    /** Returns the current UseTranspose setting. */
    bool UseTranspose() const {return(false);};

    /** Returns true if this object can provide an approximate Inf-norm, false otherwise. */
    bool HasNormInf() const{return(false);};

    /** Returns a pointer to the Epetra_Comm communicator associated
     * with this operator.
     */
    const Epetra_Comm & Comm() const{return(Epetra_CrsMatrix::Comm());};
    //const Epetra_Comm & Comm() const{return(ELP->GetOperator()->Comm());};

    /** Returns the Epetra_BlockMap object associated with the domain of
     * this matrix operator.
     */
    const Epetra_Map & OperatorDomainMap() const
    {
      if (UseTranspose()) return(Epetra_CrsMatrix::RangeMap());
      else return(Epetra_CrsMatrix::DomainMap());
      //if (UseTranspose()) return(ELP->GetOperator()->OperatorRangeMap());
      //else return(ELP->GetOperator()->OperatorDomainMap());
    }

    /** Returns the Epetra_BlockMap object associated with the range of
     * this matrix operator.
     */
    const Epetra_Map & OperatorRangeMap() const
    {
      if (UseTranspose()) return(Epetra_CrsMatrix::DomainMap());
      else return(Epetra_CrsMatrix::RangeMap());
      //if (UseTranspose()) return(ELP->GetOperator()->OperatorDomainMap());
      //else return(ELP->GetOperator()->OperatorRangeMap());
    }

 private:
    int is_real;
    Epetra_CrsMatrix* Az;   /* Imag part     */

};

#endif /* _TRILINOS_EPETRA_MATRIX_H */
