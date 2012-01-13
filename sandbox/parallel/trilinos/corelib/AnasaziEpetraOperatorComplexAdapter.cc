/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include "AnasaziEpetraOperatorComplexAdapter.hpp"
#include "AnasaziEpetraMultiVectorComplexAdapter.hpp"
#include <iostream>

/*! \file AnasaziEpetraOperatorComplexAdapter.cpp
 *   \brief Implementations of Anasazi multi-vector and operator classes using Epetra_MultiVector_Complex and Epetra_Operator_Complex classes
 *   */

namespace Anasazi {

  ///////////////////////////////////////////////////////////////////////////////
  //
  //--------Anasazi::EpetraOpComplex Implementation----------------------------
  //
  ///////////////////////////////////////////////////////////////////////////////

  //
  // AnasaziOperator constructors
  //
  EpetraOpComplex::EpetraOpComplex(const Teuchos::RefCountPtr<Epetra_Operator_Complex> &Op)
    : Epetra_Op_Complex(Op)
  {
  }

  EpetraOpComplex::~EpetraOpComplex()
  {
  }
  //
  // AnasaziOperator applications
  //
  void EpetraOpComplex::Apply ( const MultiVec<dcomplex>& X,
                               MultiVec<dcomplex>& Y ) const
  {
    //
    // This standard operator computes Y = A*X
    //
    MultiVec<dcomplex> & temp_X = const_cast<MultiVec<dcomplex> &>(X);
    Epetra_MultiVector_Complex* vec_X = dynamic_cast<Epetra_MultiVector_Complex* >(&temp_X);
    Epetra_MultiVector_Complex* vec_Y = dynamic_cast<Epetra_MultiVector_Complex* >(&Y);
    assert( vec_X!=NULL && vec_Y!=NULL );

    int info = Epetra_Op_Complex->Apply( *vec_X, *vec_Y );
    TEST_FOR_EXCEPTION( info != 0, OperatorError,
                        "Anasazi::EpetraOpComplex::Apply(): Error returned from Epetra_Operator_Complex::Apply()" );
  }
} // end namespace Anasazi
