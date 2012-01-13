/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#ifndef ANASAZI_EPETRA_OPERATOR_COMPLEX_ADAPTER_HPP
#define ANASAZI_EPETRA_OPERATOR_COMPLEX_ADAPTER_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziOperator.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "trilinos_epetra_operator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "qcomplex.h"

/** @file AnasaziEpetraOperatorComplexAdapter.hpp
 *
 *  Declarations of Anasazi multi-vector and operator classes using
 *  Epetra_MultiVector_Complex class
 */

namespace Anasazi {
  
  //-------------------------------------------------------------
  
  ///////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziEpetraOpComplex---------------
  //
  ///////////////////////////////////////////////////////////////
  
  /*! 
    \brief Basic adapter class for Anasazi::Operator that uses Epetra_Operator_Complex.
  */
  class EpetraOpComplex : public virtual Operator<dcomplex> {
  public:
    //! @name Constructor/Destructor
    //@{ 
    
    //! Basic constructor.  Accepts reference-counted pointer to an Epetra_Operator_Complex.
    EpetraOpComplex(const Teuchos::RefCountPtr<Epetra_Operator_Complex> &Op );
    
    //! Destructor
    ~EpetraOpComplex();
    //@}
    
    //! @name Operator application method
    //@{ 
    
    /*! \brief This method takes the Anasazi::MultiVec \c X and
      applies the operator to it resulting in the Anasazi::MultiVec \c Y.
    */
    void Apply ( const MultiVec<dcomplex>& X, MultiVec<dcomplex>& Y ) const;
    //@} 
    
  private:
    Teuchos::RefCountPtr<Epetra_Operator_Complex> Epetra_Op_Complex;
  };
  //-------------------------------------------------------------

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::OperatorTraits for Epetra::Operator_Complex.
  //
  ////////////////////////////////////////////////////////////////////

  /*! 
    \brief Template specialization of Anasazi::OperatorTraits class using the Epetra_Operator_Complex virtual base class and 
    Epetra_MultiVector_Complex class.

    This interface will ensure that any Epetra_Operator_Complex and Epetra_MultiVector_Complex will be accepted by the Anasazi
    templated solvers.
  */

  template <> 
  class OperatorTraits < dcomplex, Epetra_MultiVector_Complex, Epetra_Operator_Complex >
  {
  public:
    
    /*! \brief This method takes the Epetra_MultiVector_Complex \c x and
      applies the Epetra_Operator_Complex \c Op to it resulting in the Epetra_MultiVector_Complex \c y.
    */    
    static void Apply ( const Epetra_Operator_Complex& Op, 
			      const Epetra_MultiVector_Complex& x, 
			      Epetra_MultiVector_Complex& y )
    { 
      TEST_FOR_EXCEPTION( Op.Apply( x, y ) != 0, OperatorError, "Error in Epetra_Operator_Complex::Apply()!" );
    }
    
  };
  
} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_EPETRA_OPERATOR_COMPLEX_ADAPTER_HPP
