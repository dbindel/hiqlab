/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */
#include "AnasaziEpetraMultiVectorComplexAdapter.hpp"
#include "trilinos_epetra_vector.h"
#include <iostream>

/** @file AnasaziEpetraMultiVectorComplexAdapter.cpp
 *
 *  Implementations of Anasazi multi-vector and operator classes using
 *  Epetra_MultiVector_Complex class
 */

namespace Anasazi {

  ///////////////////////////////////////////////////////////////////////////////
  //
  //--------Anasazi::EpetraMultiVecComplex Implementation-------------------------------
  //
  ///////////////////////////////////////////////////////////////////////////////

  // Construction/Destruction

  EpetraMultiVecComplex::EpetraMultiVecComplex(const Epetra_BlockMap& Map, dcomplex* array,
                                 const int numvecs, const int stride)
    : Epetra_MultiVector_Complex(Map, array, stride, numvecs)
  {
  }


  EpetraMultiVecComplex::EpetraMultiVecComplex(const Epetra_BlockMap& Map, const int numvecs)
    : Epetra_MultiVector_Complex(Map, numvecs)
  {
  }


  EpetraMultiVecComplex::EpetraMultiVecComplex(Epetra_DataAccess CV, const Epetra_MultiVector_Complex& P_vec,
                                 const std::vector<int>& index )
    : Epetra_MultiVector_Complex(CV, P_vec, &(const_cast<std::vector<int> &>(index))[0], index.size())
  {
  }


  EpetraMultiVecComplex::EpetraMultiVecComplex(const Epetra_MultiVector_Complex& P_vec)
    : Epetra_MultiVector_Complex(P_vec)
  {
  }


  //
  //  member functions inherited from Anasazi::MultiVec
  //
  //
  //  Simulating a virtual copy constructor. If we could rely on the co-variance
  //  of virtual functions, we could return a pointer to EpetraMultiVec
  //  (the derived type) instead of a pointer to the pure virtual base class.
  //

  MultiVec<dcomplex>* EpetraMultiVecComplex::Clone ( const int numvecs ) const
  {
    EpetraMultiVecComplex * ptr_apv = new EpetraMultiVecComplex(Map(), numvecs);
    return ptr_apv; // safe upcast.
  }
  //
  //  the following is a virtual copy constructor returning
  //  a pointer to the pure virtual class. vector values are
  //  copied.
  //

  MultiVec<dcomplex>* EpetraMultiVecComplex::CloneCopy() const
  {
    EpetraMultiVecComplex *ptr_apv = new EpetraMultiVecComplex(*this);
    return ptr_apv; // safe upcast
  }


  MultiVec<dcomplex>* EpetraMultiVecComplex::CloneCopy ( const std::vector<int>& index ) const
  {
    EpetraMultiVecComplex * ptr_apv = new EpetraMultiVecComplex(Copy, *this, index);
    return ptr_apv; // safe upcast.
  }


  MultiVec<dcomplex>* EpetraMultiVecComplex::CloneView ( const std::vector<int>& index )
  {
    EpetraMultiVecComplex * ptr_apv = new EpetraMultiVecComplex(View, *this, index);
    return ptr_apv; // safe upcast.
  }


  void EpetraMultiVecComplex::SetBlock( const MultiVec<dcomplex>& A, const std::vector<int>& index )
  {
    EpetraMultiVecComplex temp_vec(View, *this, index);

    int numvecs = index.size();
    if ( A.GetNumberVecs() != numvecs ) {
      std::vector<int> index2( numvecs );
      for(int i=0; i<numvecs; i++)
        index2[i] = i;
      EpetraMultiVecComplex *tmp_vec = dynamic_cast<EpetraMultiVecComplex *>(&const_cast<MultiVec<dcomplex> &>(A));
      assert(tmp_vec!=NULL);
      EpetraMultiVecComplex A_vec(View, *tmp_vec, index2);
      // FIXME: Might need to change scalars
      temp_vec.MvAddMv( dcomplex(1.0,0.0), A_vec, dcomplex(0.0,0.0), A_vec );
    }
    else {
      // FIXME: Might need to change scalars
      temp_vec.MvAddMv( dcomplex(1.0,0.0), A, dcomplex(0.0,0.0), A );
    }
  }

  //-------------------------------------------------------------
  //
  // *this <- alpha * A * B + beta * (*this)
  //
  //-------------------------------------------------------------

  void EpetraMultiVecComplex::MvTimesMatAddMv ( const dcomplex alpha, const MultiVec<dcomplex>& A,
                                         const Teuchos::SerialDenseMatrix<int,dcomplex>& B, const dcomplex beta )
  {
    Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
    Epetra_MultiVector_Complex B_Pvec(LocalMap, B.values(), B.stride(), B.numCols());
//    Epetra_MultiVector B_Pvec(Copy, LocalMap, B.values(), B.stride(), B.numCols());

    EpetraMultiVecComplex *A_vec = dynamic_cast<EpetraMultiVecComplex *>(&const_cast<MultiVec<dcomplex> &>(A));
    assert(A_vec!=NULL);

    int ret = Multiply( 'N', 'N', alpha, *A_vec, B_Pvec, beta );
    assert( ret == 0 );
  }

  //-------------------------------------------------------------
  //
  // *this <- alpha * A + beta * B
  //
  //-------------------------------------------------------------

  void EpetraMultiVecComplex::MvAddMv ( const dcomplex alpha , const MultiVec<dcomplex>& A,
                                 const dcomplex beta, const MultiVec<dcomplex>& B)
  {
    EpetraMultiVecComplex *A_vec = dynamic_cast<EpetraMultiVecComplex *>(&const_cast<MultiVec<dcomplex> &>(A));
    assert(A_vec!=NULL);
    EpetraMultiVecComplex *B_vec = dynamic_cast<EpetraMultiVecComplex *>(&const_cast<MultiVec<dcomplex> &>(B));
    assert(B_vec!=NULL);

    // FIXME: Might need to change scalars
    int ret = Update( alpha, *A_vec, beta, *B_vec, dcomplex(0.0,0.0) );
    assert( ret == 0 );
  }

  //-------------------------------------------------------------
  //
  // dense B <- alpha * A^T * (*this)
  //
  //-------------------------------------------------------------

  void EpetraMultiVecComplex::MvTransMv ( const dcomplex alpha, const MultiVec<dcomplex>& A,
                                   Teuchos::SerialDenseMatrix<int,dcomplex>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                                   , ConjType conj
#endif
                                   ) const
  {
    EpetraMultiVecComplex *A_vec = dynamic_cast<EpetraMultiVecComplex *>(&const_cast<MultiVec<dcomplex> &>(A));

    if (A_vec) {
      Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
      Epetra_MultiVector_Complex B_Pvec(LocalMap, B.values(), B.stride(), B.numCols());
//      Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());

      int ret = B_Pvec.Multiply( 'T', 'N', alpha, *A_vec, *this, 0.0 );
      dcomplex* v = new dcomplex[B.stride()*B.numCols()];
      B_Pvec.ExtractCopy(v, B.stride());
      Teuchos::SerialDenseMatrix<int,dcomplex> B_temp(Teuchos::View, v, B.stride(), B.stride(), B.numCols());
      B.assign(B_temp);
      assert( ret == 0 );
      // Delete temporaries
      delete[] v;
    }
  }

  //-------------------------------------------------------------
  //
  // b[i] = A[i]^T * this[i]
  //
  //-------------------------------------------------------------

  void EpetraMultiVecComplex::MvDot ( const MultiVec<dcomplex>& A, std::vector<dcomplex>* b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                               , ConjType conj
#endif
                               ) const
  {
    EpetraMultiVecComplex *A_vec = dynamic_cast<EpetraMultiVecComplex *>(&const_cast<MultiVec<dcomplex> &>(A));
    assert(A_vec!=NULL);
    if ((A_vec!=NULL) && (b!=NULL) && ( (int)b->size() >= A_vec->NumVectors() ) ) {
      int ret = this->Dot( *A_vec, &(*b)[0] );
      assert( ret == 0 );
    }
  }

  //-------------------------------------------------------------
  //
  // this[i] = alpha[i] * this[i]
  //
  //-------------------------------------------------------------
  void EpetraMultiVecComplex::MvScale ( const std::vector<dcomplex>& alpha )
  {
    // Check to make sure the vector is as long as the multivector has columns.
    int numvecs = this->NumVectors();
    assert( (int)alpha.size() == numvecs );

    int ret = 0;
    std::vector<int> tmp_index( 1, 0 );
    for (int i=0; i<numvecs; i++) {
      Epetra_MultiVector_Complex temp_vec(View, *this, &tmp_index[0], 1);
      ret = temp_vec.Scale( alpha[i] );
      assert (ret == 0);
      tmp_index[0]++;
    }
  }
} // end namespace Anasazi
