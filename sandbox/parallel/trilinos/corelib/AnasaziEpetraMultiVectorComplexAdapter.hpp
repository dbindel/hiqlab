/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#ifndef ANASAZI_EPETRA_MULTIVECTOR_COMPLEX_ADAPTER_HPP
#define ANASAZI_EPETRA_MULTIVECTOR_COMPLEX_ADAPTER_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziMultiVec.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "trilinos_epetra_vector.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "qcomplex.h"
#include <iostream>

/** @file AnasaziEpetraMultiVectorComplexAdapter.hpp
 * 
 *  Declarations of Anasazi multi-vector and operator classes using 
 *  Epetra_MultiVector_Complex class
 */

namespace Anasazi {
  
  ///////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziEpetraMultiVecComplex----------
  //
  ///////////////////////////////////////////////////////////////
  
  /*! 
    \brief Basic adapter class for Anasazi::MultiVec that uses Epetra_MultiVector_Complex.
  */
  class EpetraMultiVecComplex : public MultiVec<dcomplex>, public Epetra_MultiVector_Complex {
  public:
    //! @name Constructors/Destructors
    //@{ 

    //! Basic EpetraMultiVecComplex constructor.
    /*! @param Map [in] An Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
      @param numvecs [in] Number of vectors in multi-vector.

      \returns Pointer to an EpetraMultiVecComplex
    */
    EpetraMultiVecComplex(const Epetra_BlockMap& Map, const int numvecs);

    //! Copy constructor.
    EpetraMultiVecComplex(const Epetra_MultiVector_Complex & P_vec);
    
    //! Create multi-vector with values from two dimensional array.
    /*! @param Map [in] An Epetra_LocalMap, Epetra_Map or Epetra_BlockMap
      @param array [in] Pointer to an array of complex double precision numbers.  The first vector starts at \c array, the
      second at \c array+stride, and so on.  This array is copied.
      @param numvecs [in] Number of vectors in the multi-vector.
      @param stride [in] The stride between vectors in memory of \c array.

      \returns Pointer to an EpetraMultiVecComplex
    */
    EpetraMultiVecComplex(const Epetra_BlockMap& Map, dcomplex* array, const int numvecs, const int stride=0);

    //! Create multi-vector from list of vectors in an existing EpetraMultiVecComplex.
    /*! @param CV [in] Enumerated type set to Copy or View.
      @param P_vec [in] An existing fully constructed Epetra_MultiVector_Complex.
      @param index [in] A integer vector containing the indices of the vectors to copy out of \c P_vec.

      \returns Pointer to an EpetraMultiVecComplex
    */
    EpetraMultiVecComplex(Epetra_DataAccess CV, const Epetra_MultiVector_Complex& P_vec, const std::vector<int>& index);

    //! Destructor
    virtual ~EpetraMultiVecComplex() {};

    //@}

    //! @name Creation methods
    //@{ 

    /*! \brief Creates a new empty EpetraMultiVecComplex containing \c numvecs columns.
      
    \returns Pointer to an EpetraMultiVecComplex
    */
    MultiVec<dcomplex> * Clone ( const int numvecs ) const;

    /*! \brief Creates a new EpetraMultiVecComplex and copies contents of \c *this into
      the new vector (deep copy).
      
      \returns Pointer to an EpetraMultiVecComplex
    */	
    MultiVec<dcomplex> * CloneCopy () const;

    /*! \brief Creates a new EpetraMultiVecComplex and copies the selected contents of \c *this 
      into the new vector (deep copy).  
      
      The copied vectors from \c *this are indicated by the \c index.size() indices in \c index.
      
      \returns Pointer to an EpetraMultiVecComplex
    */
    MultiVec<dcomplex> * CloneCopy ( const std::vector<int>& index ) const;
    
    /*! \brief Creates a new EpetraMultiVecComplex that shares the selected contents of \c *this.
      
    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.
    
    \returns Pointer to an EpetraMultiVecComplex
    */
    MultiVec<dcomplex> * CloneView ( const std::vector<int>& index );

    //@}

    //! @name Attribute methods	
    //@{ 

    //! Obtain the vector length of *this.
    int GetNumberVecs () const { return NumVectors(); }

    //! Obtain the number of vectors in *this.
    int GetVecLength () const { return GlobalLength(); }

    //@}

    //! @name Update methods
    //@{ 
    /*! \brief Update \c *this with \f$\alpha AB + \beta (*this)\f$.
     */
    void MvTimesMatAddMv ( const dcomplex alpha, const MultiVec<dcomplex>& A, 
			   const Teuchos::SerialDenseMatrix<int,dcomplex>& B, const dcomplex beta );

    /*! \brief Replace \c *this with \f$\alpha A + \beta B\f$.
     */
    void MvAddMv ( const dcomplex alpha, const MultiVec<dcomplex>& A, const dcomplex beta,
		   const MultiVec<dcomplex>& B);

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$\alpha A^T(*this)\f$.
    */
    void MvTransMv ( const dcomplex alpha, const MultiVec<dcomplex>& A, Teuchos::SerialDenseMatrix<int,dcomplex>& B 
#ifdef HAVE_ANASAZI_EXPERIMENTAL
		     , ConjType conj = Anasazi::CONJ
#endif
		     ) const;
  
    /*! \brief Compute a vector \c b where the components are the individual dot-products, i.e. \f$ b[i] = A[i]^H(this[i])\f$ where \c A[i] is the i-th column of \c A.
	*/
    void MvDot ( const MultiVec<dcomplex>& A, std::vector<dcomplex>* b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
		 , ConjType conj = Anasazi::CONJ
#endif
		 ) const;

    /*! \brief Scale each element of the vectors in \c *this with \c alpha.
     */
    void MvScale ( const dcomplex alpha ) { int ret = this->Scale( alpha ); assert(ret == 0); }
    
    /*! \brief Scale each element of the \c i-th vector in \c *this with \c alpha[i].
     */
    void MvScale ( const std::vector<dcomplex>& alpha );

    //@}
    //! @name Norm method
    //@{ 
    
    /*! \brief Compute the 2-norm of each individual vector of \c *this.  
      Upon return, \c normvec[i] holds the 2-norm of the \c i-th vector of \c *this
    */
    void MvNorm ( std::vector<double>* normvec ) const {
      if ((normvec!=NULL) && ((int)normvec->size() >= GetNumberVecs()) ) {
	int ret = Norm2(&(*normvec)[0]);
	assert( ret == 0 );
      }
    };
    //@}

    //! @name Initialization methods
    //@{ 
    /*! \brief Copy the vectors in \c A to a set of vectors in \c *this.  

    The \c numvecs vectors in \c A are copied to a subset of vectors in \c *this
    indicated by the indices given in \c index.
    */
    void SetBlock ( const MultiVec<dcomplex>& A, const std::vector<int>& index );

    /*! \brief Fill the vectors in \c *this with random numbers.
     */
    void MvRandom() { int ret = Random(); assert( ret == 0 ); };

    /*! \brief Replace each element of the vectors in \c *this with \c alpha.
     */
    void MvInit ( const dcomplex alpha ) { int ret = PutScalar( alpha ); assert( ret == 0 ); };

    //@}
    //! @name Print method
    //@{ 
    /*! \brief Print \c *this EpetraMultiVec.
     */
    void MvPrint( ostream& os ) const 
    { 
	os << "Real part of EpetraMultiVecComplex\n";
	os << *this << endl; 
	os << "Imag part of EpetraMultiVecComplex\n";
	os <<  this->view_Vz() << endl;
    };
    //@}

  private:
  };
  //-------------------------------------------------------------
  
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::MultiVecTraits for Epetra::MultiVector_Complex.
  //
  ////////////////////////////////////////////////////////////////////

  /*! 
    \brief Template specialization of Anasazi::MultiVecTraits class using the Epetra_MultiVector_Complex class.

    This interface will ensure that any Epetra_MultiVector_Complex will be accepted by the Anasazi
    templated solvers.  
  */

  template<>
  class MultiVecTraits<dcomplex, Epetra_MultiVector_Complex>
  {
  public:

    //! @name Creation methods
    //@{ 

    /*! \brief Creates a new empty Epetra_MultiVector_Complex containing \c numvecs columns.
      
    \return Reference-counted pointer to the new Epetra_MultiVector_Complex.
    */
    static Teuchos::RefCountPtr<Epetra_MultiVector_Complex> Clone( const Epetra_MultiVector_Complex& mv, const int numvecs )
    { return Teuchos::rcp( new Epetra_MultiVector_Complex(mv.Map(), numvecs) ); }

    /*! \brief Creates a new Epetra_MultiVector_Complex and copies contents of \c mv into the new vector (deep copy).
      
      \return Reference-counted pointer to the new Epetra_MultiVector_Complex.
    */
    static Teuchos::RefCountPtr<Epetra_MultiVector_Complex> CloneCopy( const Epetra_MultiVector_Complex& mv )
    { return Teuchos::rcp( new Epetra_MultiVector_Complex( mv ) ); }

    /*! \brief Creates a new Epetra_MultiVector_Complex and copies the selected contents of \c mv into the new vector (deep copy).  

      The copied vectors from \c mv are indicated by the \c indeX.size() indices in \c index.      
      \return Reference-counted pointer to the new Epetra_MultiVector_Complex.
    */
    static Teuchos::RefCountPtr<Epetra_MultiVector_Complex> CloneCopy( const Epetra_MultiVector_Complex& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector_Complex(::Copy, mv, &tmp_index[0], index.size()) ); 
    }

    /*! \brief Creates a new Epetra_MultiVector that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new Epetra_MultiVector.
    */      
    static Teuchos::RefCountPtr<Epetra_MultiVector_Complex> CloneView( Epetra_MultiVector_Complex& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector_Complex(::View, mv, &tmp_index[0], index.size()) ); 
    }

    /*! \brief Creates a new const Epetra_MultiVector that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const Epetra_MultiVector.
    */      
    static Teuchos::RefCountPtr<const Epetra_MultiVector_Complex> CloneView( const Epetra_MultiVector_Complex& mv, const std::vector<int>& index )
    { 
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector_Complex(::View, mv, &tmp_index[0], index.size()) ); 
    }

    //@}

    //! @name Attribute methods
    //@{ 

    //! Obtain the vector length of \c mv.
    static int GetVecLength( const Epetra_MultiVector_Complex& mv )
    { return mv.GlobalLength(); }

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const Epetra_MultiVector_Complex& mv )
    { return mv.NumVectors(); }
    //@}

    //! @name Update methods
    //@{ 

    /*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv( const dcomplex alpha, const Epetra_MultiVector_Complex& A, 
				 const Teuchos::SerialDenseMatrix<int,dcomplex>& B, 
				 const dcomplex beta, Epetra_MultiVector_Complex& mv )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector_Complex B_Pvec(LocalMap, B.values(), B.stride(), B.numCols());
//      Epetra_MultiVector_Complex B_Pvec(::Copy, LocalMap, B.values(), B.stride(), B.numCols());

      int ret = mv.Multiply( 'N', 'N', alpha, A, B_Pvec, beta );
      assert( ret == 0 );   
    }

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( const dcomplex alpha, const Epetra_MultiVector_Complex& A, const dcomplex beta, const Epetra_MultiVector_Complex& B, Epetra_MultiVector_Complex& mv )
    { 
      // FIXME: Might have to change 0.0 to ScalarType
      int ret = mv.Update( alpha, A, beta, B, dcomplex(0.0,0.0) );
      assert( ret == 0 );
    }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
    */
    static void MvTransMv( const dcomplex alpha, const Epetra_MultiVector_Complex& A, const Epetra_MultiVector_Complex& mv, Teuchos::SerialDenseMatrix<int,dcomplex>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
			   , ConjType conj = Anasazi::CONJ
#endif
			   )
    { 
      // FIXME: Currently B_Pvec copies B, would like to point
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector_Complex B_Pvec(LocalMap, B.values(), B.stride(), B.numCols());
//      Epetra_MultiVector B_Pvec(::View, LocalMap, B.values(), B.stride(), B.numCols());
      
      // FIXME: Might have to change 0.0 to ScalarType
      int ret = B_Pvec.Multiply( 'T', 'N', alpha, A, mv, dcomplex(0.0,0.0) );
      dcomplex* v = new dcomplex[B.stride()*B.numCols()];
      B_Pvec.ExtractCopy(v, B.stride());
      Teuchos::SerialDenseMatrix<int,dcomplex> B_temp(Teuchos::View, v, B.stride(), B.stride(), B.numCols());
      B.assign(B_temp);
      assert( ret == 0 );
      // Delete temporaries
      delete[] v;
    }
    
    /*! \brief Compute a vector \c b where the components are the individual dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^Tmv[i]\f$.
     */
    static void MvDot( const Epetra_MultiVector_Complex& mv, const Epetra_MultiVector_Complex& A, std::vector<dcomplex>* b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
		       , ConjType conj = Anasazi::CONJ
#endif
		       )
    {
      int ret = mv.Dot( A, &(*b)[0] );
      assert( ret == 0 );
    }

    //@}
    //! @name Norm method
    //@{ 

    /*! \brief Compute the 2-norm of each individual vector of \c mv.  
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
    */
    static void MvNorm( const Epetra_MultiVector_Complex& mv, std::vector<double>* normvec )
    { 
      int ret = mv.Norm2(&(*normvec)[0]);
      assert( ret == 0 );
    }

    //@}

    //! @name Initialization methods
    //@{ 
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
     */
    static void SetBlock( const Epetra_MultiVector_Complex& A, const std::vector<int>& index, Epetra_MultiVector_Complex& mv )
    { 
      // Extract the "numvecs" columns of mv indicated by the index vector.
      int numvecs = index.size();
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      Epetra_MultiVector_Complex temp_vec(::View, mv, &tmp_index[0], numvecs);

      if ( A.NumVectors() != numvecs ) {
        std::vector<int> index2( numvecs );
        for(int i=0; i<numvecs; i++)
	  index2[i] = i;
        Epetra_MultiVector_Complex A_vec(::View, A, &index2[0], numvecs);      
        // FIXME: Might have to change 0.0 to ScalarType
        int ret = temp_vec.Update( 1.0, A_vec, 0.0, A_vec, 0.0 );
	assert( ret == 0 );
      }
      else {
        // FIXME: Might have to change 0.0 to ScalarType
        int ret = temp_vec.Update( 1.0, A, 0.0, A, 0.0 );
	assert( ret == 0 );
      }
    }

    /*! \brief Scale each element of the vectors in \c *this with \c alpha.
     */
    void MvScale ( const dcomplex alpha, Epetra_MultiVector_Complex& mv ) 
    { int ret = mv.Scale( alpha ); 
      assert( ret == 0 );
    }
    
    /*! \brief Scale each element of the \c i-th vector in \c *this with \c alpha[i].
     */
    void MvScale ( const std::vector<dcomplex>& alpha, Epetra_MultiVector_Complex& mv )
    { 
      // Check to make sure the vector is as long as the multivector has columns.
      int numvecs = mv.NumVectors();
      assert( (int)alpha.size() == numvecs );

      int ret = 0;
      std::vector<int> tmp_index( 1, 0 );
      for (int i=0; i<numvecs; i++) {
	Epetra_MultiVector_Complex temp_vec(::View, mv, &tmp_index[0], 1);
        ret = temp_vec.Scale( alpha[i] );
	assert (ret == 0);
	tmp_index[0]++;
      }
    }

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( Epetra_MultiVector_Complex& mv )
    { int ret = mv.Random(); assert( ret == 0 ); }

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( Epetra_MultiVector_Complex& mv, dcomplex alpha = Teuchos::ScalarTraits<dcomplex>::zero() )
    { int ret = mv.PutScalar(alpha); assert( ret == 0 ); }

    //@}

    //! @name Print method
    //@{ 

    /*! \brief Print the \c mv multi-vector to the \c os output stream.
     */
    static void MvPrint( const Epetra_MultiVector_Complex& mv, ostream& os )
    { os << mv << endl; }

    //@}
  };        

} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_EPETRA_MULTIVECTOR_COMPLEX_ADAPTER_HPP
