// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// This test is for BlockKrylovSchur solving a standard (Ax=xl) complex Hermitian
// eigenvalue problem.
//
// The matrix used is from MatrixMarket:
// Name: MHD1280B: Alfven Spectra in Magnetohydrodynamics
// Source: Source: A. Booten, M.N. Kooper, H.A. van der Vorst, S. Poedts and J.P. Goedbloed University of Utrecht, the Netherlands
// Discipline: Plasma physics
// URL: http://math.nist.gov/MatrixMarket/data/NEP/mhd/mhd1280b.html
// Size: 1280 x 1280
// NNZ: 22778 entries

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#endif

// I/O for Harwell-Boeing files
#ifdef HAVE_ANASAZI_TRIUTILS
#include "iohb.h"
#endif

// templated multivector and sparse matrix classes
#include "AnasaziEpetraMultiVectorComplexAdapter.hpp"
#include "AnasaziEpetraOperatorComplexAdapter.hpp"
#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_DataAccess.h"
#include "trilinos_epetra_operator.h"
#include "Amesos_Operator_Complex.h"

#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"
#include "qcomplex.h"

int main(int argc, char *argv[]) 
{
  int info = 0;
  int MyPID = 0;
  bool boolret;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
  Teuchos::RefCountPtr<Epetra_MpiComm> Comm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#endif


  bool testFailed;
  bool verbose = false;
  bool debug = false;
  std::string filename("mhd1280b.cua");
  std::string which("LM");
/*
  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
*/
  verbose = true;

#ifndef HAVE_ANASAZI_TRIUTILS
  cout << "This test requires Triutils. Please configure with --enable-triutils." << endl;
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  if (verbose && MyPID == 0) {
    cout << "End Result: TEST FAILED" << endl;	
  }
  return -1;
#endif

#ifdef HAVE_COMPLEX
  typedef std::complex<double> ST;
#elif HAVE_COMPLEX_H
  typedef ::complex<double> ST;
#else
  typedef double ST;
  // no complex. quit with failure.
  if (verbose && MyPID == 0) {
    cout << "Not compiled with complex support." << endl;
    cout << "End Result: TEST FAILED" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
#endif
  typedef Teuchos::ScalarTraits<ST>          SCT;
  typedef SCT::magnitudeType                  MT;
  typedef Anasazi::MultiVec<ST>               MV;
  typedef Anasazi::Operator<ST>               OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
  ST ONE  = SCT::one();

  if (verbose && MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }
  
  // Construct the data for the matrix
  int dim  = 10;
  int dim2 = 10;
  int nnz  = 10;
  double dvals[2*nnz];
  int colptr[dim+1];
  int rowind[nnz+1];
  std::vector<ST> cvals(nnz);
  for (int i = 0; i < dim; ++i) {
      colptr[i] = i+1;
      rowind[i] = i+1;
      dvals[i*2]  = 1.0 + 1.0/dim*i; 
      dvals[i*2+1]= 0.0;
      cvals[i] = ST(1.0+1.0/dim*i,0.0);
  }
  colptr[dim]=nnz+1;
  
  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( new Epetra_Map(dim, 0, *Comm) );

  // Build the problem matrix
  Epetra_CrsMatrix InputMatr(Copy,*Map,3);
  Epetra_CrsMatrix InputMati(Copy,*Map,3);
  for (int i = 0; i < dim; ++i) {
      if (Map->MyGID(i)) {
          for (int j = 0; j < colptr[i+1]-colptr[i]; ++j) {
              int k      = colptr[i]+j-1;
              int ind[1] = { rowind[k]-1 };
              InputMatr.InsertGlobalValues( i, 1, &(dvals[k*2  ]), ind );
              InputMati.InsertGlobalValues( i, 1, &(dvals[k*2+1]), ind );
          }
      }
  }
  InputMatr.FillComplete();
  InputMati.FillComplete();
  Epetra_CrsMatrix_Complex InputMat(InputMatr,InputMati);
  Teuchos::RefCountPtr< Anasazi::EpetraOpComplex > K
    = Teuchos::rcp( new Anasazi::EpetraOpComplex( Teuchos::rcp(&InputMat,false) ),false );

  // Create initial vectors
  int blockSize = 2;
  Teuchos::RefCountPtr< Anasazi::EpetraMultiVecComplex > ivec = Teuchos::rcp( new Anasazi::EpetraMultiVecComplex(*Map,blockSize) );
  ivec->Random();
  std::cout << "Output initial vectors";
  ivec->MvPrint(std::cout);
  std::cout << "\n";

  // CHECK ONE: GetNumberVecs()
  std::cout << "CHECK: GetNumberVecs\n";
  std::cout << "Number of vecs:" << blockSize                 << "\n";
  std::cout << "--------------:" << MVT::GetNumberVecs(*ivec) << "\n";
  std::cout << "\n";

  // CHECK TWO
  std::cout << "CHECK: GetVecLength\n";
  std::cout << "Length of vecs:" << dim << "\n";
  std::cout << "--------------:" << MVT::GetVecLength(*ivec) << "\n";
  std::cout << "\n";

  // CHECK THREE
  dcomplex alpha, beta;
  alpha = dcomplex(2.0,3.0);
  std::cout << "CHECK: MvInit\n";
  std::cout << "Initialize with scalar:" << alpha << "\n";
  ivec->MvPrint(std::cout);
  MVT::MvInit(*ivec, alpha);
  ivec->MvPrint(std::cout);
  std::cout << "\n";

  // CHECK FOUR
  alpha = dcomplex(0.0,0.0);
  MVT::MvInit(*ivec, alpha);
  std::cout << "CHECK: MvRandom\n";
  std::cout << "     Before\n";
  ivec->MvPrint(std::cout);
  MVT::MvRandom(*ivec);
  std::cout << "     After\n";
  ivec->MvPrint(std::cout);
  std::cout << "\n";

  // CHECK FIVE
  alpha = dcomplex(1.0,1.0);
  beta  = dcomplex(1.0,1.0);
  MVT::MvInit(*ivec, alpha);
  std::cout << "CHECK: MvScale\n";
  std::cout << "     Before\n";
  ivec->MvPrint(std::cout);
  MVT::MvScale(*ivec, beta);
  std::cout << "     After\n";
  ivec->MvPrint(std::cout);
  std::cout << "\n";

  std::vector<dcomplex> gamma;
  for (int i = 0; i < dim; ++i)
      gamma.push_back(dcomplex( i, i));
  
/*
  // CHECK SIX
  MVT::MvInit(*ivec, alpha);
  std::cout << "CHECK: MvScale\n";
  std::cout << "     Before\n";
  std::cout << *ivec;
  MVT::MvScale(*ivec, gamma);
  std::cout << "     After\n";
  std::cout << *ivec;
  std::cout << "\n";
*/
  // CHECK SEVEN
  std::vector<double> norms;
  norms.resize(blockSize);
  MVT::MvInit(*ivec, alpha);
  std::cout << "CHECK: MvNorm\n";
  std::cout << "     Before\n";
  ivec->MvPrint(std::cout);
  MVT::MvNorm(*ivec, &norms);
  std::cout << "     Norms:\n";
  for (int i = 0; i < blockSize; ++i)
      std::cout << "Norm[" << i << "]:" << norms[i] << "\n";
  std::cout << "\n";

  // CHECK EIGHT
  MVT::MvInit(*ivec, alpha);
  std::cout << "CHECK: MvDot\n";
  std::cout << "     Before\n";
  ivec->MvPrint(std::cout);
  MVT::MvDot(*ivec, *ivec, &gamma);
  std::cout << "     DotProducts:\n";
  for (int i = 0; i < blockSize; ++i)
      std::cout << "Dot[" << i << "]:" << gamma[i] << "\n";
  std::cout << "\n";

  // CHECK NINE
  MVT::MvInit(*ivec, alpha);
  Teuchos::SerialDenseMatrix<int,dcomplex> B(blockSize,blockSize);
  std::cout << "CHECK: MvTransMv\n";
  std::cout << "     Before\n";
  ivec->MvPrint(std::cout);
  MVT::MvTransMv(alpha,*ivec,*ivec,B);
  std::cout << "     DenseMatrix\n";
  std::cout << B;
  std::cout << "\n";

  // CHECK TEN
  std::cout << "CHECK: MvAddMv\n";
  std::cout << "     Before\n";
  ivec->MvPrint(std::cout);
  MVT::MvAddMv(alpha,*ivec,alpha,*ivec,*ivec);
  std::cout << "     After\n";
  ivec->MvPrint(std::cout);
  std::cout << "\n";

  // CHECK ELEVEN
  std::cout << "CHECK: MvTimesMatAddMv\n";
  std::cout << "     Before\n";
  B.random();
  ivec->MvPrint(std::cout);
  std::cout << B;
  MVT::MvTimesMatAddMv(alpha,*ivec,B,alpha,*ivec);
  std::cout << "     After\n";
  ivec->MvPrint(std::cout);
  std::cout << "\n";

  // CHECK TWELVE
  std::cout << "CHECK: Op::Apply\n";
  std::cout << "     Before\n";
  MVT::MvInit(*ivec, alpha);
  ivec->MvPrint(std::cout);
  std::cout << InputMatr;
  std::cout << InputMati;
  Teuchos::RefCountPtr< Anasazi::EpetraMultiVecComplex > ivec2= 
	       Teuchos::rcp( new Anasazi::EpetraMultiVecComplex(*Map,blockSize) );
  OPT::Apply(*K,*ivec,*ivec2);
  ivec->MvPrint(std::cout);
  ivec2->MvPrint(std::cout);

  return 0;

}	
