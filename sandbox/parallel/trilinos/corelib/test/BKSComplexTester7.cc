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

//#define COMPUTE_SIMPLE_PROBLEM
//#define USE_MY_MULTIVEC_OP

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

  // Eigensolver settings
#ifdef COMPUTE_SIMPLE_PROBLEM
  int blockSize = 2;
  const int nev = 1;
  int numBlocks = 3;
#else
  int blockSize = 5;
  const int nev = 4;
  int numBlocks = 8;
#endif
  
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
#ifdef COMPUTE_SIMPLE_PROBLEM
  int dim  = 10;
  int dim2 = 10;
  int nnz  = 10;
  double dvals[nnz];
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
#else
  // Get the data from the HB file
  int dim,dim2,nnz;
  double *dvals;
  int *colptr,*rowind;
  nnz = -1;
  info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,
                              &colptr,&rowind,&dvals);
  if (info == 0 || nnz < 0) {
    if (verbose && MyPID == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  // Convert interleaved doubles to complex values and individual
  // arrays
  std::vector<ST> cvals(nnz);
  for (int ii=0; ii<nnz; ii++) {
    cvals[ii] = ST(dvals[ii*2],dvals[ii*2+1]);
  }
#endif
  
  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( new Epetra_Map(dim, 0, *Comm) );

#ifndef USE_MY_MULTIVEC_OP
  std::cout << "Using Epetra_MultiVector_Complex\n";
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
  Teuchos::RefCountPtr< Anasazi::EpetraMultiVecComplex > ivec = Teuchos::rcp( new Anasazi::EpetraMultiVecComplex(*Map,blockSize) );
  ivec->Random();
#else
  std::cout << "Using MyMultiVector\n";
  // Build the problem matrix
  Teuchos::RefCountPtr< MyBetterOperator<ST> > K 
    = Teuchos::rcp( new MyBetterOperator<ST>(dim,colptr,nnz,rowind,&cvals[0]) );

  // Create initial vectors
  Teuchos::RefCountPtr<MyMultiVec<ST> > ivec = Teuchos::rcp( new MyMultiVec<ST>(dim,blockSize) );
  ivec->MvRandom();
#endif

  // Create eigenproblem
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<ST,MV,OP> > problem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(K,ivec) );
  //
  // Inform the eigenproblem that the operator K is symmetric
  problem->setHermitian(false);
//  problem->setHermitian(true);
  //
  // Set the number of eigenvalues requested
  problem->setNEV( nev );
  //
  // Inform the eigenproblem that you are done passing it information
  boolret = problem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << endl
           << "End Result: TEST FAILED" << endl;	
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
  if (debug) {
    verbosity += Anasazi::Debug;
  }


  // Eigensolver parameters
  int maxRestarts = 100;
  MT tol = 1.0e-6;
  //
  // Create parameter list to pass into the solver manager
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  //
  // Create the solver manager
  Anasazi::BlockKrylovSchurSolMgr<ST,MV,OP> MySolverMgr(problem, MyPL);

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged) {
    testFailed = true;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
  Teuchos::RefCountPtr<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    ostringstream os;
    os.setf(ios::scientific, ios::floatfield);
    os.precision(6);

    // Compute the direct residual
    std::vector<MT> normV( numev );
    Teuchos::SerialDenseMatrix<int,ST> T(numev,numev);
    for (int i=0; i<numev; i++) {
      T(i,i) = sol.Evals[i].realpart;
    }
    Teuchos::RefCountPtr<MV> Kvecs = MVT::Clone( *evecs, numev );

    OPT::Apply( *K, *evecs, *Kvecs );

    MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Kvecs );
    MVT::MvNorm( *Kvecs, &normV );
  
    os << "Direct residual norms computed in BlockKrylovSchurComplex_test.exe" << endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual  " << endl
       << "----------------------------------------" << endl;
    for (int i=0; i<numev; i++) {
      if ( SCT::magnitude(sol.Evals[i].realpart) != SCT::zero() ) {
        normV[i] = SCT::magnitude(normV[i]/sol.Evals[i].realpart);
      }
      os << setw(20) << sol.Evals[i].realpart << setw(20) << normV[i] << endl;
      if ( normV[i] > tol ) {
        testFailed = true;
      }
    }
    if (verbose && MyPID==0) {
      cout << endl << os.str() << endl;
    }

  }
  
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  // Clean up.
  free( dvals );
  free( colptr );
  free( rowind );

  if (testFailed) {
    if (verbose && MyPID==0) {
      cout << "End Result: TEST FAILED" << endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    cout << "End Result: TEST PASSED" << endl;	
  } 
  return 0;

}	
