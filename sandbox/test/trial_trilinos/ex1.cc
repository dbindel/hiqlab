// @HEADER
// ***********************************************************************
// 
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
// 
// ***********************************************************************
// @HEADER

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_AMESOS) && defined(HAVE_DIDASKO_TEUCHOS) && defined(HAVE_DIDASKO_TRIUTILS)

#include "Epetra_ConfigDefs.h"
#include "Amesos_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Amesos.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_vector.h"
#include "trilinos_epetra_linearproblem.h"
#include "trilinos_amesos_base.h"

using namespace Trilinos_Util;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // initialize an Gallery object
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 4); //must be a square number
 
  // get pointers to the linear problem, containing matrix, LHS and RHS.
  // if you need to access them, you can for example uncomment the following
  // code:
  // Epetra_CrsMatrix* Matrix = Gallery.GetMatrix();
  // Epetra_MultiVector* LHS = Gallery.GetStartingSolution();
  // Epetra_MultiVector* RHS = Gallery.GetRHS();
  //
  // NOTE: StartingSolution and RHS are pointers to Gallery's internally stored
  // vectors. Using StartingSolution and RHS, we can verify the residual
  // after the solution of the linear system. However, users may define as well
  // their own vectors for solution and RHS. 
  
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  // -- Epetra_CrsMatrix_Complex Matrix
  Epetra_CrsMatrix* Ar = dynamic_cast<Epetra_CrsMatrix *>(Problem->GetMatrix());
  Epetra_CrsMatrix_Complex Az(*Ar,*Ar); 

  // -- Epetra_MultiVector_Complex RHS,LHS
  Epetra_MultiVector* Br = Problem->GetRHS();
  Epetra_MultiVector_Complex Bz(*Br, *Br);
  Epetra_MultiVector_Complex Xz(*Br, *Br);
std::cout << Bz;
std::cout << *(Bz.get_Vz());

  // -- Epetra_LinearProblem_Complex
  Epetra_LinearProblem_Complex ELP(&Az, &Xz, &Bz);

  // initialize Amesos solver:
  // `Solver' is the pointer to the Amesos solver
  // (note the use of the base class Amesos_BaseSolver)
  Amesos_BaseSolver_Complex*  ABS;
//  Amesos_BaseSolver* Solver;
  // Amesos_Factory is the function class used to create the solver.
  // This class contains no data.
  Amesos_Complex Amesos_Factory;
//  Amesos Amesos_Factory;

  // empty parameter list
  Teuchos::ParameterList List;
  
  // may also try: "Amesos_Umfpack", "Amesos_Lapack", ...
  string SolverType = "Amesos_Mumps";
//  string SolverType = "Amesos_Klu";
  
  ABS = Amesos_Factory.Create(SolverType, ELP);
//  Solver = Amesos_Factory.Create(SolverType, ELP);
//  Solver = Amesos_Factory.Create(SolverType, *Problem);
  // Amesos_Factory returns 0 is the selected solver is not
  // available
  assert(ABS);
//  assert (Solver);

  // start solving
  ABS->SymbolicFactorization();
  ABS->NumericFactorization();
  ABS->Solve();

  // verify that residual is really small  
  double residual, diff;

//  Gallery.ComputeResidual(&residual);
//  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);
std::cout << Xz;
std::cout << *(Xz.get_Vz());

  if( Comm.MyPID() == 0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }

  // delete Solver
  delete ABS;
    
  if (residual > 1e-5)
    exit(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  exit(EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure Didasko with:\n"
       "--enable-epetra\n"
       "--enable-teuchos\n"
       "--enable-triutils\n"
       "--enable-amesos");

  return 0;
}

#endif

