
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

// this example creates a tridiagonal matrix of type
//
//     |  2  -1            |
//     | -1   2   -1       | 
// A = |      ...  ... ... |
//     |            -1  2  |

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include <vector>
#include <algorithm>

#include "trilinos_super_matrix.h"
#include "trilinos_super_vector.h"
#include "trilinos_indexmap.h"

#define SIZE_OF_PROBLEM 4

class IndexMapA : public IndexMap {
 public:
   IndexMapA(int sop):sop(sop) {};
   ~IndexMapA() {};
   int inline row(int id) 
   {   
       if (0<=id && id<=sop) return 2*id;
       return -sop;
   };
   int inline col(int id)
   {
       if (0<=id && id<=sop) return 2*id;
       return -sop;
   }
 private:
   int sop;
};
class IndexMapB : public IndexMap{
 public:
   IndexMapB(int sop):sop(sop) {};
   ~IndexMapB() {};
   int inline row(int id) 
   {   
       if (0<=id && id<=sop) return 2*id+1;
       return -sop;
   };
   int inline col(int id)
   {
       if (0<=id && id<=sop) return 2*id+1;
       return -sop;
   }
 private:
   int sop;
};

class IndexMapC : public IndexMap{
 public:
   IndexMapC(int sop):sop(sop) {};
   ~IndexMapC() {};
   int inline row(int id) 
   {   
       if (id%2 == 0) return id/2;
       return -sop;
   };
   int inline col(int id)
   {
       if (id%2 == 0) return id/2;
       return -sop;
   }
 private:
   int sop;
};

int index_mapping1(int id);
int index_mapping2(int id);

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // set global dimension of the matrix to 5, could be any number
  int NumGlobalElements = SIZE_OF_PROBLEM;
  
  // create a map
  Epetra_Map Map(NumGlobalElements,0,Comm);
  
  // local number of rows
  int NumMyElements = Map.NumMyElements();
  
  // get update list
  int * MyGlobalElements = Map.MyGlobalElements( );

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor

  int * NumNz = new int[NumMyElements];

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for ( int i=0; i<NumMyElements; i++)
    if (MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalElements-1)
      NumNz[i] = 2;
    else
      NumNz[i] = 3;

  // Create a Epetra_Matrix
  Epetra_CrsMatrix A(Copy,Map,NumNz);
  Epetra_CrsMatrix B(Copy,Map,NumNz);
  // (NOTE: constructor `Epetra_CrsMatrix A(Copy,Map,3);' was ok too.)
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1, diagonal term 2

  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  double *Values2 = new double[2];
  Values2[0] = -2.0; Values2[1] = -2.0;
  int *Indices = new int[2];
  double two = 2.0;
  double four = 4.0;
  int NumEntries;

  for( int i=0 ; i<NumMyElements; ++i ) {
    if (MyGlobalElements[i]==0) {
	Indices[0] = 1;
	NumEntries = 1;
    } else if (MyGlobalElements[i] == NumGlobalElements-1) {
      Indices[0] = NumGlobalElements-2;
      NumEntries = 1;
    } else {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      NumEntries = 2;
    }
    A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    // Put in the diagonal entry
    A.InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i);

    B.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values2, Indices);
    // Put in the diagonal entry
    B.InsertGlobalValues(MyGlobalElements[i], 1, &four, MyGlobalElements+i);
  }
  
  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  A.FillComplete();
  B.FillComplete();

std::cout << A;
std::cout << B;
//

  // Full Map
  Epetra_Map FullMap(NumGlobalElements*2,0, Comm);

  // Construct super matrix
  Epetra_Super_CrsMatrix D(FullMap);
  IndexMapA ima(NumGlobalElements);
  IndexMapB imb(NumGlobalElements);
  IndexMapC imc(NumGlobalElements);
  D.InsertMatrix(&A, &ima);
  D.InsertMatrix(&B, &imb);
//  D.InsertMatrix(&A, &index_mapping1, &index_mapping2);
//  D.InsertMatrix(&A, &index_mapping1, &index_mapping1);
//  D.InsertMatrix(&B, &index_mapping2, &index_mapping1);
//  D.InsertMatrix(&B, &index_mapping2, &index_mapping2);
  D.ConstructMatrix();
  std::cout << *(D.GetCrsMatrix());

  // build up two distributed vectors q and z, and compute
  // q = A * z
  Epetra_Vector q(D.GetCrsMatrix()->RowMap());
  Epetra_Vector z(D.GetCrsMatrix()->RowMap());

  // Construct super vector
  Epetra_Super_MultiVector f(FullMap,2);
  Epetra_MultiVector fr(Map,2);
  Epetra_MultiVector fi(Map,2);
  Epetra_MultiVector fe(Map,2);
  fr.Random();
  fi.Random();
  f.InsertVector(&fr, &ima);
  f.InsertVector(&fi, &imb);
  f.ConstructVector();
  std::cout << fr;
  std::cout << fi;
  std::cout << *(f.GetMultiVector());

  fe.PutScalar( 3.14 );
  std::cout << fe;
  ExtractMultiVector(f.GetMultiVector(),&fe,&imc);
  std::cout << fe;

  // Fill z with 1's
  z.PutScalar( 1.0 );

  D.GetCrsMatrix()->Multiply(false, z, q); // Compute q = A*z

  double dotProduct;
  z.Dot( q, &dotProduct );

  if( Comm.MyPID() == 0 ) 
    cout << "q dot z = " << dotProduct << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  delete NumNz;
  
  return( EXIT_SUCCESS );

} /* main */

int index_mapping1(int id) 
{
    if (0<=id && id<=SIZE_OF_PROBLEM) return id;
    return -SIZE_OF_PROBLEM;
}

int index_mapping2(int id) 
{

    if (0<=id && id<=SIZE_OF_PROBLEM) return id+SIZE_OF_PROBLEM;
    return -SIZE_OF_PROBLEM;

/*
    if (0<=id && id<=3) {
       if (id==0) return 5;
       if (id==1) return 2;
       if (id==2) return 7;
       if (id==3) return 0;
    }
    return -4;
*/
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure Didasko with:\n"
       "--enable-epetra");

  return 0;
}
#endif
