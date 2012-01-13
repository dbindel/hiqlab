
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
#include <vector>
#include <algorithm>

#include "trilinos_super_matrix.h"

int s2f_id_mapping(Epetra_Map& Map,
                int* NumMyElements, int** MappedID,
                int (*f)(int) );
int index_mapping1(int id);
int index_mapping2(int id);
Epetra_Map s2f_map(Epetra_Map& Map, int (*f)(int) );
void NumEntriesPerRow(const Epetra_CrsMatrix& A, const Epetra_Map& TargetMap, 
                      int* MyNNZ_B);

struct map_pair {
    int l;
    int g;
};

static int compare_map_pair(const map_pair& pair1, const map_pair& pair2)
{
    if (pair1.g  < pair2.g ) return  1;
    if (pair1.g  > pair2.g ) return  0;
    if (pair1.l  < pair2.l ) return  1;
    if (pair1.l  > pair2.l ) return  0;
    return 0;
}

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // set global dimension of the matrix to 5, could be any number
  int NumGlobalElements = 4;
  
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

//

  // Full Map
  Epetra_Map FullMap(NumGlobalElements*2,0, Comm);

  // Construct shifted A
  Epetra_Map TargetMap_A = s2f_map(FullMap, &index_mapping1);
  int MyNNZ_A[TargetMap_A.NumMyElements()];
  NumEntriesPerRow(A,TargetMap_A,MyNNZ_A);
  Epetra_CrsMatrix Am(Copy,TargetMap_A,MyNNZ_A);
  Epetra_Export Exporter_A(Map,TargetMap_A);
  Am.Export(A,Exporter_A,Add);
  Am.FillComplete();

  // Construct shifted B
  Epetra_Map TargetMap_B = s2f_map(FullMap, &index_mapping2);
  int MyNNZ_B[TargetMap_B.NumMyElements()];
  NumEntriesPerRow(B,TargetMap_B,MyNNZ_B);
  Epetra_CrsMatrix Bm(Copy,TargetMap_B,MyNNZ_B);
  Epetra_Export Exporter_B(Map,TargetMap_B);
  Bm.Export(B,Exporter_B,Add);
  Bm.FillComplete();

  // Construct Full C
  int GID_A[Am.RowMap().NumMyElements()];
  int GID_B[Bm.RowMap().NumMyElements()];
  for (int i = 0; i < Am.RowMap().NumMyElements(); ++i)
      GID_A[i] = index_mapping1(Am.RowMap().GID(i));
  for (int i = 0; i < Bm.RowMap().NumMyElements(); ++i)
      GID_B[i] = index_mapping2(Bm.RowMap().GID(i));

  int MyNNZ_C[FullMap.NumMyElements()];
  memset(MyNNZ_C, 0, FullMap.NumMyElements() * sizeof(int));
  for (int i = 0; i < Am.RowMap().NumMyElements(); ++i)
      MyNNZ_C[FullMap.LID(GID_A[i])] += Am.NumMyEntries(i);
  for (int i = 0; i < Bm.RowMap().NumMyElements(); ++i)
      MyNNZ_C[FullMap.LID(GID_B[i])] += Bm.NumMyEntries(i);


  Epetra_CrsMatrix C(Copy,FullMap,MyNNZ_C);  


  for (int i = 0; i < Am.RowMap().NumMyElements(); ++i) {
      double val_t[Am.NumMyEntries(i)];
      int   id_t [Am.NumMyEntries(i)];
      int   nnz_t;
      Am.ExtractMyRowCopy(i,Am.NumMyEntries(i),nnz_t,val_t,id_t);
      for (int j = 0; j < nnz_t; ++j)
          id_t[j] = index_mapping1(Am.RowMap().GID(id_t[j]));
      C.InsertGlobalValues(GID_A[i],nnz_t,val_t,id_t);
  }

  for (int i = 0; i < Bm.RowMap().NumMyElements(); ++i) {
      double val_t[Bm.NumMyEntries(i)];
      int   id_t [Bm.NumMyEntries(i)];
      int   nnz_t;
      Bm.ExtractMyRowCopy(i,Bm.NumMyEntries(i),nnz_t,val_t,id_t);
      for (int j = 0; j < nnz_t; ++j)
          id_t[j] = index_mapping2(Bm.RowMap().GID(id_t[j]));
      //FIXME: Cannot insert twice into the same place. Also must insert
      //       first to sum.
      C.InsertGlobalValues(GID_B[i],nnz_t,val_t,id_t);
  }
  C.FillComplete();
  std::cout << C;
//

  // build up two distributed vectors q and z, and compute
  // q = A * z
  Epetra_Vector q(A.RowMap());
  Epetra_Vector z(A.RowMap());

  // Fill z with 1's
  z.PutScalar( 1.0 );

  A.Multiply(false, z, q); // Compute q = A*z

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

void NumEntriesPerRow(const Epetra_CrsMatrix& A, const Epetra_Map& TargetMap, 
                      int* MyNNZ_B)
{
    Epetra_Export Exporter(A.RowMap(),TargetMap);

    // Compute the nonzero structure of current distributed matrix
    Epetra_Vector id_A(A.RowMap());
    Epetra_Vector id_B(TargetMap);
    int NumMyElements_A = A.RowMap().NumMyElements();
    double MyNNZ_A[NumMyElements_A];
    int MyID_A[NumMyElements_A];
    for (int i = 0; i <= NumMyElements_A; ++i) {
        MyNNZ_A[i] = (double)(A.NumMyEntries(i));
        MyID_A[i] = i;
    }
    id_A.ReplaceMyValues(NumMyElements_A,MyNNZ_A,MyID_A);
    id_B.Export(id_A,Exporter,Add);

    // Compute the nonzero structure of new distributed matrix
    int NumMyElements_B = TargetMap.NumMyElements();
    double id_B_d[NumMyElements_B];
    id_B.ExtractCopy(id_B_d); 
    for (int i = 0; i <= NumMyElements_B; ++i) {
        MyNNZ_B[i] = (int)(id_B_d[i]);
    }

} 

Epetra_Map s2f_map(Epetra_Map& Map, int (*f)(int) )
{
    int* MappedID;
    int NumMyEntries;
    s2f_id_mapping(Map, &NumMyEntries, &MappedID, f);
    return Epetra_Map(-1,NumMyEntries, MappedID, 0, Map.Comm());
}

int s2f_id_mapping(Epetra_Map& Map,
                int* NumMyElements, int** MappedID,
                int (*f)(int) )
{
    int MyMinGID = Map.MinMyGID();
    int MyMaxGID = Map.MaxMyGID();

    // Find ids within this processor
    std::vector<map_pair> id_pair;
    map_pair id_pair_t;
    int numid = -(*f)(-1);
    *NumMyElements = 0;
    for (int i = 0; i < numid; ++i) {
        int gid = (*f)(i);
        if (MyMinGID <= gid && gid <= MyMaxGID) {
            (*NumMyElements)++;
            id_pair_t.l = i;
            id_pair_t.g = gid;
            id_pair.push_back(id_pair_t);
        }
    }

    // Sort (localmatrixid, globalmatrixid) pair
    sort(id_pair.begin(),id_pair.end(),compare_map_pair);

    *MappedID = new int[*NumMyElements];
    for (int i = 0; i < *NumMyElements; ++i)
        (*MappedID)[i] = id_pair[i].l;

    return 0;
}

int index_mapping1(int id) 
{
    if (0<=id && id<=3) return id;
    return -4;
}

int index_mapping2(int id) 
{

    if (0<=id && id<=3) return id+4;
    return -4;
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
