
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

// Epetra_Export classes
// This code should be run with at least two processes

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
#include "Epetra_Export.h"
#include <vector>

int foo_mapping(Epetra_Map Map,
                int* NumMyElements, int** MappedID,
                int (*f)(int) );
int index_mapping1(int id);
int index_mapping2(int id);

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumGlobalElements = 4; // global dimension of the problem

  int NumMyElements; // local nodes
  Epetra_IntSerialDenseVector MyGlobalElements;

  if( Comm.MyPID() == 0 ) {
    NumMyElements = 2;
    MyGlobalElements.Size(NumMyElements);
    MyGlobalElements[0] = 0;
    MyGlobalElements[1] = 1;
  } else {
    NumMyElements = 2;
    MyGlobalElements.Size(NumMyElements);
    MyGlobalElements[0] = 2;
    MyGlobalElements[1] = 3;
  }

  // create a map
  Epetra_Map Map(-1,MyGlobalElements.Length(),
		 MyGlobalElements.Values(),0, Comm);

  // create a vector based on map
  Epetra_Vector x(Map);
  Epetra_Vector z(Map);
  for( int i=0 ; i<NumMyElements ; ++i ) {
    x[i] = 10*( 2*Comm.MyPID()+1+i );
    z[i] = 20*( 2*Comm.MyPID()+1+i );
  }
  cout << x;
  cout << z;

  Epetra_Map FullMap(NumGlobalElements*2,0, Comm);

  int* MappedID;
  foo_mapping(FullMap, &NumMyElements, &MappedID,
               &index_mapping1);
  MyGlobalElements.Size(NumMyElements);
  for (int i = 0; i < NumMyElements; ++i)
      MyGlobalElements[i] = MappedID[i];
  Epetra_Map TargetMap1(-1,MyGlobalElements.Length(),
		 MyGlobalElements.Values(),0, Comm);

  foo_mapping(FullMap, &NumMyElements, &MappedID,
               &index_mapping2);
  MyGlobalElements.Size(NumMyElements);
  for (int i = 0; i < NumMyElements; ++i)
      MyGlobalElements[i] = MappedID[i];
  Epetra_Map TargetMap2(-1,MyGlobalElements.Length(),
		 MyGlobalElements.Values(),0, Comm);

  Epetra_Export Exporter1(Map,TargetMap1);
  Epetra_Export Exporter2(Map,TargetMap2);

  // work on vectors
  Epetra_Vector y(TargetMap1);
  Epetra_Vector w(TargetMap2);

  y.Export(x,Exporter1,Add);
  w.Export(z,Exporter2,Add);

  cout << y;
  cout << w;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);

}

int foo_mapping(Epetra_Map Map,
                int* NumMyElements, int** MappedID,
                int (*f)(int) )
{
    int MyMinGID = Map.MinMyGID();
    int MyMaxGID = Map.MaxMyGID();

    std::vector<int> mapped_id;
    int numid = -(*f)(-1);
    *NumMyElements = 0;
    for (int i = 0; i < numid; ++i) {
        int gid = (*f)(i);
        if (MyMinGID <= gid && gid <= MyMaxGID) {
            (*NumMyElements)++;
            mapped_id.push_back(i);
        }
    }

    *MappedID = new int[*NumMyElements];
    for (int i = 0; i < *NumMyElements; ++i)
        (*MappedID)[i] = mapped_id[i];

    return 0;
}

int index_mapping1(int id) 
{
    if (0<=id && id<=3) return id;
    return -4;
}
int index_mapping2(int id) 
{
    if (0<=id && id<=3) {
       if (id==0) return 5;
       if (id==1) return 0;
       if (id==2) return 7;
       if (id==3) return 2;
    }
    return -4;
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
