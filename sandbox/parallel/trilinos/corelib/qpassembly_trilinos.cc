/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>

#include "Epetra_MpiComm.h"

#include "qcomplex.h"
#include "qassembly.h"
#include "qpassembly_trilinos.h"
#include "coordmatrix.h"
#include "mesh.h"

using namespace std;
extern Epetra_MpiComm* HIQLAB_Comm;

int* Mesh_trilinos_my_nnz(Mesh* mesh, double cx, double cv, double ca,
                                      int is_reduced, Epetra_Map* Map)
{

    int* ir;
    int* jc;
    double* Ax;
    int* nnz;
    int NumGlobalElements = Map->NumGlobalElements();
    int NumMyElements     = Map->NumMyElements();
    int MinMyGID          = Map->MinMyGID();
    int MaxMyGID          = Map->MaxMyGID();
    int Base              = MinMyGID;

    nnz                   = new int[NumMyElements];

    // -- Construct assembler for local process(range restricted)
    CoordMatrix assembler(NumGlobalElements);
    assembler.restrict_range(MinMyGID, MaxMyGID);
    mesh->assemble_dR(&assembler, cx, cv, ca, is_reduced);

    // -- Compute nnz structure of local process
    assembler.pack();
    int ncoord = assembler.get_ncoord();
    ir = new int[NumGlobalElements+1];
    jc = new int[ncoord];
    Ax = new double[ncoord];
    memset(ir, 0, (NumGlobalElements+1) * sizeof(int));
    memset(jc, 0,                ncoord * sizeof(int));
    memset(Ax, 0,                ncoord * sizeof(double));
    assembler.to_sparse_row(ir, jc, Ax, NULL);
    for (int i = 0; i < NumMyElements; ++i)
        nnz[i] = ir[Base+i+1]-ir[Base+i];

    // -- Clean up
    delete[] ir;
    delete[] jc;
    delete[] Ax;

    return nnz;

}


void Mesh_assemble_my_dR_trilinos(Epetra_CrsMatrix* Ar, Epetra_CrsMatrix* Ai, Mesh* mesh,
                                double cx, double cv, double ca, int is_reduced)
{
    int* ir;
    int* jc;
    double* Ax;
    double* Az;
    const Epetra_Map& Map = Ar->RowMap();
    int NumGlobalElements = Map.NumGlobalElements();
    int NumMyElements     = Map.NumMyElements();
    int MinMyGID          = Map.MinMyGID();
    int MaxMyGID          = Map.MaxMyGID();
    int Base              = MinMyGID;

    // -- Construct assembler for local process(range restricted)
    CoordMatrix assembler(NumGlobalElements);
    assembler.restrict_range(MinMyGID, MaxMyGID);
    mesh->assemble_dR(&assembler, cx, cv, ca, is_reduced);

    // -- Compute CRS of local matrix
    assembler.pack();
    int ncoord = assembler.get_ncoord();
    ir = new int[NumGlobalElements+1];
    jc = new int[ncoord];
    Ax = new double[ncoord];
    memset(ir, 0, (NumGlobalElements+1) * sizeof(int));
    memset(jc, 0,                ncoord * sizeof(int));
    memset(Ax, 0,                ncoord * sizeof(double));
    if (Ai) {
        Az = new double[ncoord];
        memset(Az, 0,                ncoord * sizeof(double));
    } else
        Az = NULL;
    assembler.to_sparse_row(ir, jc, Ax, Az);

    // -- Construct matrix
    for (int i = 0; i < NumMyElements; ++i) {
        Ar->InsertGlobalValues( Base+i, ir[Base+i+1]-ir[Base+i],
                          &(Ax[ir[Base+i]]), &(jc[ir[Base+i]]) );
        if (Ai)
            Ai->InsertGlobalValues( Base+i, ir[Base+i+1]-ir[Base+i],
                          &(Az[ir[Base+i]]), &(jc[ir[Base+i]]) );
    }
    HIQLAB_Comm->Barrier();
    Ar->FillComplete();
    Ar->OptimizeStorage();
    if (Ai) {
        Ai->FillComplete();
        Ai->OptimizeStorage();
    }

    // -- Clean up
    delete[] ir;
    delete[] jc;
    delete[] Ax;
    if (Ai) delete[] Az;

}

Epetra_CrsMatrix_Complex* Mesh_assemble_dRz_trilinos(Mesh* mesh, double cx, double cv, double ca,
                           int is_reduced, int is_real)
{
    int* nnz;
    Epetra_CrsMatrix_Complex* A;

    int numid = mesh->get_numid();
    int ndf   = mesh->get_ndf();
    int numnp = mesh->numnp();
    int NumGlobalElements = (is_reduced) ? numid : ndf*numnp;

    // -- Construct Linear Map
    Epetra_Map* Map = new Epetra_Map(NumGlobalElements, 0, *HIQLAB_Comm);

    // -- Compute nnz structure
    nnz = Mesh_trilinos_my_nnz(mesh, cx, cv, ca, is_reduced, Map);

    // -- Construct matrix
    A  = new Epetra_CrsMatrix_Complex(Copy, *Map, nnz, is_real);
    Mesh_assemble_my_dR_trilinos(A, A->get_Az(), mesh, cx, cv, ca, is_reduced);

    delete[] nnz;
    delete Map;

    return A;
}


Epetra_CrsMatrix* Mesh_assemble_dR_trilinos(Mesh* mesh, double cx, double cv, double ca,
                           int is_reduced)
{
    int* nnz;
    Epetra_CrsMatrix* A;

    int numid = mesh->get_numid();
    int ndf   = mesh->get_ndf();
    int numnp = mesh->numnp();
    int NumGlobalElements = (is_reduced) ? numid : ndf*numnp;

    // -- Construct Linear Map
    Epetra_Map* Map = new Epetra_Map(NumGlobalElements, 0, *HIQLAB_Comm);

    // -- Compute nnz structure
    nnz = Mesh_trilinos_my_nnz(mesh, cx, cv, ca, is_reduced, Map);

    // -- Construct matrix
    A  = new Epetra_CrsMatrix(Copy, *Map, nnz);
    Mesh_assemble_my_dR_trilinos(A, NULL, mesh, cx, cv, ca, is_reduced);

    delete[] nnz;
    delete Map;

    return A;
}

Epetra_CrsMatrix* Mesh_assemble_dRi_trilinos(Mesh* mesh, double cx, double cv, double ca,
                           int is_reduced)
{
    int* nnz;
    Epetra_CrsMatrix* A;
    Epetra_CrsMatrix* Ai;

    int numid = mesh->get_numid();
    int ndf   = mesh->get_ndf();
    int numnp = mesh->numnp();
    int NumGlobalElements = (is_reduced) ? numid : ndf*numnp;

    // -- Construct Linear Map
    Epetra_Map* Map = new Epetra_Map(NumGlobalElements, 0, *HIQLAB_Comm);

    // -- Compute nnz structure
    nnz = Mesh_trilinos_my_nnz(mesh, cx, cv, ca, is_reduced, Map);

    // -- Construct matrix
    A  = new Epetra_CrsMatrix(Copy, *Map, nnz);
    Ai = new Epetra_CrsMatrix(Copy, *Map, nnz);
    Mesh_assemble_my_dR_trilinos(A, Ai, mesh, cx, cv, ca, is_reduced);

    delete[] nnz;
    delete A;

    return Ai;
}

void Mesh_assemble_my_R_trilinos(Epetra_Vector* Vr, Epetra_Vector* Vi, Mesh* mesh)
{
    double* fx;
    double* fz;

    int numid = mesh->get_numid();
    const Epetra_BlockMap& Map = Vr->Map();
    int  NumMyElements = Map.NumMyElements();
    int* MyGlobalElements = Map.MyGlobalElements();
    int  MinMyGID = Map.MinMyGID();
    int  MaxMyGID = Map.MaxMyGID();
    int  Base     = MinMyGID;

    // -- Compute local portion of vector
    //FIXME:fix when distributed vector is implemented
    mesh->assemble_R();
    fx = new double[numid];
    memset(fx, 0, numid * sizeof(double));
    mesh->get_f(fx);
    if (Vi) {
        fz = new double[numid];
        memset(fz, 0, numid * sizeof(double));
        mesh->get_fi(fz);
    }

    // -- Construct distributed vector
    Vr->SumIntoGlobalValues( NumMyElements, &(fx[Base]), MyGlobalElements);
    if (Vi)
        Vi->SumIntoGlobalValues( NumMyElements, &(fz[Base]), MyGlobalElements);
    HIQLAB_Comm->Barrier();

    // -- Clean up
    delete[] fx;
    if (Vi)
        delete[] fz;
}

Epetra_Vector_Complex* Mesh_assemble_Rz_trilinos(Mesh* mesh, int is_real)
{
    int numid = mesh->get_numid();
    Epetra_Vector_Complex* V;

    // -- Construct Linear Map
    Epetra_Map Map(numid, 0, *HIQLAB_Comm);

    // -- Construct vector
    V = new Epetra_Vector_Complex(Map, is_real);
    Mesh_assemble_my_R_trilinos(V, V->get_Vz(), mesh);

    return V;
}

Epetra_Vector* Mesh_assemble_R_trilinos(Mesh* mesh)
{
    int numid = mesh->get_numid();
    Epetra_Vector* Vr;

    // -- Construct Linear Map
    Epetra_Map Map(numid, 0, *HIQLAB_Comm);

    // -- Construct vector
    Vr = new Epetra_Vector(Map);
    Mesh_assemble_my_R_trilinos(Vr, NULL, mesh);

    return Vr;
}

Epetra_Vector* Mesh_assemble_Ri_trilinos(Mesh* mesh)
{
    int numid = mesh->get_numid();
    Epetra_Vector* Vr;
    Epetra_Vector* Vi;

    // -- Construct Linear Map
    Epetra_Map Map(numid, 0, *HIQLAB_Comm);

    // -- Construct vector
    Vr = new Epetra_Vector(Map);
    Vi = new Epetra_Vector(Map);
    Mesh_assemble_my_R_trilinos(Vr, Vi, mesh);

    delete Vr;

    return Vi;
}
