#include <iostream>
#include "trilinos_mesh.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_Export.h"

void Mesh_SetU_Epetra_MultiVector(Mesh* mesh, Epetra_MultiVector* emv, int ind)
{

    int lda;
    int NumMyElements;

    // -- Construct a map which has all nodes on root
    const Epetra_Comm& comm = emv->Comm();
    if (comm.MyPID()==0)
        NumMyElements = emv->GlobalLength();
    else
        NumMyElements = 0;
    Epetra_Map sp_map(-1, NumMyElements, 0, comm);
    Epetra_MultiVector sp_vec(sp_map, emv->NumVectors());

    // -- Export Global vector to root process vector
    Epetra_Export exporter(emv->Map(), sp_map);
    sp_vec.Export(*emv, exporter, Add);

    // -- Extract view of the vector
    double* u;
    sp_vec.ExtractView(&u, &lda);

    // -- Copy only for root
    if (comm.MyPID()==0)
        mesh->set_u( &(u[lda*ind]), NULL, NULL);

    // -- Barrier till copy finishes
    comm.Barrier();

}
