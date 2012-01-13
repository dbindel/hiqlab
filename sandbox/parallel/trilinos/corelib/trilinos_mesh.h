#ifndef _TRILINOS_MESH_H
#define _TRILINOS_MESH_H

#include "Epetra_MultiVector.h"
#include "mesh.h"

void Mesh_SetU_Epetra_MultiVector(Mesh* mesh, Epetra_MultiVector* emv, int ind=0);


#endif /* _TRILINOS_MESH_H */
