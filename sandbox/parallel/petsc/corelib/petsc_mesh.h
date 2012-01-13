#ifndef _PETSC_MESH_H
#define _PETSC_MESH_H

#include "petscvec.h"
#include "mesh.h"

void Mesh_SetU_Petsc_Vec(Mesh* mesh, Vec v, int is_reduced);

#endif /* _PETSC_MESH_H */
