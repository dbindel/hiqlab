#ifndef _QPETSC_PC_H
#define _QPETSC_PC_H

#include <string>

#include "petsc.h"
#include "petscpc.h"
#include "petscksp.h"
#include "mesh.h"

int PCSetCoordinatesFromMesh(PC pc, Mesh* mesh);

#endif /* _QPETSC_PC_H */
