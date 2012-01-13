#ifndef _QPETSC_H
#define _QPETSC_H

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "mesh.h"

inline PetscViewer get_PETSC_VIEWER_STDOUT_WORLD()
{
    return PETSC_VIEWER_STDOUT_WORLD;
}

const char* qPetscErrorMessage(int errnum);

PetscViewer qPetscViewerASCIIOpen(const char* name);
PetscViewer qPetscViewerBinaryOpen(const char* name, int type);
int qPetscViewerSetFormat(PetscViewer viewer, int format);

Vec qVecCreate();
int qVecNorm(Vec x, int type, double* val = 0);
Vec qVecLoad(PetscViewer viewer, int outtype);
Mat qMatLoad(PetscViewer viewer, int outtype);

KSP qKSPCreate();
int qKSPSetOperators(KSP ksp, Mat A, Mat B, int sameflag);
PC  qKSPGetPC(KSP ksp);
int qKSPTrueResidualNorm(KSP ksp, int type, double* rnorm);

int PCSetCoordinatesFromMesh(PC pc, Mesh* mesh);

#endif /* _QPETSC_H */
