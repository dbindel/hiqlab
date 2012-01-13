#ifndef _MATSHELLAB_H
#define _MATSHELLAB_H

#include "petscksp.h"
PETSC_EXTERN_CXX_BEGIN

/*
     Define a context of the user-provided MATSHELL 

     MatShellAB
*/
typedef struct {
    Mat A;
    Mat B;
    Vec v;
    PetscScalar alpha;
    PetscScalar beta;
} MatShellAB;

EXTERN PetscErrorCode PETSCKSP_DLLEXPORT MatShellABCreate(MatShellAB** matshell);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT MatShellABSetOperators(MatShellAB* shell, Mat A, Mat B);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT MatShellABSetShift(MatShellAB* shell, PetscScalar alpha, PetscScalar beta);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT MatShellABApply(Mat op, Vec x, Vec y);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT MatShellABApplyTranspose(Mat op, Vec x, Vec y);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT MatShellABDestroy(MatShellAB* shell);

PETSC_EXTERN_CXX_END
#endif /* MATSHELLAB_H */
