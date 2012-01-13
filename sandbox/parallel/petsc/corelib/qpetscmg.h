#ifndef QPETSCMG_H
#define QPETSCMG_H

#include "petscksp.h"
PETSC_EXTERN_CXX_BEGIN

// -- PCMG setters
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetStationarySmoothers(PC pc, PCType pctype, 
                     PetscReal scale, PetscInt git, PetscInt lit, PetscInt ndown, PetscInt nup);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetDirectCoarseSolve(PC pc, MatType mtype);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetInterpolations(PC pc, Mat* P, PetscInt nump);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetOperators(PC pc, Mat* A, PetscInt nump);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetOperator(PC pc, PetscInt i, Mat A);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCMGGetOperator(PC pc, PetscInt i, Mat* A);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetKrylovSmoothers(PC pc, KSPType ksptype, PCType pctype,
                                PetscInt ndown, PetscInt nup, PetscReal scale, PetscInt restart, PetscInt git, PetscInt lit);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetSmootherDiagonalMassPreconditioning(PC pc, Mat M, PetscInt type);

// -- ProjectedPC
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCRegisterProjectedPC();
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetPC(PC pc, PC spc);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetSpaces(PC pc, Vec* Q, Vec* Z, Vec* KiZ, PetscInt maxn);
//EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetQ(PC pc, Vec* Q);
//EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetZ(PC pc, Vec* Z);
//EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetKiZ(PC pc, Vec* KiZ);
//EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetMaxSize(PC pc, PetscInt maxn);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetCurrentSize(PC pc, PetscInt curn);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetComputedSize(PC pc, PetscInt caln);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetMaxSize(PC pc, PetscInt* maxn);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetCurrentSize(PC pc, PetscInt* curn);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetComputedSize(PC pc, PetscInt* caln);

// -- Kaczmarz
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCRegisterKaczmarz();
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCKaczmarzSetIterations(PC pc, PetscInt git, PetscInt lit);
EXTERN PetscErrorCode PETSCKSP_DLLEXPORT PCKaczmarzSetOmega(PC pc, PetscReal omega);

#endif /* QPETSCMG_H */
