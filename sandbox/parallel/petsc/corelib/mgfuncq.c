#define PETSCKSP_DLL

#include <cstring>
#include "src/ksp/pc/impls/mg/mgimpl.h"      /*I "petscksp.h I*/
                                             /*I "petscmg.h" I*/
#include "petscmat.h"
#include "petscmg.h"
#include "qpetscmg.h"

#undef __FUNCT__
#define __FUNCT__ "PCMGSetStationarySmoothers"
/*@C
    PCMGSetSmoothers - Sets stationary smoother for all levels except coarse level.

    Collective

    Input Parameters:

@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetStationarySmoothers(PC pc, PCType pctype, PetscReal scale, PetscInt git, PetscInt lit, PetscInt ndown, PetscInt nup)
{
    PetscErrorCode ierr;
    PC_MG **mg = (PC_MG**)pc->data;
    PetscInt levels;
    KSP ksp;
    PC  ksppc;
    PetscInt i;


    PetscFunctionBegin;
    ierr = PCMGGetLevels(pc, &levels);CHKERRQ(ierr);

    // -- For each level set the type of down smoother and iterations
    for (i = 1; i < levels; ++i) {

        ierr = PCMGGetSmootherDown(pc,i,&ksp);CHKERRQ(ierr);
        ierr = KSPSetType(ksp,KSPRICHARDSON);CHKERRQ(ierr);
        ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,ndown);
//        ierr = KSPRichardsonSetScale(ksp,scale);  // -- Set damping factor
        ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

        ierr = KSPGetPC(ksp,&ksppc);CHKERRQ(ierr);
        ierr = PCSetType(ksppc,pctype);CHKERRQ(ierr);

        // set iterations
        if (!strcmp(pctype,PCSOR)) {
            ierr = PCSORSetIterations(ksppc,git,lit);CHKERRQ(ierr);
        } else if (!strcmp(pctype,"kaczmarz")) {
            ierr = PCKaczmarzSetIterations(ksppc,git,lit);CHKERRQ(ierr);
            ierr = PCKaczmarzSetOmega(ksppc,scale);CHKERRQ(ierr);
        }

    }

    // -- For each level set the type of up smoother and iterations
    for (i = 1; i < levels; ++i) {

        ierr = PCMGGetSmootherUp(pc,i,&ksp);CHKERRQ(ierr);
        ierr = KSPSetType(ksp,KSPRICHARDSON);CHKERRQ(ierr);
        ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,nup);
//        ierr = KSPRichardsonSetScale(ksp,scale);  // -- Set damping factor
        ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

        ierr = KSPGetPC(ksp,&ksppc);CHKERRQ(ierr);
        ierr = PCSetType(ksppc,pctype);CHKERRQ(ierr);

        // set iterations
        if (!strcmp(pctype,PCSOR)) {
            ierr = PCSORSetIterations(ksppc,git,lit);CHKERRQ(ierr);
        } else if (!strcmp(pctype,"kaczmarz")) {
            ierr = PCKaczmarzSetIterations(ksppc,git,lit);CHKERRQ(ierr);
            ierr = PCKaczmarzSetOmega(ksppc,scale);CHKERRQ(ierr);
        }
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCMGSetKrylovSmoothers"
/*@C
    PCMGSetSmoothers - Sets Krylov smoother for all levels except coarse level.

    Collective

    Input Parameters:

@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetKrylovSmoothers(PC pc, KSPType ksptype, PCType pctype, PetscInt ndown, PetscInt nup, PetscReal scale, PetscInt restart, PetscInt git, PetscInt lit)
{
    PetscErrorCode ierr;
    PC_MG **mg = (PC_MG**)pc->data;
    PetscInt levels;
    KSP ksp;
    PC  ksppc;
    PetscInt i;


    PetscFunctionBegin;
    ierr = PCMGGetLevels(pc, &levels);CHKERRQ(ierr);

    // -- For each level set the type of down smoother and iterations
    for (i = 1; i < levels; ++i) {

        ierr = PCMGGetSmootherDown(pc,i,&ksp);CHKERRQ(ierr);
        ierr = KSPSetType(ksp,ksptype);CHKERRQ(ierr);
        ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,ndown);
        ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

        ierr = KSPGetPC(ksp,&ksppc);CHKERRQ(ierr);
        ierr = PCSetType(ksppc,pctype);CHKERRQ(ierr);

        // set restart length if GMRES
        if (!strcmp(ksptype,KSPGMRES) || !strcmp(ksptype,KSPFGMRES) || !strcmp(ksptype,KSPLGMRES)) {
            ierr = KSPGMRESSetRestart(ksp,restart);CHKERRQ(ierr);
        }

        // set iterations
        if (!strcmp(pctype,PCSOR)) {
            ierr = PCSORSetIterations(ksppc,git,lit);CHKERRQ(ierr);
        } else if (!strcmp(pctype,"kaczmarz")) {
            ierr = PCKaczmarzSetIterations(ksppc,git,lit);CHKERRQ(ierr);
            ierr = PCKaczmarzSetOmega(ksppc,scale);CHKERRQ(ierr);
        }

    }

    // -- For each level set the type of up smoother and iterations
    for (i = 1; i < levels; ++i) {

        ierr = PCMGGetSmootherUp(pc,i,&ksp);CHKERRQ(ierr);
        ierr = KSPSetType(ksp,ksptype);CHKERRQ(ierr);
        ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,nup);
        ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

        ierr = KSPGetPC(ksp,&ksppc);CHKERRQ(ierr);
        ierr = PCSetType(ksppc,pctype);CHKERRQ(ierr);

        // set restart length if GMRES
        if (!strcmp(ksptype,KSPGMRES) || !strcmp(ksptype,KSPFGMRES) || !strcmp(ksptype,KSPLGMRES)) {
            ierr = KSPGMRESSetRestart(ksp,restart);CHKERRQ(ierr);
        }

        // set iterations
        if (!strcmp(pctype,PCSOR)) {
            ierr = PCSORSetIterations(ksppc,git,lit);CHKERRQ(ierr);
        } else if (!strcmp(pctype,"kaczmarz")) {
            ierr = PCKaczmarzSetIterations(ksppc,git,lit);CHKERRQ(ierr);
            ierr = PCKaczmarzSetOmega(ksppc,scale);CHKERRQ(ierr);
        }
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCMGConstructDiagonalMassPreconditioner"
PetscErrorCode PETSCKSP_DLLEXPORT PCMGConstructDiagonalMassPreconditioner(Mat M, PetscInt dtype, Mat* Md)
{
    PetscErrorCode ierr;
    PetscInt       Istart, Iend1, i;
    MatType        mtype;
    Vec            vec;

    PetscFunctionBegin;

    // -- Extract info
    ierr = MatGetVecs(M,PETSC_NULL,&vec);CHKERRQ(ierr);
    if (dtype == 0) {
        ierr = MatGetDiagonal(M,vec);CHKERRQ(ierr);
    } else {
        ierr = MatGetRowSum(M,vec);CHKERRQ(ierr);
    }
//    ierr = VecSqrt(vec);CHKERRQ(ierr);
//    ierr = VecReciprocal(vec);CHKERRQ(ierr);
    ierr = MatGetType(M,&mtype);CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(M,&Istart,&Iend1);CHKERRQ(ierr);

    // -- Construct matrix
    ierr = MatCreate(PETSC_COMM_WORLD,Md);CHKERRQ(ierr);
    ierr = MatSetSizes(*Md, Iend1-Istart, Iend1-Istart, PETSC_DETERMINE, PETSC_DETERMINE);CHKERRQ(ierr);
    ierr = MatSetType(*Md, mtype);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(*Md, 1, PETSC_NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(*Md, 1, PETSC_NULL, 0, PETSC_NULL);CHKERRQ(ierr);
    for (i = Istart; i < Iend1; ++i) {
        ierr = MatSetValue(*Md, i, i, 1.0, INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(*Md,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*Md,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    ierr = MatDiagonalScale(*Md, vec, PETSC_NULL);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PCMGSetSmootherDiagonalMassPreconditioning"
/*@C
    PCMGSetSmoothers - Sets diagonal mass preconditioner for smoother for all levels except coarse level.

    Collective

    Input Parameters:

@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetSmootherDiagonalMassPreconditioning(PC pc, Mat M, PetscInt dtype)
{
    PetscErrorCode ierr;
    PC_MG          **mg = (PC_MG**)pc->data;
    PetscInt          n;
    Mat            dA,dB,B,Md;
    PC             pcm;
    MatStructure   uflag;
    PetscInt       i;
    PetscTruth     opsset;

    PetscFunctionBegin;

    ierr = PCMGGetLevels(pc,&n);CHKERRQ(ierr);

    // -- Extract operator and project
    ierr = PCMGConstructDiagonalMassPreconditioner(M, dtype, &Md);CHKERRQ(ierr);

    ierr = KSPGetPC(mg[n-1]->smoothd,&pcm);CHKERRQ(ierr);
    ierr = PCGetOperators(pcm,&dA,PETSC_NULL,&uflag);CHKERRQ(ierr);
    ierr = PCSetOperators(pcm,dA,Md,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = PCSetType(pcm,PCMAT);CHKERRQ(ierr);
    ierr = KSPSetPC(mg[n-1]->smoothd,pcm);CHKERRQ(ierr);

    ierr = KSPGetPC(mg[n-1]->smoothu,&pcm);CHKERRQ(ierr);
    ierr = PCGetOperators(pcm,&dA,PETSC_NULL,&uflag);CHKERRQ(ierr);
    ierr = PCSetOperators(pcm,dA,Md,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = PCSetType(pcm,PCMAT);CHKERRQ(ierr);
    ierr = KSPSetPC(mg[n-1]->smoothu,pcm);CHKERRQ(ierr);

    dB = M;

    for (i = n-2; i > 0; i--) {
        ierr = MatPtAP(dB,mg[i+1]->interpolate,MAT_INITIAL_MATRIX,1.0,&B);CHKERRQ(ierr);
        ierr = PCMGConstructDiagonalMassPreconditioner(B, dtype, &Md);CHKERRQ(ierr);

        ierr = KSPGetPC(mg[i]->smoothd,&pcm);CHKERRQ(ierr);
        ierr = PCGetOperators(pcm,&dA,PETSC_NULL,&uflag);CHKERRQ(ierr);
        ierr = PCSetOperators(pcm,dA,Md,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
        ierr = PCSetType(pcm,PCMAT);CHKERRQ(ierr);
        ierr = KSPSetPC(mg[i]->smoothd,pcm);CHKERRQ(ierr);

        ierr = KSPGetPC(mg[i]->smoothu,&pcm);CHKERRQ(ierr);
        ierr = PCGetOperators(pcm,&dA,PETSC_NULL,&uflag);CHKERRQ(ierr);
        ierr = PCSetOperators(pcm,dA,Md,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
        ierr = PCSetType(pcm,PCMAT);CHKERRQ(ierr);
        ierr = KSPSetPC(mg[i]->smoothu,pcm);CHKERRQ(ierr);

        if  (i!=n-2) {ierr = MatDestroy(dB);CHKERRQ(ierr);}
        dB = B;
    }
    if  (n!=2) {ierr = MatDestroy(dB);CHKERRQ(ierr);}

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PCMGComputeGalerkinOperators"
/*@C
    PCMGComputeGalerkinOperators - Computes the coarse grid operators by
                                   Galerkin projection using interpolation
                                   operator.

    Collective

    Input Paramters:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCMGComputeGalerkinOperators(PC pc, MatType mtype)
{
    PetscErrorCode ierr;
    PC_MG          **mg = (PC_MG**)pc->data;
    PetscInt          n;
    Mat            dA,dB,B;
    MatStructure   uflag;
    PetscInt       i;
    PetscTruth     opsset;

    PetscFunctionBegin;

    ierr = PCMGGetLevels(pc,&n);CHKERRQ(ierr);

    // -- If operator is not set for finest smoother use the one set in PC
    ierr = KSPGetOperatorsSet(mg[n-1]->smoothd,PETSC_NULL,&opsset);CHKERRQ(ierr);
    if (!opsset) {
        ierr = PetscInfo(pc,"Using outer operators to define finest grid operator \n  because PCMGGetSmoother(pc,nlevels-1,&ksp);KSPSetOperators(ksp,...); was not called.\n");CHKERRQ(ierr);
        ierr = KSPSetOperators(mg[n-1]->smoothd,pc->mat,pc->pmat,pc->flag);CHKERRQ(ierr);
        ierr = KSPSetOperators(mg[n-1]->smoothu,pc->mat,pc->pmat,pc->flag);CHKERRQ(ierr);
    }

    // -- Extract operator and project
    ierr = KSPGetOperators(mg[n-1]->smoothd,&dA,&dB,&uflag);CHKERRQ(ierr);
    for (i = n-2; i > -1; i--) {
        ierr = MatPtAP(dB,mg[i+1]->interpolate,MAT_INITIAL_MATRIX,1.0,&B);CHKERRQ(ierr);
        ierr = KSPSetOperators(mg[i]->smoothd,B,B,uflag);CHKERRQ(ierr);
        ierr = KSPSetOperators(mg[i]->smoothu,B,B,uflag);CHKERRQ(ierr);
        if  (i!=n-2) {ierr = PetscObjectDereference((PetscObject)dB);CHKERRQ(ierr);}
        dB = B;
    }
    ierr = PetscObjectDereference((PetscObject)dB);CHKERRQ(ierr);

    // -- Take coarsest operator and change to the desired matrix type
    ierr = KSPGetOperators(mg[0]->smoothd,&dA,&dB,&uflag);CHKERRQ(ierr);
    ierr = MatConvert(dB,mtype,MAT_INITIAL_MATRIX,&B);CHKERRQ(ierr);
    ierr = KSPSetOperators(mg[0]->smoothd,B,B,uflag);CHKERRQ(ierr);

    ierr = PetscObjectDereference((PetscObject)B);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCMGSetDirectCoarseSolve"
/*@C
    PCMGSetDirectCoarseSolve - Sets a direct solve for the coarsest grid

    Collective

    Input Parameters:

@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetDirectCoarseSolve(PC pc, MatType mtype)
{
    PetscErrorCode ierr;
    PC_MG          **mg = (PC_MG**)pc->data;
    PetscInt          n = mg[0]->levels;
    KSP ksp;
    PC  cpc;

    PetscFunctionBegin;

    if (!mg[n-1]->interpolate) SETERRQ(PETSC_ERR_ARG_WRONGSTATE,"Must set interpolates before calling");

    ierr = PCMGGetCoarseSolve(pc,&ksp);CHKERRQ(ierr);
    ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&cpc);CHKERRQ(ierr);
    ierr = PCSetType(cpc,"lu");CHKERRQ(ierr);

    ierr = PCMGComputeGalerkinOperators(pc, mtype);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCMGSetInterpolations"
/*@C
    PCMGSetInterpolations - Sets the prolongators between the grids
                            they are ordered from coarsest to finest,

    Collective

    Input Parameters:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetInterpolations(PC pc, Mat* P, PetscInt nump)
{
    PetscErrorCode ierr;
    PetscInt i;

    PetscFunctionBegin;

    for (i = 0; i < nump; ++i) {
      ierr = PCMGSetInterpolation(pc,i+1,P[i]);CHKERRQ(ierr);
//      ierr = PCMGSetInterpolate(pc,i+1,P[i]);CHKERRQ(ierr); //2.3.2
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCMGSetOperators"
/*@C
    PCMGSetOperators - Sets the operator for each grid
                            they are ordered from coarsest to finest,

    Collective

    Input Parameters:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetOperators(PC pc, Mat* A, PetscInt nump)
{
    PetscErrorCode ierr;
    PC_MG          **mg = (PC_MG**)pc->data;
    PetscInt i;
    Mat dA;
    MatType mtype;
    MatStructure   uflag;

    PetscFunctionBegin;

    // -- Take coarsest operator and change to the desired matrix type
    ierr = KSPGetOperators(mg[0]->smoothd,PETSC_NULL,&dA,&uflag);CHKERRQ(ierr);
    ierr = MatGetType(dA,&mtype);CHKERRQ(ierr);
    ierr = MatConvert(A[0],mtype,MAT_INITIAL_MATRIX,&dA);CHKERRQ(ierr);
    ierr = KSPSetOperators(mg[0]->smoothd,dA,dA,uflag);CHKERRQ(ierr);
    ierr = MatDestroy(dA);CHKERRQ(ierr);
    for (i = 1; i < nump; ++i) {
      ierr = KSPSetOperators(mg[i]->smoothd,A[i],A[i],DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCMGSetOperator"
/*@C
    PCMGSetOperator - Sets the operator for a grid

    Collective

    Input Parameters:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCMGSetOperator(PC pc, PetscInt igrid, Mat A)
{
    PetscErrorCode ierr;
    PC_MG          **mg = (PC_MG**)pc->data;
    PetscInt i;
    Mat dA;
    MatType mtype;
    MatStructure   uflag;

    PetscFunctionBegin;

    if (igrid==0) {
        // -- Take coarsest operator and change to the desired matrix type
        ierr = KSPGetOperators(mg[0]->smoothd,PETSC_NULL,&dA,&uflag);CHKERRQ(ierr);
        ierr = MatGetType(dA,&mtype);CHKERRQ(ierr);
        ierr = MatConvert(A,mtype,MAT_INITIAL_MATRIX,&dA);CHKERRQ(ierr);
        ierr = KSPSetOperators(mg[0]->smoothd,dA,dA,uflag);CHKERRQ(ierr);
    } else {
        ierr = KSPSetOperators(mg[igrid]->smoothd,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCMGGetOperator"
/*@C
    PCMGGetOperator - Gets the operator for a grid

    Collective

    Input Parameters:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCMGGetOperator(PC pc, PetscInt igrid, Mat* A)
{
    PetscErrorCode ierr;
    PC_MG          **mg = (PC_MG**)pc->data;
    PetscInt i;
    Mat dA;
    MatStructure   uflag;

    PetscFunctionBegin;

    ierr = KSPGetOperators(mg[igrid]->smoothd,PETSC_NULL,A,&uflag);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
