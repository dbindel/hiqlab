#define PETSCKSP_DLL

/* ----------------------------------------------------------------------------

    This file implements a projected PC preconditioner in PTESc as part of PC.
    The preconditioner K is applied on the subspace projected by two subspaces
    spanned by the columns of matrices Q and Z.

    \tilde{K} = (I - ZZ^*)K(I - QQ^I)

    The following basic routines are provided for the preconditioner.

      PCCreate_ProjectedPC()          - Creates the preconditioner context
      PCSetFromOptions_ProjectedPC()  - Sets runtime options
      PCApply_ProjectedPC()           - Applies the preconditioner
      PCDestroy_ProjectedPC()         - Destroys the preconditioner context

    These routines are actually called via the common user interaface routines
    PCCreate(), PCSetFromOptions(), PCApply(), and PCDestroy(), so the application
    code interface remains identical for all preconditioners.

    Another key routine is:

      PCSetUp_ProjectedPC()           - Prepares the use of the preconditioner by
                                        setting data structures and options. The interface
                                        routine PCSetUp() is no usually called directly by
                                        the user, but instead by PCApply() if necessary.

    Additional basic routines are:

      PCView_ProjectedPC()            - Prints details of runtime options that have actually
                                        been used. This is called by the interface routine
                                        PCView().

   ----------------------------------------------------------------------------*/

/*
    Include files needed for the ProjectedPC preconditioner:

      pcimpl.h - private include file intended for use by all preconditioners
*/

#include "private/pcimpl.h"   /*I "petscpc.h" I*/

/*
    Private context (data structure) for the ProjectedPC preconditioner.
*/
typedef struct {
    PC               pc;      /* the preconditioner context */
    Vec*              Q;      /* vectors of the right projection space to project out 
                                 they must be orthogonal */
    Vec*              Z;      /* vectors of the left  projection space to project out
                                 they must be orthogonal */
    Vec*            KiZ;      /* vectors of the left  projection space to project out 
                                 with PC applied */
    PetscInt       maxn;      /* maximum number of vectors */
    PetscInt       curn;      /* current number of vectors */
    PetscInt       caln;      /* computed number of vectors of KiZ */
    PetscScalar*  QhKiZ;      /* the matrix Q^*K^{-1}Z in column major format 
                                 the lda is given as maxn */
    PetscScalar*  QhKix;      /* the vector Q^*K^{-1}x in column major format
                                 the lda is given as maxn */ 
} PC_ProjectedPC;

EXTERN_C_BEGIN
EXTERN PetscErrorCode PCCreate_ProjectedPC(PC);
EXTERN PetscErrorCode PCProjectedPCSetPC_ProjectedPC(PC,PC);
EXTERN PetscErrorCode PCProjectedPCSetSpaces_ProjectedPC(PC,Vec*,Vec*,Vec*,PetscInt);
EXTERN PetscErrorCode PCProjectedPCSetQ_ProjectedPC(PC,Vec*);
EXTERN PetscErrorCode PCProjectedPCSetZ_ProjectedPC(PC,Vec*);
EXTERN PetscErrorCode PCProjectedPCSetKiZ_ProjectedPC(PC,Vec*);
EXTERN PetscErrorCode PCProjectedPCSetMaxSize_ProjectedPC(PC,PetscInt);
EXTERN PetscErrorCode PCProjectedPCSetCurrentSize_ProjectedPC(PC,PetscInt);
EXTERN PetscErrorCode PCProjectedPCSetComputedSize_ProjectedPC(PC,PetscInt);
EXTERN PetscErrorCode PCProjectedPCGetPC_ProjectedPC(PC,PC*);
EXTERN PetscErrorCode PCProjectedPCGetQ_ProjectedPC(PC,Vec**);
EXTERN PetscErrorCode PCProjectedPCGetZ_ProjectedPC(PC,Vec**);
EXTERN PetscErrorCode PCProjectedPCGetKiZ_ProjectedPC(PC,Vec**);
EXTERN PetscErrorCode PCProjectedPCGetMaxSize_ProjectedPC(PC,PetscInt*);
EXTERN PetscErrorCode PCProjectedPCGetCurrentSize_ProjectedPC(PC,PetscInt*);
EXTERN PetscErrorCode PCProjectedPCGetComputedSize_ProjectedPC(PC,PetscInt*);
EXTERN_C_END
static PetscErrorCode PCSetUp_ProjectedPC(PC);
static PetscErrorCode PCApply_ProjectedPC(PC,Vec,Vec);
static PetscErrorCode PCDestroy_ProjectedPC(PC);

EXTERN_C_BEGIN
#ifdef PETSC_USE_COMPLEX
int zgesv_(const PetscInt& N, const PetscInt& nrhs,
           PetscScalar* A, const PetscInt& ldA, PetscInt* ipiv,
           PetscScalar* B, const PetscInt& ldB, PetscInt& info);
#else
int dgesv_(const PetscInt& N, const PetscInt& nrhs,
           PetscScalar* A, const PetscInt& ldA, PetscInt* ipiv,
           PetscScalar* B, const PetscInt& ldB, PetscInt& info);
#endif
EXTERN_C_END


/* ---------------------------------------------------------------------------- */
/*
    PCCreate_ProjectedPC - Creates a ProjectedPC context, PC_ProjectedPC,
    and sets this as the private data within the generic preconditioning context, PC,
    that was created within PCCreate().

    Input Parameter:

      pc - the preconditioner context

    Application Interface Routine: PCCreate()
*/

/*MC
    PCPROJECTEDPC - Projected version of a given preconditioner PC

    Options Database Key:

    Level: intermediate

    Concepts: preconditioners

    Notes:

.seealso: PCCreate(), PCSetType(), PCType() (for a list of available types), PC
M*/

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCCreate_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCCreate_ProjectedPC(PC pc)
{
    PC_ProjectedPC     *ppc;
    PetscErrorCode      ierr;

    PetscFunctionBegin;

    // -- Creates the private data structure for this preconditioner and 
    //    attach it to the PC object.
    ierr     = PetscNew(PC_ProjectedPC,&ppc);CHKERRQ(ierr);
    pc->data = (void*)ppc;

    // -- Logs the memory usage; not required but allows PETSc to monitor
    //    how much memory is being used for various purposes.
    ierr = PetscLogObjectMemory(pc,sizeof(PC_ProjectedPC)); CHKERRQ(ierr);

    // -- Initialize the pointers to vectors to zero
    ppc->pc     = 0;
    ppc->Q      = 0;
    ppc->Z      = 0;
    ppc->KiZ    = 0;

    ppc->maxn   = 0;
    ppc->curn   = 0;
    ppc->caln   = 0;
    ppc->QhKiZ  = 0;
    ppc->QhKix  = 0;

    // -- Set the pointers for the functions that are provided above.
    pc->ops->apply                 = PCApply_ProjectedPC;
    pc->ops->applytranspose        = PCApply_ProjectedPC;
    pc->ops->setup                 = PCSetUp_ProjectedPC;
    pc->ops->destroy               = PCDestroy_ProjectedPC;
    pc->ops->setfromoptions        = 0; //PCSetFromOptions_ProjectedPC;
    pc->ops->view                  = 0; //PCView_ProjectedPC;
    pc->ops->applyrichardson       = 0;
    pc->ops->applysymmetricleft    = 0;
    pc->ops->applysymmetricright   = 0;

    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCSetPC_C",
					     "PCProjectedPCSetPC_ProjectedPC",
					     PCProjectedPCSetPC_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCSetSpaces_C",
					     "PCProjectedPCSetSpaces_ProjectedPC",
					     PCProjectedPCSetSpaces_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCSetQ_C",
					     "PCProjectedPCSetQ_ProjectedPC",
					     PCProjectedPCSetQ_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCSetZ_C",
					     "PCProjectedPCSetZ_ProjectedPC",
					     PCProjectedPCSetZ_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCSetKiZ_C",
					     "PCProjectedPCSetKiZ_ProjectedPC",
					     PCProjectedPCSetKiZ_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCSetMaxSize_C",
					     "PCProjectedPCSetMaxSize_ProjectedPC",
					     PCProjectedPCSetMaxSize_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCSetCurrentSize_C",
					     "PCProjectedPCSetCurrentSize_ProjectedPC",
					     PCProjectedPCSetCurrentSize_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCSetComputedSize_C",
					     "PCProjectedPCSetComputedSize_ProjectedPC",
					     PCProjectedPCSetComputedSize_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCGetPC_C",
					     "PCProjectedPCGetPC_ProjectedPC",
					     PCProjectedPCGetPC_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCGetQ_C",
					     "PCProjectedPCGetQ_ProjectedPC",
					     PCProjectedPCGetQ_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCGetZ_C",
					     "PCProjectedPCGetZ_ProjectedPC",
					     PCProjectedPCGetZ_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCGetKiZ_C",
					     "PCProjectedPCGetKiZ_ProjectedPC",
					     PCProjectedPCGetKiZ_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCGetMaxSize_C",
					     "PCProjectedPCGetMaxSize_ProjectedPC",
					     PCProjectedPCGetMaxSize_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCGetCurrentSize_C",
					     "PCProjectedPCGetCurrentSize_ProjectedPC",
					     PCProjectedPCGetCurrentSize_ProjectedPC);CHKERRQ(ierr);
    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc,"PCProjectedPCGetComputedSize_C",
					     "PCProjectedPCGetComputedSize_ProjectedPC",
					     PCProjectedPCGetComputedSize_ProjectedPC);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
EXTERN_C_END

/* ---------------------------------------------------------------------------- */
/*
    PCDestroy_ProjectedPC - Destroys the private context for the ProjectedPC preconditioner
    that was created with PCCreate_ProjectedPC().

    Input parameter:
      pc - the preconditioner context

    Application Interface Routine: PCDestroy()
*/
#undef __FUNCT__
#define __FUNCT__ "PCDestroy_ProjectedPC"
static PetscErrorCode PCDestroy_ProjectedPC(PC pc)
{
    PC_ProjectedPC     *ppc = (PC_ProjectedPC*)pc->data;
    PetscErrorCode      ierr;
    PetscInt            i;
    PetscObject         po;
    PetscInt            trefct;

    PetscFunctionBegin;
    if (ppc->pc) {
        ierr = PCDestroy(ppc->pc);CHKERRQ(ierr);
    }
    if (ppc->Q) {

        trefct = 0;
        for (i = 0; i < ppc->maxn; ++i) {
            if (ppc->Q[i]) {
                po = (PetscObject)ppc->Q[i];
                ierr = VecDestroy(ppc->Q[i]);CHKERRQ(ierr);
                trefct += po->refct;
            }
        }

        if (trefct <= 0) {
            ierr = PetscFree(ppc->Q);CHKERRQ(ierr);
        }
    }
    if (ppc->Z) {

        trefct = 0;
        for (i = 0; i < ppc->maxn; ++i) {
            if (ppc->Z[i]) {
                po = (PetscObject)ppc->Z[i];
                ierr = VecDestroy(ppc->Z[i]);CHKERRQ(ierr);
                trefct += po->refct;
            }
        }    

        if (trefct <= 0) {
            ierr = PetscFree(ppc->Z);CHKERRQ(ierr);
        }
    }
    if (ppc->KiZ) {

        trefct = 0;
        for (i = 0; i < ppc->maxn; ++i) {
            if (ppc->KiZ[i]) {
                po = (PetscObject)ppc->KiZ[i];
                ierr = VecDestroy(ppc->KiZ[i]);CHKERRQ(ierr);
                trefct += po->refct;
            }
        }

        if (trefct <= 0) {
            ierr = PetscFree(ppc->KiZ);CHKERRQ(ierr);
        }
    }
    if (ppc->QhKiZ) {
        ierr = PetscFree(ppc->QhKiZ); CHKERRQ(ierr);
    }
    if (ppc->QhKix) {
        ierr = PetscFree(ppc->QhKix); CHKERRQ(ierr);
    }
 
    // -- Free the private data structure that was hanign off the PC
    ierr = PetscFree(ppc); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

/* ---------------------------------------------------------------------------- */
/*
    PCSetUp_ProjectedPC - Prepares the use of the Jacobi preconditioner
                          by setting data structures and options.

    Input Parameter:
      pc - the preconditioner context

    Application Interface Routine: PCSetUp()

    Notes:
    The interface routine PCSetUp() is not usually called directly by 
    the user, but instead is called by PCApply() if necessary.
*/
#undef __FUNCT__
#define __FUNCT__ "PCSetUp_ProjectedPC"
static PetscErrorCode PCSetUp_ProjectedPC(PC pc)
{
    PC_ProjectedPC    *ppc = (PC_ProjectedPC*)pc->data;
    PetscErrorCode     ierr;
 
    PetscFunctionBegin;   
    if (ppc->maxn>0) {
      ierr = PetscMalloc(ppc->maxn*ppc->maxn*sizeof(PetscScalar),&ppc->QhKiZ);CHKERRQ(ierr);
      ierr = PetscMalloc(          ppc->maxn*sizeof(PetscScalar),&ppc->QhKix);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------- */
/*
    PCApply_ProjectedPC - Applies the ProjectedPC preconditioner to a vector.

    Input:
      pc - the preconditioner context
      x  - input vector

    Output:
      y  - output vector

    Application Interface Routine: PCApply()
*/
#undef __FUNCT__
#define __FUNCT__ "PCApply_ProjectedPC"
static PetscErrorCode PCApply_ProjectedPC(PC pc, Vec x, Vec y)
{
    PetscErrorCode     ierr;
    PC_ProjectedPC    *ppc = (PC_ProjectedPC*)pc->data;
    Vec*               Q   = ppc->Q;
    Vec*               Z   = ppc->Z;
    Vec*             KiZ   = ppc->KiZ;

    PetscScalar*    QhKiZ  = ppc->QhKiZ;
    PetscScalar*    QhKix  = ppc->QhKix;
    PetscScalar*    QhKiZt;
    PetscInt*       ipiv;
    PetscInt        info;

    PetscInt           i;
    PetscInt          curn = ppc->curn;
    PetscInt          caln = ppc->caln;
    PetscInt          lda  = ppc->maxn;
    PetscInt          maxn = ppc->maxn;
    PetscScalar       val;

    PetscFunctionBegin;
    
    // -- Apply the pc context in this preconditioner
    ierr = PCApply(ppc->pc,x,y);CHKERRQ(ierr);

    if (Q && Z && KiZ) {

        // -- Compute -QhKix
        ierr = VecMDot(y,curn,Q,QhKix);CHKERRQ(ierr);
        for (i = 0; i < curn; ++i)
            QhKix[i] = -QhKix[i];

        // -- Update KiZ and QhKiZ if not all computed
        while (caln < curn) {
            ierr = PCApply(ppc->pc,Z[caln],KiZ[caln]);CHKERRQ(ierr);
  
            for (i = 0; i < caln; ++i) {
	      ierr = VecDot(KiZ[i],Q[caln],&val);CHKERRQ(ierr);
              QhKiZ[caln + lda*i] = val;
            }
            ierr = VecMDot(KiZ[caln],caln+1,Q,&QhKiZ[lda*caln]);

            caln++;
        }
        ppc->caln = caln;

        // -- Compute alpha = (QhKiZ)\-QhKix with zgesv(dgesv)
        ierr = PetscMalloc(maxn*maxn*sizeof(PetscScalar),&QhKiZt);CHKERRQ(ierr);
        ierr = PetscMalloc(     caln*sizeof(PetscInt),   &ipiv);  CHKERRQ(ierr);
        for (i = 0; i < maxn*maxn; ++i)
            QhKiZt[i] = QhKiZ[i];
#ifdef PETSC_USE_COMPLEX
        ierr = zgesv_(caln,1,QhKiZt,lda,ipiv,QhKix,lda,info);
#else
        ierr = dgesv_(caln,1,QhKiZt,lda,ipiv,QhKix,lda,info);
#endif
        ierr = PetscFree(QhKiZt); CHKERRQ(ierr);
        ierr = PetscFree(ipiv); CHKERRQ(ierr);

        // -- Compute Kix - Z * alpha
        ierr = VecMAXPY(y,curn,QhKix,KiZ);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetPC_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetPC_ProjectedPC(PC pc, PC spc)
{
    PC_ProjectedPC *ppc;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    ppc->pc = spc;
    ierr = PetscObjectReference((PetscObject)spc);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetSpaces_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetSpaces_ProjectedPC(PC pc, Vec* Q,Vec* Z, Vec* KiZ, PetscInt maxn)
{
    PC_ProjectedPC *ppc;
    PetscErrorCode ierr;
    PetscInt i;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    ppc->Q   = Q;
    ppc->Z   = Z;
    ppc->KiZ = KiZ;
    ppc->maxn= maxn;

    for (i = 0; i < maxn; ++i) {
        ierr = PetscObjectReference((PetscObject)ppc->Q[i]  );CHKERRQ(ierr);
        ierr = PetscObjectReference((PetscObject)ppc->Z[i]  );CHKERRQ(ierr);
        ierr = PetscObjectReference((PetscObject)ppc->KiZ[i]);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetQ_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetQ_ProjectedPC(PC pc, Vec* Q)
{
    PC_ProjectedPC *ppc;
    PetscErrorCode ierr;
    PetscInt i;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    ppc->Q = Q;

    for (i = 0; i < ppc->maxn; ++i) {
        ierr = PetscObjectReference((PetscObject)ppc->Q[i]);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetZ_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetZ_ProjectedPC(PC pc, Vec* Z)
{
    PC_ProjectedPC *ppc;
    PetscErrorCode ierr;
    PetscInt i;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    ppc->Z = Z;

    for (i = 0; i < ppc->maxn; ++i) {
        ierr = PetscObjectReference((PetscObject)ppc->Z[i]);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetKiZ_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetKiZ_ProjectedPC(PC pc, Vec* KiZ)
{
    PC_ProjectedPC *ppc;
    PetscErrorCode ierr;
    PetscInt i;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    ppc->KiZ = KiZ;

    for (i = 0; i < ppc->maxn; ++i) {
        ierr = PetscObjectReference((PetscObject)ppc->KiZ[i]);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetMaxSize_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetMaxSize_ProjectedPC(PC pc, PetscInt maxn)
{
    PC_ProjectedPC *ppc;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    ppc->maxn = maxn;
    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetCurrentSize_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetCurrentSize_ProjectedPC(PC pc, PetscInt curn)
{
    PC_ProjectedPC *ppc;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    ppc->curn = curn;
    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetComputedSize_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetComputedSize_ProjectedPC(PC pc, PetscInt caln)
{
    PC_ProjectedPC *ppc;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    ppc->caln = caln;
    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCGetPC_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetPC_ProjectedPC(PC pc, PC *spc)
{
    PC_ProjectedPC *ppc;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    *spc = ppc->pc;
    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCGetQ_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetQ_ProjectedPC(PC pc, Vec** Q)
{
    PC_ProjectedPC *ppc;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    *Q = ppc->Q;
    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCGetZ_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetZ_ProjectedPC(PC pc, Vec** Z)
{
    PC_ProjectedPC *ppc;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    *Z = ppc->Z;
    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCGetKiZ_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetKiZ_ProjectedPC(PC pc, Vec** KiZ)
{
    PC_ProjectedPC *ppc;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    *KiZ = ppc->KiZ;
    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCGetMaxSize_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetMaxSize_ProjectedPC(PC pc, PetscInt* maxsize)
{
    PC_ProjectedPC *ppc;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    *maxsize = ppc->maxn;
    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCGetCurrentSize_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetCurrentSize_ProjectedPC(PC pc, PetscInt* curn)
{
    PC_ProjectedPC *ppc;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    *curn = ppc->curn;
    PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCGetComputedSize_ProjectedPC"
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetComputedSize_ProjectedPC(PC pc, PetscInt* caln)
{
    PC_ProjectedPC *ppc;

    PetscFunctionBegin;
    ppc = (PC_ProjectedPC*)pc->data;
    *caln = ppc->caln;
    PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetMaxSize"
/*@
    PCProjectedPCSetMaxSize - Sets the maximum size of the projected out subspaces. This
                              function MUST be called BEFORE PCSetUp, for proper memory
                              allocation.

    Not Collective

    Input Parameters:
    pc   - the preconditioner context
    maxn - the maximum size of the projected out subspace

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetMaxSize(PC pc, PetscInt maxn)
{
    PetscErrorCode ierr, (*f)(PC,PetscInt);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCSetMaxSize_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,maxn);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetCurrentSize"
/*@
    PCProjectedPCSetCurrentSize - Sets the current size of the projected out subspaces. 

    Not Collective

    Input Parameters:
    pc   - the preconditioner context
    curn - the current size of the projected out subspace

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetCurrentSize(PC pc, PetscInt curn)
{
    PetscErrorCode ierr, (*f)(PC,PetscInt);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCSetCurrentSize_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,curn);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetComputedSize"
/*@
    PCProjectedPCSetComputedSize - Sets the computed size of the projected out subspaces. 

    Not Collective

    Input Parameters:
    pc   - the preconditioner context
    caln - the computed size of the projected out subspace

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetComputedSize(PC pc, PetscInt caln)
{
    PetscErrorCode ierr, (*f)(PC,PetscInt);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCSetComputedSize_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,caln);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetPC"
/*@
    PCProjectedPCSetPC - Sets the pc object for this preconditioner.

    Not Collective

    Input Parameters:
    pc    - the preconditioner context
    spc   - the pc object used in this preconditioner

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetPC(PC pc, PC spc)
{
    PetscErrorCode ierr, (*f)(PC,PC);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCSetPC_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,spc);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetSpaces"
/*@
    PCProjectedPCSetSpaces - Sets the vectors of the projected out subspaces. 

    Not Collective

    Input Parameters:
    pc   - the preconditioner context
    Q    - the vectors of the projected out subspace
    Z    - the vectors of the projected out subspace
    KiZ  - the vectors of the projected out subspace
    maxn - size of spaces

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetSpaces(PC pc, Vec* Q, Vec* Z, Vec* KiZ, PetscInt maxn)
{
    PetscErrorCode ierr, (*f)(PC,Vec*,Vec*,Vec*,PetscInt);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCSetSpaces_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,Q,Z,KiZ,maxn);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetQ"
/*@
    PCProjectedPCSetQ - Sets the vectors of the projected out subspaces. 

    Not Collective

    Input Parameters:
    pc   - the preconditioner context
    Q    - the vectors of the projected out subspace

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetQ(PC pc, Vec* Q)
{
    PetscErrorCode ierr, (*f)(PC,Vec*);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCSetQ_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,Q);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetZ"
/*@
    PCProjectedPCSetZ - Sets the vectors of the projected out subspaces. 

    Not Collective

    Input Parameters:
    pc   - the preconditioner context
    Z    - the vectors of the projected out subspace

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetZ(PC pc, Vec* Z)
{
    PetscErrorCode ierr, (*f)(PC,Vec*);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCSetZ_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,Z);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCSetKiZ"
/*@
    PCProjectedPCSetKiZ - Sets the KiZ vectors of the projected out subspaces. 

    Not Collective

    Input Parameters:
    pc    - the preconditioner context
    KiZ   - the vectors of the projected out subspace

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCSetKiZ(PC pc, Vec* KiZ)
{
    PetscErrorCode ierr, (*f)(PC,Vec*);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCSetKiZ_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,KiZ);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCGetMaxSize"
/*@
    PCProjectedPCGetMaxSize - Gets the maximum size of the projected out subspaces. 

    Not Collective

    Input Parameters:
    pc   - the preconditioner context

    Output Parameters:
    maxn - the maximum size of the projected out subspace

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetMaxSize(PC pc, PetscInt* maxn)
{
    PetscErrorCode ierr, (*f)(PC,PetscInt*);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCGetMaxSize_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,maxn);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCGetCurrentSize"
/*@
    PCProjectedPCGetMaxSize - Gets the current size of the projected out subspaces. 

    Not Collective

    Input Parameters:
    pc   - the preconditioner context

    Output Parameters:
    curn - the current size of the projected out subspace

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetCurrentSize(PC pc, PetscInt* curn)
{
    PetscErrorCode ierr, (*f)(PC,PetscInt*);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCGetCurrentSize_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,curn);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCProjectedPCGetComputedSize"
/*@
    PCProjectedPCGetComputedSize - Gets the computed size of the projected out subspaces. 

    Not Collective

    Input Parameters:
    pc   - the preconditioner context

    Output Parameters:
    caln - the computed size of the projected out subspace

    Options Database Key:

    Levels: beginner

    Concepts: PCProjectedPC preconditioner

.seealso:
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCProjectedPCGetComputedSize(PC pc, PetscInt* caln)
{
    PetscErrorCode ierr, (*f)(PC,PetscInt*);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCProjectedPCGetComputedSize_C",(void (**)(void))&f);CHKERRQ(ierr);
    if (f) {
      ierr = (*f)(pc,caln);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} 

#undef __FUNCT__
#define __FUNCT__ "PCRegisterProjectedPC"
/*@
    PCRegisterProjectedPC - Registers the ProjectedPC

    Not Collective

    Input Parameters:

    Options Database Key:

    Levels: beginner

    Concepts: ProjectedPC preconditioner

.seealso: 
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCRegisterProjectedPC()
{
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = PCRegister("projectedpc",0,"PCCreate_ProjectedPC",PCCreate_ProjectedPC);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
