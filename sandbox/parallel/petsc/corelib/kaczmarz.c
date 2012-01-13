#define PETSCKSP_DLL

/* -------------------------------------------------------------------
  This file implements a Kaczmarz preconditioner in PETSc as part of PC
   
  The following basic routines are required for each preconditioner.

    PCCreate_Kaczmarz()        - Creates a preconditioner context
    PCSetFromOptions_Kaczmarz   - Sets runtimes options
    PCApply_Kaczmarz            - Applies the preconditioner context
    PCDestroy_Kaczmarz          - Destroys the preconditioner context

  These routines are actually called via the common user interface 
  routines 

    PCCreate()
    PCSetFromOptions()
    PCApply()
    PCDestroy()

  so the application code interface remains identical for all 
  preconditioners.

  Another key routine is:

    PCSetUp_Kaczmarz()          - Prepares for the use of the preconditioner
                                  by setting data structures and options. The
                                  interface routine PCSetUp() is not usually
                                  called directly by the user, but instead is 
                                  called by PCApply() if necessary.

  Additional basic routines are:

    PCView_Kaczmarz()           - Prints details of runtime options that have
                                  actually been used. This is called via the
                                  interface routine PCView()

  -------------------------------------------------------------------- */

/*
    Include files needed for the Kaczmarz preconditioner:

      pcimpl.h - private include file intended for use by all preconditioners
*/

#include "private/pcimpl.h"   /*I "petscpc.h"  I*/
#include "private/matimpl.h"
#include "src/mat/impls/aij/seq/aij.h"
/*
    Private context (data structure) for the Kaczmarz preconditioner.
*/
typedef struct {
    Vec        diagAAt;     /* vector containing the reciprocals of the 
                             diagonal elements of A*A' */
    Vec        nnzcol;      /* vector containing number of nonzeros per
                             column */
    Vec        lvec;        /* vector to use in scatter context */
    VecScatter vctx;        /* vector scatter context */  
    Mat        Aloc;        /* local portion of matrix in Mat_SeqAIJ format */
    PetscInt   git;         /* Number of global iterations */
    PetscInt   lit;         /* Number of local  iterations */
    PetscInt   nghost;      /* Number of ghost nodes required */
    PetscInt*  garray;      /* Global ids of ghost nodes */
    PetscReal  omega;       /* relaxation parameter */
} PC_Kaczmarz;


EXTERN_C_BEGIN
EXTERN PetscErrorCode PCCreate_Kaczmarz(PC);
EXTERN PetscErrorCode PCKaczmarzSetIterations_Kaczmarz(PC,PetscInt,PetscInt);
EXTERN PetscErrorCode PCKaczmarzSetOmega_Kaczmarz(PC,PetscReal);
EXTERN_C_END
static PetscErrorCode PCDestroy_Kaczmarz(PC);
static PetscErrorCode PCSetUp_Kaczmarz(PC);
static PetscErrorCode PCApply_Kaczmarz(PC,Vec,Vec);
EXTERN PetscErrorCode PCKaczmarzSetIterations(PC,PetscInt,PetscInt);

/* ------------------------------------------------------------------- */
/*
    PCCreate_Kaczmarz - Creates a Kaczmarz preconditioner context, PC_Kaczmarz,
    and sets this as the private data within the generic preconditioning context, 
    PC, that was created within PCCreate().

    Input Parameter:
    pc - the preconditioner context

    Application Interface Routine: PCCreate()
*/

/*MC
      PCKACZMARZ - Kaczmarz

    Options Database Key:

    Levels: beginner

    Concepts: Kaczmarz, preconditioners

    Notes:

.seealso: PCCreate(), PCSetType(), PCType() (for list of available types), PC,

M*/

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCCreate_Kaczmarz"
PetscErrorCode PETSCKSP_DLLEXPORT PCCreate_Kaczmarz(PC pc)
{
    PC_Kaczmarz     *kz;
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr     = PetscNew(PC_Kaczmarz,&kz);CHKERRQ(ierr);
    pc->data = (void*)kz;

    /*
      Logs the memory usage;
    */
    ierr = PetscLogObjectMemory(pc,sizeof(PC_Kaczmarz));CHKERRQ(ierr);

    /*
      Initialize the pointers to vectors to ZERO;
    */
    kz->diagAAt = 0;
    kz->nnzcol  = 0;
    kz->Aloc    = 0;
    kz->lvec    = 0;
    kz->garray  = 0;
    kz->vctx    = 0;

    kz->git     = 1;
    kz->lit     = 1;
    kz->nghost  = 0;
    kz->omega   = 1;

    /*
      Set the pointers for the functions that are provided above.
      Now when the user-level routines (sucha as PCApply(), PCDestroy(), etc.)
      are called, they will automatically call these functions. Note we choose not
      to provide a couple of these functions since they are not needed.
    */
    pc->ops->apply              = PCApply_Kaczmarz;
    pc->ops->applytranspose     = PCApply_Kaczmarz;
    pc->ops->setup              = PCSetUp_Kaczmarz;
    pc->ops->destroy            = PCDestroy_Kaczmarz;
    pc->ops->setfromoptions     = 0;//PCSetFromOptions_Kaczmarz();
    pc->ops->view               = 0;//PCView_Kaczmarz();
    pc->ops->applyrichardson    = 0;
    pc->ops->applysymmetricleft = 0;
    pc->ops->applysymmetricright= 0;

    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc, "PCKaczmarzSetIterations_C",
                                                              "PCKaczmarzSetIterations_Kaczmarz",
                                                               PCKaczmarzSetIterations_Kaczmarz);
                                                               CHKERRQ(ierr);

    ierr = PetscObjectComposeFunctionDynamic((PetscObject)pc, "PCKaczmarzSetOmega_C",
                                                              "PCKaczmarzSetOmega_Kaczmarz",
                                                               PCKaczmarzSetOmega_Kaczmarz);
                                                               CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
EXTERN_C_END

/* ------------------------------------------------------------------- */
/*
    PCDestroy_Kaczmarz - Destroys the private context for the Kaczmarz preconditioner
    that was created with PCCreate_Kaczmarz().

    Input Parameter:
    pc - the preconditioner context

    Application Interface Routine: PCDestroy()
*/
#undef __FUNCT__
#define __FUNCT__ "PCDestroy_Kaczmarz"
static PetscErrorCode PCDestroy_Kaczmarz(PC pc)
{
    PC_Kaczmarz    *kz = (PC_Kaczmarz*)pc->data;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    if (kz->diagAAt)    {ierr = VecDestroy(kz->diagAAt); CHKERRQ(ierr);}
    if (kz->nnzcol)     {ierr = VecDestroy(kz->nnzcol);  CHKERRQ(ierr);}
    if (kz->Aloc)       {ierr = MatDestroy(kz->Aloc); CHKERRQ(ierr);}
    if (kz->lvec)       {ierr = VecDestroy(kz->lvec); CHKERRQ(ierr);}
    if (kz->vctx)       {ierr = VecScatterDestroy(kz->vctx); CHKERRQ(ierr);}
    /*
        Free the private data structure that was hanging off the PC
    */
    ierr = PetscFree(kz);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
    PCSetUp_Kaczmarz - Prepares for the use of the Kaczmarz preconditioner
                       by setting data structures and options.

    Input Parameter:
    pc - the preconditioner context

    Application Interface Routine: PCSetUp()

    Notes:
    The interface routine PCSetUp() is not called directly by the user,
    but instead is called by PCApply() if necessary.
*/
#undef __FUNCT__
#define __FUNCT__ "PCSetUp_Kaczmarz"
static PetscErrorCode PCSetUp_Kaczmarz(PC pc)
{
    PC_Kaczmarz     *kz = (PC_Kaczmarz*)pc->data;
    Vec             diagAAt, nnzcol;
    PetscErrorCode  ierr;
    PetscInt        nghost = 0;
    PetscInt*       garray;
    PetscInt        i,j;
    PetscInt        Msize, Nsize;
    PetscInt        mymin, mymax1;
    PetscInt        irow, ncols;
    const PetscInt*       cols;
    const PetscScalar*    vals;
    PetscScalar     sqval;
    PetscScalar*    sumsquare;
    PetscScalar*    ones;
    PetscInt*       indices;
    IS              ISr,ISc;

    PetscFunctionBegin;

    /* 
        Allocate space first time the function is called
    */
    if (pc->setupcalled == 0) {
        ierr = MatGetVecs(pc->pmat,&kz->diagAAt,PETSC_NULL);CHKERRQ(ierr);
        ierr = MatGetVecs(pc->pmat,&kz->nnzcol ,PETSC_NULL);CHKERRQ(ierr);
        PetscLogObjectParent(pc,kz->diagAAt);
        PetscLogObjectParent(pc,kz->nnzcol );
    }

    diagAAt = kz->diagAAt;
    nnzcol  = kz->nnzcol;
    ierr = VecSet(nnzcol,0);CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(pc->pmat,&mymin,&mymax1);CHKERRQ(ierr);
    ierr = MatGetSize(pc->pmat,&Msize,&Nsize);           CHKERRQ(ierr);
    ierr = VecGetArray(diagAAt,&sumsquare);              CHKERRQ(ierr);


    if (diagAAt && nnzcol) {

        // --  Make an array as long as the number of columns 
        //     and mark those columns that are off the diagonal portion
        ierr  = PetscMalloc(Nsize*sizeof(PetscInt),&indices);CHKERRQ(ierr);
        ierr  = PetscMemzero(indices,Nsize*sizeof(PetscInt));CHKERRQ(ierr);

        for (i = mymin; i < mymax1; ++i) {

            ierr = MatGetRow(pc->pmat,i,&ncols,&cols,&vals);CHKERRQ(ierr);

            sqval = 0;
            for (j = 0; j < ncols; ++j) {
                sqval += PetscConj(vals[j])*vals[j];

                // -- If off processor record
                if (mymin > cols[j] || mymax1 <= cols[j]) {
                    if (!indices[cols[j]]) nghost++;
                    indices[cols[j]] = 1;
                }
            }
            if (PetscRealPart(sqval)==0 && PetscImaginaryPart(sqval)==0){
                ierr = PetscInfo(pc,"All zero row detected in matrix, assuming 1 on the diagonal at these rows\n");
                       CHKERRQ(ierr);
                sqval = 1;
            }
            sumsquare[i-mymin] = 1/PetscRealPart(sqval);

            ierr = MatRestoreRow(pc->pmat,i,&ncols,&cols,&vals);CHKERRQ(ierr);
        }

        // -- Set indices of ghost nodes
        ierr  = PetscMalloc(nghost*sizeof(PetscInt),&garray);CHKERRQ(ierr);
        j = 0;
        for (i = 0; i < Nsize; ++i)
            if (indices[i]) garray[j++] = i;

        if (kz->garray){
            ierr = PetscFree(kz->garray);CHKERRQ(ierr);
        }
        kz->nghost = nghost;
        kz->garray = garray;

        ierr = PetscFree(indices);CHKERRQ(ierr);

        // -- Construct index array and ones array to set in nnzcol vector
        ierr = PetscMalloc((mymax1-mymin+nghost)*sizeof(PetscScalar),&ones);CHKERRQ(ierr);
        ierr = PetscMalloc((mymax1-mymin+nghost)*sizeof(PetscInt),&indices);CHKERRQ(ierr);
        for (i = 0; i < mymax1-mymin; ++i) {
            ones[i] = 1;
            indices[i] = mymin + i;
        }
        for (i = 0; i < nghost; ++i) {
            ones[i+mymax1-mymin]    = 1;
            indices[i+mymax1-mymin] = garray[i];
        }
        ierr = VecSetValues(nnzcol,mymax1-mymin+nghost,indices,ones,ADD_VALUES);
        ierr = VecAssemblyBegin(nnzcol);
        ierr = VecAssemblyEnd(nnzcol);

        // -- Clean up
        ierr = PetscFree(indices);CHKERRQ(ierr);
        ierr = PetscFree(ones);CHKERRQ(ierr);

    }

    ierr = VecRestoreArray(diagAAt,&sumsquare);           CHKERRQ(ierr);

    // -- Allocate space for indices below
    ierr = PetscMalloc((mymax1-mymin+nghost)*sizeof(PetscInt),&indices);CHKERRQ(ierr);

    // -- Construct local submatrix
    if (kz->Aloc) {
        ierr = MatDestroy(kz->Aloc);CHKERRQ(ierr);
    }
    for (i = 0; i < mymax1-mymin; ++i)
        indices[i] = mymin + i;
    for (i = 0; i < nghost; ++i)
        indices[i+mymax1-mymin] = garray[i];
    ierr = ISCreateGeneral(pc->pmat->comm,mymax1-mymin+nghost,indices,&ISc);CHKERRQ(ierr);
    ierr = ISCreateStride (pc->pmat->comm,mymax1-mymin       ,mymin,1,&ISr);CHKERRQ(ierr);
    ierr = MatGetLocalMatCondensed(pc->pmat,MAT_INITIAL_MATRIX,&ISr,&ISc,&kz->Aloc);CHKERRQ(ierr);
    ierr = ISDestroy(ISr);CHKERRQ(ierr);   
    ierr = ISDestroy(ISc);CHKERRQ(ierr);   
//    ierr = MatGetLocalMat(pc->pmat,MAT_INITIAL_MATRIX,&kz->Aloc);CHKERRQ(ierr);

    // -- Construct local submatrix and scatter context
    if (kz->lvec) {
        ierr = VecDestroy(kz->lvec);CHKERRQ(ierr);
    }
    if (kz->vctx) {
        ierr = VecScatterDestroy(kz->vctx);CHKERRQ(ierr);
    }
    ierr = VecCreateSeq(PETSC_COMM_SELF,mymax1-mymin+nghost,&kz->lvec);CHKERRQ(ierr);
    ierr = ISCreateGeneral(pc->pmat->comm ,mymax1-mymin+nghost,indices,&ISc);CHKERRQ(ierr);
    ierr = ISCreateStride (PETSC_COMM_SELF,mymax1-mymin+nghost,0,1,&ISr);CHKERRQ(ierr);
    ierr = VecScatterCreate(nnzcol,ISc,kz->lvec,ISr,&kz->vctx);CHKERRQ(ierr);
    ierr = ISDestroy(ISr);CHKERRQ(ierr);   
    ierr = ISDestroy(ISc);CHKERRQ(ierr);   

    // -- Clean up
    ierr = PetscFree(indices);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCRegisterKaczmarz"
/*@
    PCRegisterKaczmarz - Registers the Kaczmarz 

    Not Collective

    Input Parameters:

    Options Database Key:

    Levels: beginner

    Concepts: Kaczmarz preconditioner

.seealso: 
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCRegisterKaczmarz()
{
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = PCRegister("kaczmarz",0,"PCCreate_Kaczmarz",PCCreate_Kaczmarz);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCKaczmarzSetIterations_Kaczmarz"
PetscErrorCode PETSCKSP_DLLEXPORT PCKaczmarzSetIterations_Kaczmarz(PC pc, PetscInt git, PetscInt lit)
{
    PC_Kaczmarz *kz;

    PetscFunctionBegin;

    kz  = (PC_Kaczmarz*)pc->data;
    if(git) kz->git = git;
    if(lit) kz->lit = lit;
    PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__
#define __FUNCT__ "PCKaczmarzSetIterations"
/*@
    PCKaczmarzSetIterations - Sets the number of local intra-processor iterations
    and global inter-processor iterations. 

    Not Collective

    Input Parameters:
    pc  - the preconditioner context
    git - number of global iterations
    lit - number of local iterations

    Options Database Key:

    Levels: beginner

    Concepts: Kaczmarz preconditioner

.seealso: 
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCKaczmarzSetIterations(PC pc, PetscInt git, PetscInt lit)
{
    PetscErrorCode ierr, (*f)(PC,PetscInt,PetscInt);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCKaczmarzSetIterations_C",(void (**)(void))&f); CHKERRQ(ierr);
    if (f) {
        ierr = (*f)(pc,git,lit);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "PCKaczmarzSetOmega_Kaczmarz"
PetscErrorCode PETSCKSP_DLLEXPORT PCKaczmarzSetOmega_Kaczmarz(PC pc, PetscReal omega)
{
    PC_Kaczmarz *kz;

    PetscFunctionBegin;

    kz  = (PC_Kaczmarz*)pc->data;
    kz->omega = omega;
    PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__
#define __FUNCT__ "PCKaczmarzOmega"
/*@
    PCKaczmarzSetOmega - Sets omega the damping parameter in iteration

    Not Collective

    Input Parameters:
    pc    - the preconditioner context
    omega - the damping parameter

    Options Database Key:

    Levels: beginner

    Concepts: Kaczmarz preconditioner

.seealso: 
@*/
PetscErrorCode PETSCKSP_DLLEXPORT PCKaczmarzSetOmega(PC pc, PetscReal omega)
{
    PetscErrorCode ierr, (*f)(PC,PetscReal);

    PetscFunctionBegin;
    PetscValidHeaderSpecific(pc,PC_COOKIE,1);
    ierr = PetscObjectQueryFunction((PetscObject)pc,"PCKaczmarzSetOmega_C",(void (**)(void))&f); CHKERRQ(ierr);
    if (f) {
        ierr = (*f)(pc,omega);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
    PCApply_Kaczmarz - Applies Kaczmarz preconditioner to a vector.

    Input Parameters:
    pc - the preconditioner context
    x  - input vector

    Output Parameter:
    y  - output vector

    Application Interface Routine: PCApply()
*/
#undef __FUNCT__
#define __FUNCT__ "PCApply_Kaczmarz"
static PetscErrorCode PCApply_Kaczmarz(PC pc, Vec x, Vec y)
{
    PC_Kaczmarz    *kz = (PC_Kaczmarz*)pc->data;
    Mat_SeqAIJ*    Aloc=  (Mat_SeqAIJ*)(kz->Aloc->data);
    PetscErrorCode ierr;
    PetscInt       il,ig,ir,ic;
    PetscInt       numrows,mymin,mymax1;
    PetscInt       col;
    PetscScalar*   yvals;
    PetscScalar*   xvals;
    PetscScalar*   diagAAt;
    PetscScalar    ztemp;
    PetscScalar    aval;
    PetscReal      omega = kz->omega;

    PetscFunctionBegin;

    // -- Get data
    ierr = MatGetOwnershipRange(pc->pmat,&mymin,&mymax1);
    numrows = mymax1 - mymin;
    
    // -- Store right hand side locally
    ierr = VecGetArray(x,&xvals); CHKERRQ(ierr);

    // -- Get diagAAt inverse
    ierr = VecGetArray(kz->diagAAt,&diagAAt); CHKERRQ(ierr);

    // -- Initialize initial solution to zero
    ierr = VecSet(y,0);CHKERRQ(ierr);

    // -- Global iteration
    for (ig = 0; ig < kz->git; ++ig) { 

        // -- Get global solution
        ierr = VecScatterBegin(kz->vctx,y,kz->lvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
//        ierr = VecScatterBegin(y,kz->lvec,INSERT_VALUES,SCATTER_FORWARD,kz->vctx);CHKERRQ(ierr);//2.3.2
        ierr = VecScatterEnd(kz->vctx,y,kz->lvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
//        ierr = VecScatterEnd(y,kz->lvec,INSERT_VALUES,SCATTER_FORWARD,kz->vctx);CHKERRQ(ierr); //2.3.2
        ierr = VecGetArray(kz->lvec,&yvals); CHKERRQ(ierr);

        // -- Processor iteration
        for (il = 0; il < kz->lit; ++il) {

            for (ir = 0; ir < numrows; ++ir) {
                ztemp = 0;
                for (ic = 0; ic < Aloc->i[ir+1]-Aloc->i[ir]; ++ic) {
                    col = Aloc->j[Aloc->i[ir]+ic];
                    aval=Aloc->a[Aloc->i[ir]+ic];
                    ztemp += yvals[col]*aval;
                }
                for (ic = 0; ic < Aloc->i[ir+1]-Aloc->i[ir]; ++ic) {
                    col = Aloc->j[Aloc->i[ir]+ic];
                    aval= PetscConj(Aloc->a[Aloc->i[ir]+ic]);
//                    aval= Aloc->a[Aloc->i[ir]+ic];
                    yvals[col] += aval * (xvals[ir] - ztemp)*diagAAt[ir]*omega;
                }
            }
        }

        // -- Return solution and average
        ierr = VecRestoreArray(kz->lvec,&yvals); CHKERRQ(ierr);
        ierr = VecSet(y,0); CHKERRQ(ierr);
        ierr = VecScatterBegin(kz->vctx,kz->lvec,y,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
//        ierr = VecScatterBegin(kz->lvec,y,ADD_VALUES,SCATTER_REVERSE,kz->vctx);CHKERRQ(ierr);
        ierr = VecScatterEnd(kz->vctx,kz->lvec,y,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
//        ierr = VecScatterEnd(kz->lvec,y,ADD_VALUES,SCATTER_REVERSE,kz->vctx);CHKERRQ(ierr);
        ierr = VecPointwiseDivide(y,y,kz->nnzcol);CHKERRQ(ierr);

    }

    ierr = VecRestoreArray(x,&xvals); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
