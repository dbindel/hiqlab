#define PETSCKSP_DLL

#include "private/matimpl.h"
//#include "src/mat/matimpl.h" //2.3.2
#include "matshellab.h"
#include "petscmat.h"

/*
    MatShellABCreate - Creates a MatShellAB context.

    Output Parameter:
    abshell - MatShellAB context.
*/

#undef __FUNCT__
#define __FUNCT__ "MatShellABCreate"
PetscErrorCode PETSCKSP_DLLEXPORT MatShellABCreate(MatShellAB** matshell)
{
    MatShellAB* newctx;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = PetscNew(MatShellAB,&newctx);CHKERRQ(ierr);
    newctx->A     = 0;
    newctx->B     = 0;
    newctx->v     = 0;
    newctx->alpha = 0;
    newctx->beta  = 0;
    *matshell     = newctx;

    PetscFunctionReturn(0);
}

/*
    MatShellABSetOperators - Set operators A and B of the shell matrix

    Input Parameters:
    A - Operator A
    B - Operator B

    Output Parameter:
    shell - 
*/
#undef __FUNCT__
#define __FUNCT__ "MatShellABSetOperators"
PetscErrorCode PETSCKSP_DLLEXPORT MatShellABSetOperators(MatShellAB* shell, Mat A, Mat B)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;
    if (shell->A) { ierr = MatDestroy(shell->A);CHKERRQ(ierr);}
    shell->A = A;
    if (shell->B) { ierr = MatDestroy(shell->B);CHKERRQ(ierr);}
    shell->B = B;
    if (shell->v) { ierr = VecDestroy(shell->v);CHKERRQ(ierr);}
    ierr = MatGetVecs(shell->A,&shell->v,PETSC_NULL);CHKERRQ(ierr);

    ierr = PetscObjectReference((PetscObject)A);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)B);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/*
    MatShellABSetShift - Set shift pair of the shell matrix.

    Input Parameter:
    alpha - 
    beta  -

    Output Parameter:
    shell -
*/
#undef __FUNCT__
#define __FUNCT__ "MatShellABSetShift"
PetscErrorCode PETSCKSP_DLLEXPORT MatShellABSetShift(MatShellAB* shell, PetscScalar alpha, PetscScalar beta)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;
    shell->alpha = alpha;
    shell->beta  = beta;
    PetscFunctionReturn(0);
}

/*
    MatShellABApply - Apply this operator

    Input Parameters:
    op - A MatShell
    x  - input vector
    
    Output Parameters:
    y  - output vector
*/
#undef __FUNCT__
#define __FUNCT__ "MatShellABApply"
PetscErrorCode PETSCKSP_DLLEXPORT MatShellABApply(Mat op, Vec x, Vec y)
{
    PetscErrorCode ierr;
//    Mat_Shell   *shell  = (Mat_Shell*)op->data;
    void*       ctx;
    MatShellAB *shellab;
    Mat A,B;
    Vec v;
    PetscScalar alpha,beta;

    ierr = MatShellGetContext(op,&ctx);CHKERRQ(ierr);
    shellab= (MatShellAB*)ctx;
    A      = shellab->A;
    B      = shellab->B;
    v      = shellab->v;
    alpha  = shellab->alpha;
    beta   = shellab->beta;

    PetscFunctionBegin;

    // -- v = A*x
    ierr = MatMult(A,x,v);CHKERRQ(ierr);

    // -- y = B*x
    ierr = MatMult(B,x,y);CHKERRQ(ierr);

    // -- y = beta*v - alpha*y
    ierr = VecAXPBY(y, beta, -alpha,v);

    PetscFunctionReturn(0);
}

/*
    MatShellABApplyTranspose - Apply the transpose of this operator

    Input Parameters:
    op - A MatShell
    x  - input vector
    
    Output Parameters:
    y  - output vector
*/
#undef __FUNCT__
#define __FUNCT__ "MatShellABApplyTranspose"
PetscErrorCode PETSCKSP_DLLEXPORT MatShellABApplyTranspose(Mat op, Vec x, Vec y)
{
    PetscErrorCode ierr;
//    Mat_Shell   *shell  = (Mat_Shell*)op->data;
    void*       ctx;
    MatShellAB *shellab;
    Mat A,B;
    Vec v;
    PetscScalar alpha,beta;

    ierr = MatShellGetContext(op,&ctx);CHKERRQ(ierr);
    shellab= (MatShellAB*)ctx;
    A      = shellab->A;
    B      = shellab->B;
    v      = shellab->v;
    alpha  = shellab->alpha;
    beta   = shellab->beta;

    PetscFunctionBegin;

    // -- v = At*x
    ierr = MatMultTranspose(A,x,v);CHKERRQ(ierr);

    // -- y = Bt*conj(x)
    ierr = MatMultTranspose(B,x,y);CHKERRQ(ierr);

    // -- y = beta*v - alpha*y
    ierr = VecAXPBY(y, beta, -alpha,v);

    PetscFunctionReturn(0);
}

/*
    MatShellABDestroy - Destroy this operator

    Input Parameters:
    shell - MatABShell context
*/
#undef __FUNCT__
#define __FUNCT__ "MatShellABDestroy"
PetscErrorCode PETSCKSP_DLLEXPORT MatShellABDestroy(MatShellAB* shell)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    if (shell->A) { ierr = MatDestroy(shell->A);CHKERRQ(ierr);}
    if (shell->B) { ierr = MatDestroy(shell->B);CHKERRQ(ierr);}
    if (shell->v) { ierr = VecDestroy(shell->v);CHKERRQ(ierr);}

    ierr = PetscFree(shell);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
