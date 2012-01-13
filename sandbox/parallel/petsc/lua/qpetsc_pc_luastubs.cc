#include "qpetsc_pc_luastubs.h"

int qPCMGSetType(PC pc, int form)
{
    PCMGType pc_types[] = {
        PC_MG_MULTIPLICATIVE,
        PC_MG_ADDITIVE,
        PC_MG_FULL,
        PC_MG_KASKADE
    };

    PetscErrorCode ierr;
    ierr = PCMGSetType(pc,pc_types[form]);CHKERRQ(ierr);

    return 0;
}

int qPCMGSetCycleType(PC pc, int ctype)
{
    PCMGCycleType cycle_types[] = {
        PC_MG_CYCLE_V,
        PC_MG_CYCLE_W
    };

    PetscErrorCode ierr;
    ierr = PCMGSetCycleType(pc,cycle_types[ctype]);CHKERRQ(ierr);

    return 0;
}

Mat PCMGGetOperator(PC pc, int igrid)
{
    PetscErrorCode ierr;
    Mat A;

    ierr = PCMGGetOperator(pc,igrid,&A);

    return A;
}
