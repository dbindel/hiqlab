#ifndef _QQPETSC_PC_H
#define _QQPETSC_PC_H

#include "petscmg.h"
#include "qpetscmg.h"

int qPCMGSetType(PC pc, int form);
int qPCMGSetCycleType(PC pc, int ctype);
Mat PCMGGetOperator(PC pc, int igrid);

#endif /* _QQPETSC_PC_H */
