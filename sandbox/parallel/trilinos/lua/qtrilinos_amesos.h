#ifndef _QTRILINOS_AMESOS_H
#define _QTRILINOS_AMESOS_H

#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "trilinos_epetra_linearproblem.h"
#include "trilinos_amesos_base.h"

// -- Constructor for Amesos_BaseSolver
Amesos_BaseSolver* qAmesos_BaseSolverCreate
                         (Epetra_LinearProblem* ELP, int s_type);

// -- Constructor for Amesos_BaseSolver_Complex
Amesos_BaseSolver_Complex* qAmesos_BaseSolver_ComplexCreate
                         (Epetra_LinearProblem_Complex* ELP, int s_type);

#endif /* _QTRILINOS_AMESOS_H */
