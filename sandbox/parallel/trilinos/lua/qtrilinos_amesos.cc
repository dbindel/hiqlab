#include "qtrilinos_amesos.h"

Amesos_BaseSolver* qAmesos_BaseSolverCreate(Epetra_LinearProblem* ELP, int s_type)
{
    string solver[7] = { "Lapack", "Klu", "Umfpack", "Superlu", "Superludist", "Scalapack", "Mumps" };

    Amesos_BaseSolver*  ABS;
    Amesos Amesos_Factory;

    ABS = Amesos_Factory.Create(solver[s_type], *ELP);
    assert(ABS);

    return ABS;
}

Amesos_BaseSolver_Complex* qAmesos_BaseSolver_ComplexCreate
                     (Epetra_LinearProblem_Complex* ELP, int s_type)
{
    string solver[7] = { "Lapack", "Klu", "Umfpack", "Superlu", "Superludist", "Scalapack", "Mumps" };

    Amesos_BaseSolver_Complex*  ABS;
    Amesos_Complex Amesos_Factory;

    ABS = Amesos_Factory.Create(solver[s_type], *ELP);
    assert(ABS);

    return ABS;
}
