#include <iostream>
#include "qtrilinos_aztecoo.h"

AztecOO* qAztecOOCreate(Epetra_LinearProblem* ELP)
{
    AztecOO* AztecOOSolver;
    AztecOOSolver = new AztecOO(*ELP);

    return AztecOOSolver;
}

AztecOptionsParams::AztecOptionsParams()
{
    options = new int[AZ_OPTIONS_SIZE];
    params  = new double[AZ_OPTIONS_SIZE];
    AZ_defaults(options,params);
}

AztecOptionsParams::~AztecOptionsParams()
{
    delete[] options;
    delete[] params;
}

void AztecOptionsParams::SetPL(Teuchos::ParameterList* PL)
{
    PL->set("smoother: Aztec options", options);
    PL->set("smoother: Aztec params" , params );
}
