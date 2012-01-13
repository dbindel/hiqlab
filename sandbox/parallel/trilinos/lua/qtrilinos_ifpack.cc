#include "qtrilinos_ifpack.h"

Ifpack_Preconditioner* qCreate_Ifpack_Preconditioner
    (Epetra_CrsMatrix* SM, const string& name,
     Teuchos::ParameterList* IfpackList, int overlap)
{
    Ifpack Factory;
    Ifpack_Preconditioner* Prec;

    // -- Create preconditioner and set parameters
    Prec = Factory.Create(name, SM, overlap);
    Prec->SetParameters(*IfpackList);

    // -- Set up initialization
    Prec->Initialize();
    Prec->Compute();

    return Prec;
}
