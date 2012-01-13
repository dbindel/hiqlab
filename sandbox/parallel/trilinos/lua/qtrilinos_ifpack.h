#include "Teuchos_ParameterList.hpp"
#include "Epetra_CrsMatrix.h"

#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"

Ifpack_Preconditioner* qCreate_Ifpack_Preconditioner
    (Epetra_CrsMatrix* SM, const string& name,
     Teuchos::ParameterList* IfpackList, int overlap);
