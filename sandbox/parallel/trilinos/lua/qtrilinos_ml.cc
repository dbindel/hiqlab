#include "qtrilinos_ml.h"
#include <iostream>


void ML_Epetra_SetDefaults(Teuchos::ParameterList* MLList, const string& name)
{
     ML_Epetra::SetDefaults(name , *MLList);
}


ML_Epetra::MultiLevelPreconditioner* qCreate_ML_Epetra_MultiLevelPreconditioner
                       (Epetra_CrsMatrix* SM, Teuchos::ParameterList* MLList)
{
    ML_Epetra::MultiLevelPreconditioner* MLPrec =
        new ML_Epetra::MultiLevelPreconditioner(*SM, *MLList);

    return MLPrec;
}
