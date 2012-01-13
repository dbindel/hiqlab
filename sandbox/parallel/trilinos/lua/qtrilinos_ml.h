#include "Teuchos_ParameterList.hpp"
#include "Epetra_CrsMatrix.h"

#include "ml_include.h"
#include "ml_MultiLevelOperator.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

void ML_Epetra_SetDefaults(Teuchos::ParameterList* MLList, const string& name);


ML_Epetra::MultiLevelPreconditioner* qCreate_ML_Epetra_MultiLevelPreconditioner
                       (Epetra_CrsMatrix* SM, Teuchos::ParameterList* MLList);
