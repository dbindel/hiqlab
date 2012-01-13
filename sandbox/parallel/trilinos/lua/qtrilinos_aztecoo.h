#ifndef _QTRILINOS_AZTECOO_H
#define _QTRILINOS_AZTECOO_H

#include "AztecOO.h"
#include "Epetra_LinearProblem.h"

#include "ml_include.h"
#include "ml_MultiLevelOperator.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

#include "Ifpack_Preconditioner.h"

#include "Teuchos_ParameterList.hpp"

// -- Constructor for AztecOO
AztecOO* qAztecOOCreate(Epetra_LinearProblem* ELP);

int* qConstructAztecOptionsArray(double* op);

// -- Aztec Options and Params object
class AztecOptionsParams
{
  public:
    AztecOptionsParams();
    virtual ~AztecOptionsParams();

    // --  Set options and params
    void SetOption(int option, int value) {options[option]=value;};
    void SetParam(int param, double value){params  [param]=value;};

    // --  Get options and params
    int     GetOption(int option) {return options[option];};
    int*    GetOptions() {return options;};
    double  GetParam(int param)   {return params[param];};
    double* GetParams() {return params;};

    void    SetPL(Teuchos::ParameterList* PL);

  private:

    int* options;
    double* params;

};

#endif /* _QTRILINOS_AZTECOO_H */
