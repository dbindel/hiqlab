#include "qtrilinos_teuchos.h"
#include <iostream>

Teuchos::RefCountPtr<Epetra_CrsMatrix>   qCreateRCP_Epetra_CrsMatrix(Epetra_CrsMatrix* ECM)
{
    return Teuchos::rcp(ECM);
}

Teuchos::RefCountPtr<Epetra_MultiVector> qCreateRCP_Epetra_MultiVector(Epetra_MultiVector* EMV)
{
    return Teuchos::rcp(EMV);
}

Teuchos::RefCountPtr<Epetra_Vector>      qCreateRCP_Epetra_Vector(Epetra_Vector* EV)
{
    return Teuchos::rcp(EV);
}

Teuchos::ParameterList* qCreateParameterList()
{
    Teuchos::ParameterList* PL = new Teuchos::ParameterList();
    return PL;
}

void qParameterList_set_double(Teuchos::ParameterList* PL, const string& name, double param)
{
    PL->set(name, param);
}

double qParameterList_get_double(Teuchos::ParameterList* PL, const string& name)
{
    return PL->get<double>(name);
}

void qParameterList_set_doubles(Teuchos::ParameterList* PL, const string& name, double* param)
{
    PL->set(name, param);
}

double* qParameterList_get_doubles(Teuchos::ParameterList* PL, const string& name)
{
    return PL->get<double *>(name);
}

void qParameterList_set_int(Teuchos::ParameterList* PL, const string& name, int param)
{
    PL->set(name, param);
}

int qParameterList_get_int(Teuchos::ParameterList* PL, const string& name)
{
    return PL->get<int>(name);
}

void qParameterList_set_ints(Teuchos::ParameterList* PL, const string& name, int* param)
{
    PL->set(name, param);
}

int* qParameterList_get_ints(Teuchos::ParameterList* PL, const string& name)
{
    return PL->get<int *>(name);
}

void   qParameterList_set_string(Teuchos::ParameterList* PL, const string& name, const string& param)
{
    PL->set(name, param);
}

string qParameterList_get_string(Teuchos::ParameterList* PL, const string& name)
{
    return PL->get<string>(name);
}
