#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

Teuchos::RefCountPtr<Epetra_CrsMatrix>   qCreateRCP_Epetra_CrsMatrix(Epetra_CrsMatrix* ECM);
Teuchos::RefCountPtr<Epetra_MultiVector> qCreateRCP_Epetra_MultiVector(Epetra_MultiVector* EMV);
Teuchos::RefCountPtr<Epetra_Vector>      qCreateRCP_Epetra_Vector(Epetra_Vector* EV);

Teuchos::ParameterList*                  qCreateParameterList();
void   qParameterList_set_double(Teuchos::ParameterList* PL, const string& name, double param);
void   qParameterList_set_int   (Teuchos::ParameterList* PL, const string& name, int    param);
void   qParameterList_set_string(Teuchos::ParameterList* PL, const string& name, const string& param);
void   qParameterList_set_ints(Teuchos::ParameterList* PL, const string& name, int* param);
void   qParameterList_set_doubles(Teuchos::ParameterList* PL, const string& name, double* param);
double qParameterList_get_double(Teuchos::ParameterList* PL, const string& name);
int    qParameterList_get_int   (Teuchos::ParameterList* PL, const string& name);
string qParameterList_get_string(Teuchos::ParameterList* PL, const string& name);
double* qParameterList_get_doubles(Teuchos::ParameterList* PL, const string& name);
int* qParameterList_get_ints(Teuchos::ParameterList* PL, const string& name);
