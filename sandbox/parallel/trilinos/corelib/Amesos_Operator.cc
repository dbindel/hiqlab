#include <iostream>
#include "Amesos_Operator.h"
#include "Epetra_Vector.h"

Amesos_Operator::Amesos_Operator(Epetra_RowMatrix* A, Epetra_RowMatrix* B, int s_type)
  : Label_(0), is_factorized(0), M(B)
{
    string solver[7] = { "Lapack", "Klu", "Umfpack", "Superlu", "Superludist", "Scalapack", "Mumps" };

    Label_ = "Amesos Operator";
    std::cout << "Amesos Operator:" << solver[s_type] << "\n";

    ELP = new Epetra_LinearProblem();
    ELP->SetOperator(A);

    Amesos Amesos_Factory;
    ABS = Amesos_Factory.Create(solver[s_type], *ELP);
    assert(ABS);

    factor();
}

Amesos_Operator::~Amesos_Operator()
{
    delete ABS;
    delete ELP;
}

void Amesos_Operator::factor()
{
    ABS->SymbolicFactorization();
    ABS->NumericFactorization();
}

int Amesos_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

    if (!X.Map().SameAs(OperatorDomainMap())) EPETRA_CHK_ERR(-1);
    if (!Y.Map().SameAs(OperatorRangeMap())) EPETRA_CHK_ERR(-2);
    if (Y.NumVectors()!=X.NumVectors()) EPETRA_CHK_ERR(-3);

    // -- Construct temporary vector
    Epetra_MultiVector xtmp(OperatorDomainMap(), Y.NumVectors());

    // -- Initialize Y
    Y.PutScalar(0.0);

    // -- Multiply by M if exists
    if (M)
       M->Multiply(false, X, xtmp);

    // -- Set RHS and LHS
    ELP->SetLHS(&Y);
    ELP->SetRHS(&xtmp);

    // -- Finally solve
    int ierr = ABS->Solve();

    return(ierr);
}
