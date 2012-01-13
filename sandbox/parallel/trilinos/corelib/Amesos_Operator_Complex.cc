#include <iostream>
#include "Amesos_Operator_Complex.h"
#include "Epetra_Vector.h"

Amesos_Operator_Complex::Amesos_Operator_Complex(
                Epetra_CrsMatrix_Complex* A,
                Epetra_CrsMatrix_Complex* B, int s_type)
  : Label_(0), is_factorized(0), M(B)
{
    string solver[7] = { "Lapack", "Klu", "Umfpack", "Superlu", "Superludist", "Scalapack", "Mumps" };

    Label_ = "Amesos Operator";
    std::cout << "Amesos Operator:" << solver[s_type] << "\n";

    ELP = new Epetra_LinearProblem_Complex();
    ELP->SetOperator(A);
    ELP->SetOperator_i(A->get_Az());

    Amesos_Complex Amesos_Factory;
    ABS = Amesos_Factory.Create(solver[s_type], *ELP);
    assert(ABS);

    factor();
}

Amesos_Operator_Complex::~Amesos_Operator_Complex()
{
    delete ABS;
    delete ELP;
}

void Amesos_Operator_Complex::factor()
{
    ABS->SymbolicFactorization();
    ABS->NumericFactorization();
}

int Amesos_Operator_Complex::Apply(
                 const Epetra_MultiVector_Complex& X,
                       Epetra_MultiVector_Complex& Y) const {
    if (!X.Map().SameAs(OperatorDomainMap())) EPETRA_CHK_ERR(-1);
    if (!Y.Map().SameAs(OperatorRangeMap())) EPETRA_CHK_ERR(-2);
    if (Y.NumVectors()!=X.NumVectors()) EPETRA_CHK_ERR(-3);

    // -- Construct temporary vector
    Epetra_MultiVector_Complex* xtmp;

    // -- Initialize Y
    Y.PutScalar(0.0, 0.0);

    // -- Multiply by M if exists and set RHS
    if (M) {
       xtmp = new Epetra_MultiVector_Complex(OperatorDomainMap(),Y.NumVectors());
       M->Multiply(false, X, *xtmp);
    } else {
       xtmp = new Epetra_MultiVector_Complex(X);
    }

    // -- Set RHS and LHS
    ELP->SetRHS  (xtmp);
    ELP->SetRHS_i(xtmp->get_Vz());
    ELP->SetLHS  (&Y);
    ELP->SetLHS_i(Y.get_Vz());

    // -- Finally solve
    int ierr = ABS->Solve();

    // -- Delete temporaries
    delete xtmp;

    return(ierr);
}
