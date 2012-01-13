#include <iostream>
#include "AztecOO_Operator_Komplex.h"
#include "Epetra_Vector.h"
#include "trilinos_komplex.h"
#include "Epetra_Comm.h"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_VectorOut.h"

#include "ml_include.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelOperator.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

AztecOO_Operator_Komplex::AztecOO_Operator_Komplex(
                AztecOO* AZS, Epetra_Map* dm, Epetra_Map* rm,
                Epetra_CrsMatrix_Complex* B)
  : Label_("AztecOO Operator_Komplex"), AZS(AZS), M(B), maxiter(1000), tol(1e-9)
{
    dmap   = new Epetra_Map(*dm);
    rmap   = new Epetra_Map(*rm);
}

AztecOO_Operator_Komplex::AztecOO_Operator_Komplex(
                AztecOO* AZS, Epetra_Map* dm,
                Epetra_CrsMatrix_Complex* B)
  : Label_("AztecOO Operator_Komplex"), AZS(AZS), M(B), maxiter(1000), tol(1e-9)
{
    dmap   = new Epetra_Map(*dm);
    rmap   = new Epetra_Map(*dm);
}

AztecOO_Operator_Komplex::~AztecOO_Operator_Komplex()
{
    delete dmap;
    delete rmap;
}

int AztecOO_Operator_Komplex::Apply(
                 const Epetra_MultiVector_Complex& X,
                       Epetra_MultiVector_Complex& Y) const {
//    if (!X.Map().SameAs(OperatorDomainMap())) EPETRA_CHK_ERR(-1);
//    if (!Y.Map().SameAs(OperatorRangeMap())) EPETRA_CHK_ERR(-2);
    if (Y.NumVectors()!=X.NumVectors()) EPETRA_CHK_ERR(-3);

    // -- Construct temporary vector
    Epetra_MultiVector_Complex* xtmp;
    Epetra_Super_MultiVector* rhs;
    Epetra_MultiVector* lhs;

    // -- Initialize Y
    Y.PutScalar(0.0, 0.0);

    // -- Multiply by M if exists and set RHS
    if (M) {
       xtmp = new Epetra_MultiVector_Complex(X.Map(),X.NumVectors());
       M->Multiply(false, X, *xtmp);
    } else {
       xtmp = new Epetra_MultiVector_Complex(X);
    }

    // -- Construct super vector
    rhs = Complex2SuperMultiVector(xtmp);
    lhs = new Epetra_MultiVector(AZS->GetUserOperator()->OperatorRangeMap(),X.NumVectors());

    // -- Set RHS and LHS
    AZS->SetRHS  (rhs->GetMultiVector());
    AZS->SetLHS  (lhs);

    // -- Finally solve
    int ierr = AZS->recursiveIterate(maxiter,tol);

    // -- Copy into Y
    Multi2ComplexMultiVector(lhs, &Y);

    // -- Show results
    if( X.Comm().MyPID() == 0 ) {
        cout << "Solver performed " << AZS->NumIters()
             << "iterations.\n";
        cout << "Norm of the true residual = " << AZS->TrueResidual() << endl;
    }

    // -- Delete temporaries
    delete lhs;
    delete rhs;
    delete xtmp;

    return(ierr);
}
