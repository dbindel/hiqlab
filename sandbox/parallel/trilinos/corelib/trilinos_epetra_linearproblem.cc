/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "trilinos_epetra_linearproblem.h"

Epetra_LinearProblem_Complex::Epetra_LinearProblem_Complex(void)
  : Epetra_LinearProblem(), Operatori_(0), Ai_(0), Xi_(0), Bi_(0)
{
}

Epetra_LinearProblem_Complex::Epetra_LinearProblem_Complex
                                (Epetra_RowMatrix* A,   Epetra_RowMatrix* Ai,
                                 Epetra_MultiVector* X, Epetra_MultiVector* Xi,
                                 Epetra_MultiVector* B, Epetra_MultiVector* Bi)
  : Epetra_LinearProblem(A, X, B), Operatori_(0), Ai_(Ai), Xi_(Xi), Bi_(Bi)
{
  Operatori_ = dynamic_cast<Epetra_Operator *>(Ai_); // Try to make matrix an operator
}

Epetra_LinearProblem_Complex::Epetra_LinearProblem_Complex
                                (Epetra_Operator* A,    Epetra_Operator* Ai,
                                 Epetra_MultiVector* X, Epetra_MultiVector* Xi,
                                 Epetra_MultiVector* B, Epetra_MultiVector* Bi)
  : Epetra_LinearProblem(A, X, B), Operatori_(Ai), Ai_(0), Xi_(Xi), Bi_(Bi)
{
  Ai_ = dynamic_cast<Epetra_RowMatrix *>(Operatori_); // Try to make operator a matrix
}

Epetra_LinearProblem_Complex::Epetra_LinearProblem_Complex(Epetra_CrsMatrix_Complex* A,
                                 Epetra_MultiVector_Complex* X,
                                 Epetra_MultiVector_Complex* B)
  : Epetra_LinearProblem(A, X, B), Operatori_(0), Ai_(A->get_Az()),
                                   Xi_(X->get_Vz()), Bi_(B->get_Vz())
{
  Operatori_ = dynamic_cast<Epetra_Operator *>(Ai_); // Try to make matrix an operator
}

Epetra_LinearProblem_Complex::Epetra_LinearProblem_Complex(Epetra_CrsMatrix_Complex* A,
                                 Epetra_Vector_Complex* X,
                                 Epetra_Vector_Complex* B)
  : Epetra_LinearProblem(A, X, B), Operatori_(0), Ai_(A->get_Az()),
                                   Xi_(X->get_Vz()), Bi_(B->get_Vz())
{
  Operatori_ = dynamic_cast<Epetra_Operator *>(Ai_); // Try to make matrix an operator
}

Epetra_LinearProblem_Complex::~Epetra_LinearProblem_Complex(void)
{
}

int Epetra_LinearProblem_Complex::CheckInput() const {

  Epetra_LinearProblem::CheckInput();

  int ierr = 0;
  if (Operatori_==0) ierr = -1;
  if (Xi_==0) ierr = -2;
  if (Bi_==0) ierr = -3;

  EPETRA_CHK_ERR(ierr);  // Return now if any essential objects missing

  if (Ai_==0) EPETRA_CHK_ERR(1); // Return warning error because this problem has no matrix (just an operator)

  if (!Ai_->OperatorDomainMap().SameAs(Xi_->Map())) ierr = -4;
  if (!Ai_->OperatorRangeMap().SameAs(Bi_->Map())) ierr = -5;

  EPETRA_CHK_ERR(ierr);
  return(0);

}
