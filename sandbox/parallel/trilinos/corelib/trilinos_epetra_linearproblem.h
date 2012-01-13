#ifndef _TRILINOS_EPETRA_H
#define _TRILINOS_EPETRA_H

/** @file trilinos_epetra_linearproblem.h
 *
 *  Classes and functions related to Epetra
 */

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Operator.h"

#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_vector.h"

/** Epetra_LinearProblem_Complex:  The Epetra Linear Problem Class for
 *                                 solving complex-valued problems
 */

class Epetra_LinearProblem_Complex : public Epetra_LinearProblem {
 public:

    /**  Epetra_LinearProblem Default Constructor.
     */
    Epetra_LinearProblem_Complex(void);

    Epetra_LinearProblem_Complex(Epetra_RowMatrix* A,   Epetra_RowMatrix* Ai,
                                 Epetra_MultiVector* X, Epetra_MultiVector* Xi,
                                 Epetra_MultiVector* B, Epetra_MultiVector* Bi);

    Epetra_LinearProblem_Complex(Epetra_Operator* A,    Epetra_Operator* Ai,
                                 Epetra_MultiVector* X, Epetra_MultiVector* Xi,
                                 Epetra_MultiVector* B, Epetra_MultiVector* Bi);

    Epetra_LinearProblem_Complex(Epetra_CrsMatrix_Complex* A,
                                 Epetra_MultiVector_Complex* X,
                                 Epetra_MultiVector_Complex* B);

    Epetra_LinearProblem_Complex(Epetra_CrsMatrix_Complex* A,
                                 Epetra_Vector_Complex* X,
                                 Epetra_Vector_Complex* B);

    virtual ~Epetra_LinearProblem_Complex(void);

    int CheckInput() const;

    /** Setters for imaginary parts
     */
    void SetOperator_i(Epetra_RowMatrix* Ai) {Ai_ = Ai; Operatori_ = Ai;};
    void SetOperator_i(Epetra_Operator* Ai) {Ai_ = dynamic_cast<Epetra_RowMatrix *>(Ai); Operatori_ = Ai;};
    void SetLHS_i(Epetra_MultiVector * Xi) {Xi_ = Xi;};
    void SetRHS_i(Epetra_MultiVector * Bi) {Bi_ = Bi;};

    /** Getters for imaginary parts
     */
    Epetra_Operator * GetOperator_i() const {return(Operatori_);};
    Epetra_RowMatrix * GetMatrix_i() const {return(Ai_);};
    Epetra_MultiVector * GetLHS_i() const {return(Xi_);};
    Epetra_MultiVector * GetRHS_i() const {return(Bi_);};

 private:

    Epetra_Operator * Operatori_;
    Epetra_RowMatrix * Ai_;
    Epetra_MultiVector * Xi_;
    Epetra_MultiVector * Bi_;

    Epetra_LinearProblem_Complex & operator=(const Epetra_LinearProblem_Complex& Problem);

};

#endif /* _TRILINOS_EPETRA_H */
