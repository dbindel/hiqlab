#ifndef _AMESOS_OPERATOR_COMPLEX_H_
#define _AMESOS_OPERATOR_COMPLEX_H_

#include <iostream>

#include "Epetra_Map.h"

#include "trilinos_epetra_vector.h"
#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_operator.h"
#include "trilinos_epetra_linearproblem.h"
#include "trilinos_amesos_base.h"

class Epetra_BlockMap;
class Epetra_Comm;

/** Amesos_Operator_Complex: An implementation of the
 *  Epetra_Operator_Complex class.
 */

class Amesos_Operator_Complex: public virtual Epetra_Operator_Complex {

 public:

    /** Constructor. */
    Amesos_Operator_Complex(
                    Epetra_CrsMatrix_Complex* A,
                    Epetra_CrsMatrix_Complex* B, int s_type=6);
    virtual ~Amesos_Operator_Complex();

    /** Multiply by the transpose
     * (This method has no effect and returns -1 as error code.)
     */
    int SetUseTranspose(bool UseTranspose){return(-1);};

    /** Apply inverse to a Epetra_MultiVector_Complex X
     *  (This method has no effect and returns -1 as error code.)
     */
    int ApplyInverse(const Epetra_MultiVector_Complex& X,
                           Epetra_MultiVector_Complex& Y) const{return 0;};

    /** Apply to a Epetra_MultiVector_Complex X
     */
    int Apply(const Epetra_MultiVector_Complex& X,
                    Epetra_MultiVector_Complex& Y) const;

    /** Returns the infinity norm of the global matrix.
     *warning This method must not be called unless HasNormInf() returns true.
     */
    double NormInf() const {return(0.0);};

    /** Returns a character string describing the operator */
    const char * Label() const {return(Label_);};

    /** Returns the current UseTranspose setting. */
    bool UseTranspose() const {return(false);};

    /** Returns true if this object can provide an approximate Inf-norm, false otherwise. */
    bool HasNormInf() const{return(false);};

    /** Returns a pointer to the Epetra_Comm communicator associated
     * with this operator.
     */
    const Epetra_Comm & Comm() const{return(ELP->GetOperator()->Comm());};

    /** Returns the Epetra_BlockMap object associated with the domain of
     * this matrix operator.
     */
    const Epetra_Map & OperatorDomainMap() const
    {
      if (UseTranspose()) return(ELP->GetOperator()->OperatorRangeMap());
      else return(ELP->GetOperator()->OperatorDomainMap());
    }

    /** Returns the Epetra_BlockMap object associated with the range of
     * this matrix operator.
     */
    const Epetra_Map & OperatorRangeMap() const
    {
      if (UseTranspose()) return(ELP->GetOperator()->OperatorDomainMap());
      else return(ELP->GetOperator()->OperatorRangeMap());
    }

 protected:

  Amesos_BaseSolver_Complex*  ABS;
  Epetra_LinearProblem_Complex* ELP;
  Epetra_CrsMatrix_Complex* M;
  char * Label_;
  int is_factorized;

  void factor();
};

#endif /* _AMESOS_OPERATOR_COMPLEX_H_ */

