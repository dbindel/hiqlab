#ifndef _AZTECOO_OPERATOR_KOMPLEX_H_
#define _AZTECOO_OPERATOR_KOMPLEX_H_

#include <iostream>

#include "Epetra_Map.h"
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"

#include "trilinos_epetra_vector.h"
#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_operator.h"
#include "trilinos_epetra_linearproblem.h"
#include "trilinos_super_matrix.h"
#include "trilinos_super_vector.h"

class Epetra_BlockMap;
class Epetra_Comm;

/** AztecOO_Operator_Komplex: An implementation of the
 *  Epetra_Operator_Complex class.
 */

class AztecOO_Operator_Komplex: public virtual Epetra_Operator_Complex {

 public:

    /** Constructor. */
    AztecOO_Operator_Komplex(
                    AztecOO* AZS, Epetra_Map* dmap,
                    Epetra_CrsMatrix_Complex* B);
    AztecOO_Operator_Komplex(
                    AztecOO* AZS, Epetra_Map* dmap, Epetra_Map* rmap,
                    Epetra_CrsMatrix_Complex* B);
    virtual ~AztecOO_Operator_Komplex();

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
    const Epetra_Comm & Comm() const{return(AZS->GetUserOperator()->Comm());};

    /** Returns the Epetra_BlockMap object associated with the domain of
     * this matrix operator.
     */
    const Epetra_Map & OperatorDomainMap() const
    {
      if (UseTranspose()) return(*rmap);
      else return(*dmap);
    }

    /** Returns the Epetra_BlockMap object associated with the range of
     * this matrix operator.
     */
    const Epetra_Map & OperatorRangeMap() const
    {
      if (UseTranspose()) return(*dmap);
      else return(*rmap);
    }

    void SetMaxIter(int n) {maxiter=n;};
    void SetTol(double t)  {tol=t;};

 protected:

  AztecOO*  AZS;
  Epetra_CrsMatrix_Complex* M;
  Epetra_Map* dmap;
  Epetra_Map* rmap;
  char * Label_;
  int maxiter;
  double tol;

};

#endif /* _AZTECOO_OPERATOR_KOMPLEX_H_ */

