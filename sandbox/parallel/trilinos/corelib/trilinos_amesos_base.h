#ifndef _TRILINOS_AMESOS_H
#define _TRILINOS_AMESOS_H
#include "trilinos_epetra_linearproblem.h"

#include "Teuchos_ParameterList.hpp"
class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Comm;

/** Amesos_BaseSolver_Complex: A pure virtual class for direct
 *  solution of complex-valued double-precision operators.
 *  See Amesos_BaseSolver.h for details
 */

class Amesos_BaseSolver_Complex {

 public:

  virtual ~Amesos_BaseSolver_Complex() {};

  /** Performs SymbolicFactorization on the matrix A.*/
  virtual int SymbolicFactorization() = 0;

  /** Performs NumericFactorization on the matrix A. */
  virtual int NumericFactorization() = 0;

  /** Solve the problem */
  virtual int Solve() = 0;

  /** If set true, X will be set to the solution of
   *  Return a value of -1 for no support
   *  Return a value of  0 if successful.
   */
  virtual int SetUseTranspose(bool UseTranspose) = 0;

  /** Returns the current UseTranspose setting. */
  virtual bool UseTranspose() const = 0;

  /** Set parameters */
  virtual int SetParameters( Teuchos::ParameterList &ParameterList ) = 0 ;

  /** Returns the Epetra_LinearProblem. */
  virtual const Epetra_LinearProblem_Complex* GetProblem() const = 0;

  /** Returns true if the solver can handle this matrix shape.  */
  virtual bool MatrixShapeOK() const = 0;

  /** Returns a pointer to the Epetra_Comm communicator with this operator */
  virtual const Epetra_Comm & Comm() const = 0;

};

/** Amesos_Complex:  A method for binding a third party direct solver
 *                   to an Epetra_LinearProblem_Complex.
 */

class Amesos_Complex {

 public:
    /** Creates an instance of the Amesos_BaseSolver class
     *  specified by ClassType.
     */
    Amesos_BaseSolver_Complex* Create(const char *ClassType,
                      const Epetra_LinearProblem_Complex& LinearProblem );

    Amesos_BaseSolver_Complex* Create(const string CT,
                      const Epetra_LinearProblem_Complex& LinearProblem );

    /** Queries whether a given interface is avaiable or not.
     */
    bool Query(const char * ClassType);
    bool Query(const string CT);

};

#endif /* _TRILINOS_AMESOS_H */
