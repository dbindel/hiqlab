#ifndef EPETRA_OPERATOR_COMPLEX_H
#define EPETRA_OPERATOR_COMPLEX_H

class Epetra_MultiVector_Complex;
class Epetra_Map;
class Epetra_Comm;

/** Epetra_Operator_Complex: A pure virtual class for using complex-valued
 *  double-precision operators.
 *  The Epetra_Operator_Complex class is a pure virtual class
 *  (specifies interface only) that enable the use of complex-valued
 *  double-precision operators.
 */

class Epetra_Operator_Complex {

 public:

    virtual ~Epetra_Operator_Complex() {};

    /** This flag allows the transpose of the given operator to be used
     *  implicitly.  Setting this flag affects only the Apply() and
     *  ApplyInverse() methods.  If the implementation of this interface
     *  does not support transpose use, this method should return a
     *  value of -1.
     */
    virtual int SetUseTranspose(bool UseTranspose) = 0;

    /** Returns the result of a Epetra_Operator_Complex applied to a
     *  Epetra_MultiVector_Complex X in Y.
     *  Return Integer error code, set to 0 if successful.
     */
    virtual int Apply(const Epetra_MultiVector_Complex& X,
                            Epetra_MultiVector_Complex& Y) const = 0;

    /** Returns the result of a Epetra_Operator_Complex inverse applied
     *  to an Epetra_MultiVector_Complex X in Y.
     *  Return Integer error code, set to 0 if successful.
     */
    virtual int ApplyInverse(const Epetra_MultiVector_Complex& X,
                                   Epetra_MultiVector_Complex& Y) const = 0;

    /** Returns the infinity norm of the global matrix.
     *  Warning This method must not be called unless HasNormInf() returns true.
     */
    virtual double NormInf() const = 0;

    /** Returns a character string describing the operator */
    virtual const char * Label() const = 0;

    /** Returns the current UseTranspose setting. */
    virtual bool UseTranspose() const = 0;

    /** Returns true if the this object can provide an approximate
     *  Inf-norm, false otherwise.
     */
    virtual bool HasNormInf() const = 0;

    /** Returns a pointer to the Epetra_Comm communicator associated
     *  with this operator.
     */
    virtual const Epetra_Comm & Comm() const = 0;

    /** Returns the Epetra_Map object associated with the domain of
     *  this operator.
     */
    virtual const Epetra_Map & OperatorDomainMap() const = 0;

    /** Returns the Epetra_Map object associated with the range of
     *  this operator.
     */
    virtual const Epetra_Map & OperatorRangeMap() const = 0;

};

#endif /* EPETRA_OPERATOR_COMPLEX_H */
