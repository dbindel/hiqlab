/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef FIELD_EVAL_H
#define FIELD_EVAL_H

class Mesh;

/*@T ---------------------------------------
 * \section{Evaluating fields}
 *
 * There are a variety of functions that we might want to evaluate
 * inside elements: stress, heat flux, power flux, potential gradient,
 * etc.  We might also want some interesting combination of these
 * variables.  The [[FieldEval]] interface provides a uniform way
 * of calling these functions, so we can (for example) write a generic
 * loop to evaluate a function $f(X)$ at all the Gauss points in the
 * mesh and then compute an $L^2$ projection.
 *
 * At present, the only implementation of this interface is
 * [[FieldEvalPower]], which makes a call back to the element
 * [[mean_power]] method.
 *@c*/

class FieldEval {
 public:
    virtual ~FieldEval();
    virtual void operator()(Mesh* mesh, int eltid, double* X, double* fX) = 0;
};

class FieldEvalPower : public FieldEval {
 public:
    ~FieldEvalPower();
    void operator()(Mesh* mesh, int eltid, double* X, double* fX);
};

/*@q*/
#endif /* FIELD_EVAL_H */
