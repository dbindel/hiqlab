/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef QLAPACK_H
#define QLAPACK_H

#include "qarray.h"
#include "qiarray.h"

void dgemm(char* trans, double alpha, QArray* A, QArray* B,
           double beta, QArray* C);
int dgetrf(QArray* A, QIArray* ipiv);
int dgetrs(char* trans, QArray* A, QIArray* ipiv, QArray* b);
int dgeev(QArray* A, QArray* wr, QArray* vl, QArray* vr);

// dgees
// dpotrf / dpotrs
// dposv
// dsyev

#endif /* QLAPACK_H */
