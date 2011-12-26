/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef STUBS_H
#define STUBS_H

#include <mex.h>
#include "cscmatrix.h"

mxArray* CSCMatrix_to_mx(CSCMatrix* K);
CSCMatrix* mx_to_CSCMatrix(const mxArray* Kmat);

#endif /* STUBS_H */
