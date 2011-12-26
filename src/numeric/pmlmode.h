/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLMODE_H
#define PMLMODE_H

/** @file pmlmode.h
 *
 * Functions associated with ARPACK calls.
 */

#include <string.h>

#include "time.h"
#include "qcomplex.h"
#include "cscmatrix.h"
#include "mesh.h"

/** Shift and invert Arnoldi eigensolver
 *
 * @param Kshift   Shifted stiffness matrix
 * @param M        Mass matrix
 * @param nev      Number of eigenvalues desired
 * @param ncv      Maximum dimension of Arnoldi decomposition (before restart)
 * @param d        Array of eigenvalues (nev)
 * @param v        Array of eigenvectors (n-by-nev)
 * @param ldv      Leading dimension of v
 */
int compute_eigs(CSCMatrix& Kshift, CSCMatrix& M, int nev, int ncv,
                 dcomplex* d, dcomplex* v, int ldv);

/** Shift and invert Arnoldi eigensolver
 *
 * @param mesh   Mesh to use in assembling mass and stiffness
 * @param w0     Initial shift (in Hz)
 * @param nev    Number of frequencies desired
 * @param ncv    Maximum dimension of Arnoldi decomposition (before restart)
 * @param dr     Real part of computed frequencies in Hz (nev)
 * @param di     Imag part of computed frequencies in Hz (nev)
 * @param vr     Real part of eigenvector                (n*nev)
 * @param vi     Imag part of eigenvector                (n*nev)
 */
int compute_eigs(Mesh* mesh, double w0, int nev, int ncv,
                 double* dr, double* di,
                 double* vr = NULL, double* vi = NULL);

#endif /* PMLMODE_H */
