#ifndef _QSLEPC_H
#define _QSLEPC_H

/** @file qslepc.h
 *
 * Functions associated with SLEPC
 */

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "slepceps.h"

class Mesh;

/** Slepc eigensolver
 *
 * @param Kshift   Shifted stiffness matrix
 * @param M        Mass matrix
 * @param nev      Number of eigenvalues desired
 * @param d        Array of eigenvalues (nev)
 * @param v        Array of eigenvectors (n-by-nev)
 */
int compute_eigs_slepc(Mat Kshift, Mat M, int nev,
                double* dr, double* di, Vec vz);

int compute_eigs_slepc(Mesh* mesh, Mat Kshift, Mat M, int nev,
                double* dr, double* di, Vec vz);

// -- Compute with mesh directly

int compute_eigs_slepc(Mesh* mesh, double w0,
                 int nev, double* dr, double* di,
                 Vec vz);

#endif /* _QSLEPC_H */
