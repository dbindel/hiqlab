/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#ifndef TRILINOS_ARPACK_H
#define TRILINOS_ARPACK_H

/** @file trilinos_arpack.h
 *
 * Functions associated with ARPACK calls.
 */

#include <string.h>

#include "time.h"
#include "qcomplex.h"
#include "cscmatrix.h"
#include "mesh.h"

#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"

#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_vector.h"

int compute_eigs_arpack(Epetra_RowMatrix* Kshift, Epetra_RowMatrix* M,
                        int nev, int ncv,  double* dr, double* di,
                        Epetra_MultiVector* vri);

int compute_eigs_arpack(Epetra_RowMatrix* Kshift, Epetra_RowMatrix* M,
                        int nev, int ncv,  double* d,
                        Epetra_MultiVector* vr);

int compute_eigs_arpack(Epetra_CrsMatrix_Complex* Kshift, Epetra_CrsMatrix_Complex* M,
                        int nev, int ncv,  double* dr, double* di,
                        Epetra_MultiVector_Complex* vz);

int compute_eigs_arpack(Epetra_CrsMatrix_Complex* Kshift, Epetra_CrsMatrix_Complex* M,
                        int nev, int ncv, dcomplex* dz,
                        Epetra_MultiVector_Complex* vz);


int compute_eigs_arpack(Mesh* mesh, double w0, int nev, int ncv,
                        double* dr, double* di,
                        Epetra_MultiVector* vri);

int compute_eigs_arpack(Mesh* mesh, double w0, int nev, int ncv,
                        double* dr, double* di,
                        Epetra_MultiVector_Complex* vz);

#endif /* TRILINOS_ARPACK_H */
