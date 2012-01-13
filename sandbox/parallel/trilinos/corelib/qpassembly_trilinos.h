/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#ifndef QPASSEMBLY_TRILINOS_H
#define QPASSEMBLY_TRILINOS_H

#include <vector>

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "trilinos_epetra_vector.h"
#include "trilinos_epetra_matrix.h"

class Mesh;

/** Compute number of non-zeros in local process
 */
int* Mesh_trilinos_my_nnz(Mesh* mesh, double cx, double cv, double ca,
                                      int is_reduced, Epetra_Map* Map);

/** Compute matrix of local process
 */
void Mesh_assemble_my_dR_trilinos(Epetra_CrsMatrix* Ar, Epetra_CrsMatrix* Ai, Mesh* mesh,
                                double cx, double cv, double ca, int is_reduced);

/** Compute matrix
 */
Epetra_CrsMatrix_Complex* Mesh_assemble_dRz_trilinos(Mesh* mesh, double cx, double cv, double ca,
                           int is_reduced, int is_real=0);

/** Compute matrix (real part)
 */
Epetra_CrsMatrix* Mesh_assemble_dR_trilinos
    (Mesh* mesh, double cx, double cv, double ca, int is_reduced);

/** Compute matrix (imag part)
 */
Epetra_CrsMatrix* Mesh_assemble_dRi_trilinos
    (Mesh* mesh, double cx, double cv, double ca, int is_reduced);

//FIXME: Whole vector assembled locally and distributed
/** Compute residual vector of local process
 */
void Mesh_assemble_my_R_trilinos(Epetra_Vector* Vr, Epetra_Vector* Vi, Mesh* mesh);

/** Compute residual vector
 */
Epetra_Vector_Complex* Mesh_assemble_Rz_trilinos(Mesh* mesh, int is_real=0);

/** Compute residual vector (real part)
 */
Epetra_Vector* Mesh_assemble_R_trilinos(Mesh* mesh);

/** Compute residual vector (imag part)
 */
Epetra_Vector* Mesh_assemble_Ri_trilinos(Mesh* mesh);

#endif /* QPASSEMBLY_TRILINOS_H */
