#ifndef _TRILINOS_ANASAZI_H
#define _TRILINOS_ANASAZI_H

/** @file trilinos_anasazi.h
 *
 * Functions associated with ANASAZI
 */

#include <string.h>

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "AnasaziBasicEigenproblem.hpp"
#include "Teuchos_ParameterList.hpp"
#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_vector.h"

class Mesh;
/** Anasazi eigensolver
 *
 * @param Kshift   Shifted stiffness matrix
 * @param M        Mass matrix
 * @param nev      Number of eigenvalues desired
 * @param d        Array of eigenvalues (nev)
 * @param v        Array of eigenvectors (n-by-nev)
 * @param pl       A parameter list
 */
template<class ST, class MV, class OP>
void create_eigenproblem(
    Teuchos::RefCountPtr< Anasazi::BasicEigenproblem<ST,MV,OP> >& eigenproblem,
    Teuchos::RefCountPtr<OP>& invKshiftM_r,
    Teuchos::RefCountPtr<MV>& ivec, int nev, Teuchos::ParameterList* pl);

template<class ST, class MV, class OP>
int compute_eigs_anasazi_general(OP* invKshiftM,
                                 MV* ivec,
                 int nev, Teuchos::ParameterList* pl,
                 std::vector<Anasazi::Value<ST> >& d, MV* v);


/** Implementations of Anasazi Eigensolver for Real Problem*/
int compute_eigs_anasazi(Epetra_Operator* Op,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector* v);

int compute_eigs_anasazi(Epetra_CrsMatrix* Kshift,
                         Epetra_CrsMatrix* M,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector* v);

int compute_eigs_anasazi(Epetra_CrsMatrix* Kshift,
                         Epetra_CrsMatrix* M, double w0, int form,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector* v);

int compute_eigs_anasazi(Mesh* mesh, double w0,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector* v);

/** Implementations of Anasazi Eigensolver for Complex Problem*/
int compute_eigs_anasazi(Epetra_Operator_Complex* Op,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector_Complex* v);

int compute_eigs_anasazi(Epetra_CrsMatrix_Complex* Kshift,
                         Epetra_CrsMatrix_Complex* M,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector_Complex* v);

int compute_eigs_anasazi(Epetra_CrsMatrix_Complex* Kshift,
                         Epetra_CrsMatrix_Complex* M, double w0, int form,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector_Complex* v);

int compute_eigs_anasazi(Mesh* mesh, double w0,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector_Complex* v);

void undo_spectral_trans(int nev, double sr, double si, int form, double* dr, double*di);

#endif /* _TRILINOS_ANASAZI_H */
