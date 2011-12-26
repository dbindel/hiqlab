/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PZLINEAR_H
#define PZLINEAR_H

#include "qcomplex.h"
#include "mesh.h"

/** Relabel the ID array so all mechanical degrees of freedom come first.
 */
int pz_block_mesh(Mesh* mesh);
int em_block_mesh(Mesh* mesh);
int em_circuit_block_mesh(Mesh* mesh);

/** Compute damped frequencies for a linearized problem.
 */
int pz_compute_eigs(Mesh* mesh, double w0, int nev, int ncv,
                    dcomplex* d, dcomplex* v = NULL);
int pz_compute_eigs(Mesh* mesh, double w0, int nev, int ncv,
                    double* dr, double* di,
                    double* vr = NULL, double* vi = NULL);

/** Compute frequencies for the purely mechanical problem.
 */
int pz_compute_eigs_mech(Mesh* mesh, double w0, int nev, int ncv,
                         dcomplex* d, dcomplex* v = NULL);
int pz_compute_eigs_mech(Mesh* mesh, double w0, int nev, int ncv,
                         double* dr, double* di,
                         double* vr = NULL, double* vi = NULL);

#endif /* PZLINEAR_H */
