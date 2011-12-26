/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef TEDLINEAR_H
#define TEDLINEAR_H

#include "qcomplex.h"
#include "mesh.h"
#include "coordmatrix.h"

/** Relabel the ID array so all mechanical degrees of freedom come first.
 */
int ted_block_mesh(Mesh* mesh);


/** Form a combination of the A and B matrices in a linearized problem.
 *
 * if T0 != 0,
 *     B = diag(Muu, Muu, Ctt)
 *     A = [0,   -Muu,  0;
 *          Kuu,  0,    Kut;
 *          0,    Ctu,  Ktt];
 * else
 *     B = diag(-Kuu, Muu, -Ctt/T0)
 *     A = [0,    Kuu,     0;
 *          Kuu,  0,       Kut;
 *          0,   -Ctu/T0, -Ktt/T0];
 */
void ted_assemble(CoordMatrix* AB, Mesh* mesh, double T0, double cT,
                  dcomplex ca, dcomplex cb);
CSCMatrix* ted_assemble(Mesh* mesh, double T0, double cT,
                        dcomplex ca, dcomplex cb);


/** Compute damped frequencies for a linearized problem.
 */
int ted_compute_eigs(Mesh* mesh, double w0, double T0, double cT,
                     int nev, int ncv,
                     dcomplex* d, dcomplex* v = NULL);
int ted_compute_eigs(Mesh* mesh, double w0, double T0, double cT,
                     int nev, int ncv,
                     double* dr, double* di,
                     double* vr = NULL, double* vi = NULL);

/** Compute damped frequencies for a linearized problem (pert method)
 */
int ted_compute_eigsp(Mesh* mesh, double w0, int nev, int ncv,
                      dcomplex* d, dcomplex* v = NULL);
int ted_compute_eigsp(Mesh* mesh, double w0, int nev, int ncv,
                      double* dr, double* di,
                      double* vr = NULL, double* vi = NULL);

/** Compute frequencies for the purely mechanical problem.
 */
int ted_compute_eigs_mech(Mesh* mesh, double w0, int nev, int ncv,
                          dcomplex* d, dcomplex* v = NULL);
int ted_compute_eigs_mech(Mesh* mesh, double w0, int nev, int ncv,
                          double* dr, double* di,
                          double* vr = NULL, double* vi = NULL);

#endif /* TEDLINEAR_H */
