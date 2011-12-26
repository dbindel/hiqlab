/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef MESHSTUBS_H
#define MESHSTUBS_H

#include <mex.h>
#include "mesh.h"
#include "pzlinear.h"

mxArray* Mesh_assemble_struct1(Mesh* m, int reduced=1);
mxArray* Mesh_assemble_dR1(Mesh* m, const mxArray* Kmat, 
                           double cx, double cv, double ca, 
                           int reduced=1);
mxArray* Mesh_assemble_dR1(Mesh* m, double cx, double cv, double ca,
                           int reduced=1);
mxArray* Mesh_element_dR1(Mesh* m, int eltid,
                          double cx, double cv, double ca);

void Mesh_get_e1 (Mesh* m, int*    e);
void Mesh_get_bc1(Mesh* m, double* bc);
void Mesh_get_u1 (Mesh* m, double* u, double* ui);
void Mesh_get_v1 (Mesh* m, double* u);
void Mesh_get_a1 (Mesh* m, double* u);
void Mesh_get_f1 (Mesh* m, double* f, double* fi);
void Mesh_get_nen_elt(Mesh* m, int* result);

#endif /* MESHSTUBS_H */
