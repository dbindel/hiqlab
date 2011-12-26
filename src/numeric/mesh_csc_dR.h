/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef MESH_CSC_DR_H
#define MESH_CSC_DR_H

#include "mesh.h"
#include "coordmatrix.h"
#include "cscmatrix.h"

CSCMatrix* assemble_dR(Mesh* mesh, 
                       double cx=1, double cv=0, double ca=0, int reduced=1);

CSCMatrix* element_dR(Mesh* mesh, int eltid,
                      double cx=1, double cv=0, double ca=0, int reduced=1);

#endif /* MESH_CSC_DR_H */
