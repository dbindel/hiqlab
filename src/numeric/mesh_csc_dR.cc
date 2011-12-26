/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include "mesh_csc_dR.h"
#include "coordmatrix.h"


CSCMatrix* assemble_dR(Mesh* mesh, double cx, double cv, double ca, 
                       int reduced)
{
    CoordMatrix assembler(mesh->get_numid());
    mesh->assemble_dR(&assembler, cx, cv, ca, reduced);
    return assembler.to_sparse();
}


CSCMatrix* element_dR(Mesh* mesh, int eltid,
                      double cx, double cv, double ca, int reduced)
{
    CoordMatrix mat(mesh->get_numid());
    if (reduced) {
        QReduceAssembler rassembler(mesh, mat);
        mesh->etype(eltid)->assemble_dR(mesh, eltid, &rassembler, cx, cv, ca);
    } else {
        mesh->etype(eltid)->assemble_dR(mesh, eltid, &mat, cx, cv, ca);
    }
    return mat.to_sparse();
}

