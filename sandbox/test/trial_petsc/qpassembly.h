/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 * $Id: qpassembly.h,v 1.4 2006/05/15 17:03:42 dbindel Exp $
 */
 
#ifndef QPASSEMBLY_H
#define QPASSEMBLY_H

#include "qassembly.h"
#include "qcomplex.h"
#include "qmatrix.h"
#include <vector>

#include "petscksp.h"

class Mesh;

//const int QP_ASSEMBLE_REMOVE_DIRICHLET   = 1;
//const int QP_ASSEMBLE_IDENTITY_DIRICHLET = 2;

class QPETScAssembler : public QAssembler {
 public:
    QPETScAssembler(Mat matrix, Mesh* mesh, 
                    int reduced);
    ~QPETScAssembler();
    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double*   Ke);
//    void add_identity_dirichlet();
 private:
    std::vector<int> myid;
    std::vector<int> myid_l;
    Mat matrix;
    Mesh* mesh;
    int reduced;
    void map_idvec(int* eltid, int n);
    int get_Vlocal(QMatrix<double>& Vl, int* eltid, int n);
//    void apply_identity_dirichlet(int* eltid, int n, double* Ke);
};

Mat Mesh_assemble_dR_petsc(Mesh* mesh, double cx, double cv, double ca,
                           int mat_type_code, int is_reduced);
Vec Mesh_assemble_R_petsc(Mesh* mesh);

#endif /* QPASSEMBLY_H */
