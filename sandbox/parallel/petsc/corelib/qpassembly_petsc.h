/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */
 
#ifndef QPASSEMBLY_PETSC_H
#define QPASSEMBLY_PETSC_H

#include <vector>

#include "petscksp.h"

#include "qassembly.h"

class Mesh;
class QAddBlockProlongator;

using std::vector;

const int QP_ASSEMBLE_NO_DIRICHLET       = 0;
const int QP_ASSEMBLE_REMOVE_DIRICHLET   = 1;
const int QP_ASSEMBLE_IDENTITY_DIRICHLET = 2;


class QPETScStructAssembler : public QStructAssembler {

  public:
    QPETScStructAssembler(Mesh* mesh, PetscInt Istart, PetscInt Iend1, QStructAssembler& assembler, int is_reduced)
                               : mesh(mesh), assembler(assembler), Istart(Istart), Iend1(Iend1), is_reduced(is_reduced) {}
    virtual ~QPETScStructAssembler();

    void add(int i, int j);
    void add(int* eltid, int n);
    void add(int* eltidm, int m, int* eltidn, int n);

  private:
    Mesh* mesh;
    QStructAssembler& assembler;
    PetscInt Istart, Iend1;
    int is_reduced;

    void map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn);
    void filter_eltid(int* eltid, int n, int sid, int eid1, int offset);
};


class QPETScAssembler : public QAssembler {

  public:
    QPETScAssembler(Mesh* mesh, Mat A, PetscInt Istart, PetscInt Iend1, int is_reduced)
                         : mesh(mesh), A(A), Istart(Istart), Iend1(Iend1), is_assemble_real(1), is_reduced(is_reduced) {}
    virtual ~QPETScAssembler();

    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double*   Ke);
    void add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, double*   Ke);

    void assemble_real() {is_assemble_real=1;}
    void assemble_imag() {is_assemble_real=0;}

    template<class T>
    void add_identity_dirichlet(T diag);

  private:
    Mesh* mesh;
    Mat A;
    PetscInt Istart, Iend1;
    int is_assemble_real;
    int is_reduced;

    template<class T>
    void add(int* eltidm, int m, int* eltidn, int n, T* Ke);
    template<class T>
    void map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn, T* Ke);
    void filter_eltid(int* eltid, int n, int sid, int eid1, int offset);
    template<class T>
    void apply_identity_dirichlet(int* eltid, int n, T* Ke);

};


class QPETScProlongatorStructAssembler : public QStructAssembler {

  public:
    QPETScProlongatorStructAssembler(QAddBlockProlongator* prolongator, PetscInt Istart, PetscInt Iend1, QStructAssembler& assembler, int is_reduced)
                               : prolongator(prolongator), assembler(assembler), Istart(Istart), Iend1(Iend1), is_reduced(is_reduced) {}
    virtual ~QPETScProlongatorStructAssembler();

    void add(int i, int j);
    void add(int* eltid, int n);
    void add(int* eltidm, int m, int* eltidn, int n);

  private:
    QAddBlockProlongator* prolongator;
    QStructAssembler& assembler;
    PetscInt Istart, Iend1;
    int is_reduced;

    void map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn);
    void filter_eltid(int* eltid, int n, int sid, int eid1, int offset);
};


class QPETScProlongatorAssembler : public QAssembler {

  public:
    QPETScProlongatorAssembler(QAddBlockProlongator* prolongator, Mat A, PetscInt Istart, PetscInt Iend1, int is_reduced)
                         : prolongator(prolongator), A(A), Istart(Istart), Iend1(Iend1), is_reduced(is_reduced) {}
    virtual ~QPETScProlongatorAssembler();

    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double*   Ke);
    void add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, double*   Ke);

  private:
    QAddBlockProlongator* prolongator;
    Mat A;
    PetscInt Istart, Iend1;
    int is_reduced;

    template<class T>
    void add(int* eltidm, int m, int* eltidn, int n, T* Ke);
    template<class T>
    void map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn, T* Ke);
    void filter_eltid(int* eltid, int n, int sid, int eid1, int offset);
};


Mat Mesh_assemble_dR_petsc(Mesh* mesh, double cxr, double cvr, double car,
                                       int mat_type_code, int is_reduced);
Mat Mesh_assemble_dR_petsc(Mesh* mesh, double cxr, double cxi, 
                                       double cvr, double cvr,
                                       double car, double cai,
                           int mat_type_code, int is_reduced);
Vec Mesh_assemble_R_petsc(Mesh* mesh, int is_reduced);
Mat Mesh_assemble_P_petsc(Mesh* from, Mesh* to, int* nemap, int mat_type_code, int is_reduced);
Mat MatRemoveZeros(Mat A, double atol);

#endif /* QPASSEMBLY_PETSC_H */
