/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */
 
#ifndef QPASSEMBLY_PART_PETSC_H
#define QPASSEMBLY_PART_PETSC_H

#include <vector>

#include "petscksp.h"

#include "qassembly.h"

class Mesh_Partition;
class QAddBlockProlongator_Partition;

extern const int QP_ASSEMBLE_NO_DIRICHLET;
extern const int QP_ASSEMBLE_REMOVE_DIRICHLET;
extern const int QP_ASSEMBLE_IDENTITY_DIRICHLET;

class QPETScPartitionStructAssembler : public QStructAssembler {

  public:
    QPETScPartitionStructAssembler(Mesh_Partition* mesh, PetscInt Istart, PetscInt Iend1, QStructAssembler& assembler, int is_reduced)
                               : mesh(mesh), assembler(assembler), Istart(Istart), Iend1(Iend1), is_reduced(is_reduced) {}
    virtual ~QPETScPartitionStructAssembler();

    void add(int i, int j);
    void add(int* eltid, int n);
    void add(int* eltidm, int m, int* eltidn, int n);

  private:
    Mesh_Partition* mesh;
    QStructAssembler& assembler;
    PetscInt Istart, Iend1;
    int is_reduced;

    void map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn);
    void filter_eltid(int* eltid, int n, int sid, int eid1, int offset);
};


class QPETScPartitionAssembler : public QAssembler {

  public:
    QPETScPartitionAssembler(Mesh_Partition* mesh, Mat A, PetscInt Istart, PetscInt Iend1, int is_reduced)
                         : mesh(mesh), A(A), Istart(Istart), Iend1(Iend1), is_assemble_real(1), is_reduced(is_reduced) {}
    virtual ~QPETScPartitionAssembler();

    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double*   Ke);
    void add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, double*   Ke);

    void assemble_real() {is_assemble_real=1;}
    void assemble_imag() {is_assemble_real=0;}

    template<class T>
    void add_identity_dirichlet(T diag);

  private:
    Mesh_Partition* mesh;
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


class QPETScPartitionProlongatorStructAssembler : public QStructAssembler {

  public:
    QPETScPartitionProlongatorStructAssembler(QAddBlockProlongator_Partition* prolongator, PetscInt Istart, PetscInt Iend1, QStructAssembler& assembler, int is_reduced)
                               : prolongator(prolongator), assembler(assembler), Istart(Istart), Iend1(Iend1), is_reduced(is_reduced) {}
    virtual ~QPETScPartitionProlongatorStructAssembler();

    void add(int i, int j);
    void add(int* eltid, int n);
    void add(int* eltidm, int m, int* eltidn, int n);

  private:
    QAddBlockProlongator_Partition* prolongator;
    QStructAssembler& assembler;
    PetscInt Istart, Iend1;
    int is_reduced;

    void map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn);
    void filter_eltid(int* eltid, int n, int sid, int eid1, int offset);
};


class QPETScPartitionProlongatorAssembler : public QAssembler {

  public:
    QPETScPartitionProlongatorAssembler(QAddBlockProlongator_Partition* prolongator, Mat A, PetscInt Istart, PetscInt Iend1, int is_reduced)
                         : prolongator(prolongator), A(A), Istart(Istart), Iend1(Iend1), is_reduced(is_reduced) {}
    virtual ~QPETScPartitionProlongatorAssembler();

    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double*   Ke);
    void add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, double*   Ke);

  private:
    QAddBlockProlongator_Partition* prolongator;
    Mat A;
    PetscInt Istart, Iend1;
    int is_reduced;

    template<class T>
    void add(int* eltidm, int m, int* eltidn, int n, T* Ke);
    template<class T>
    void map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn, T* Ke);
    void filter_eltid(int* eltid, int n, int sid, int eid1, int offset);
};


Mat Mesh_Partition_assemble_dR_petsc(Mesh_Partition* mesh, double cxr, double cvr, double car,
                                       int mat_type_code, int is_reduced);
Mat Mesh_Partition_assemble_dR_petsc(Mesh_Partition* mesh, double cxr, double cxi, 
                                       double cvr, double cvr,
                                       double car, double cai,
                           int mat_type_code, int is_reduced);
Vec Mesh_Partition_assemble_R_petsc(Mesh_Partition* mesh, int is_reduced);
Mat Mesh_Partition_assemble_P_petsc(Mesh_Partition* from, Mesh_Partition* to, int* nemap, int mat_type_code, int is_reduced);

#endif /* QPASSEMBLY_PART_PETSC_H */
