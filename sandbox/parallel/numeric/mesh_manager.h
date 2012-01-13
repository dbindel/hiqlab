/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef MESH_MANAGER_H
#define MESH_MANAGER_H

#include <vector>

#include "mesh_partition.h"
#include "cscmatrix.h"
#include "mesh_add_block.h"
#include "mesh_partitioner.h"


/*@t ----------
 * \section{Mesh_Manager object}
 *
 * A class for testing Mesh objects partitioned into Mesh_Partition
 * objects. It is used mainly to assemble the matrix and rhs.
 *
 *@c*/

class Mesh_Manager {

  public:
    Mesh_Manager();
    virtual ~Mesh_Manager();

    /** Add meshes to the manager
     */
    int add_mesh(Mesh_Partition* mesh);

    /** Compute total number of ids of all meshes added
     */
    int compute_numid();

    /** Set partition mesh with Mesh_Partitioner
     */
    void set_partition(Mesh_Partitioner* meshp);

    /** Initialize
     */
    void initialize();

    /** Set ids
     */
    void set_ids(Mesh_Partitioner* meshp);

    /** Assemble the residual on each partition
     */
    void assemble_R();

    Mesh_Partition* get_mesh(int i) {return meshes[i];}
    int nummesh() {return meshes.size();}

  private:

    std::vector<Mesh_Partition*> meshes;

};


/*@t ----------
 * \section{Building CSCMatrices with Mesh_Manager}
 *
 *@c*/

/** Computing the sparsity structure of the matrix
 */
void build_csc_matrix(Mesh_Manager* mm, int n, 
                      std::vector<int>& jc, std::vector<int>& ir,
                      int reduced);
CSCMatrix* build_csc_matrix(Mesh_Manager* mm, int reduced=1);

/** Assembling the stiffness and mass matrices
 */
CSCMatrix* assemble_dR(Mesh_Manager* mm, int cx, int cv, int ca, int reduced=1);

/** 
 */
void       get_f(Mesh_Manager* mm, double* vr, double* vi=NULL);


/*@t ----------
 * \section{Building CSCMatrix Prolongators with QAddBlockProlongator}
 *
 *@c*/

/** Build the non-zero structure as a CSCMatrix
 */
CSCMatrix* build_csc_matrix(QAddBlockProlongator* prolongator, int m, int n, int reduced=1);

/** Build the prolongator as a CSCMatrix
 */
CSCMatrix* assemble_P(Mesh* from, Mesh* to, int* nemap, int reduced=1);


/*@t ----------
 * \section{Prolongator_Mesh_Manager object}
 * 
 *@c*/

class Prolongator_Mesh_Manager {

  public:
    Prolongator_Mesh_Manager(Mesh_Manager* from, Mesh_Manager* to);
    ~Prolongator_Mesh_Manager();

    void set_nemap(int numnp, int* nmp);
    void construct_partition_nemap(int n, int numnpn, int* nemap_global);

    Mesh_Manager* get_from(){return from;}
    Mesh_Manager* get_to()  {return   to;}

  private:
    Mesh_Manager* from;
    Mesh_Manager* to;
    std::vector<int> NEMAP;
};


/*@t ----------
 * \section{Building CSCMatrix Prolongators with Prolongator_Mesh_Manager}
 *
 *@c*/

/** Build the non-zero structure as a CSCMatrix
 */
CSCMatrix* build_csc_matrix(Prolongator_Mesh_Manager* pmm, int reduced);

/** Build the prolongator as a CSCMatrix
 */
CSCMatrix* assemble_P(Prolongator_Mesh_Manager* pmm, int reduced);

#endif /* MESH_MANAGER_H */
