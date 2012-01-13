#ifndef MESH_PARTITIONERP_H
#define MESH_PARTITIONERP_H

#include <vector>

#include "mesh.h"
#include "mesh_partitioner.h"
#include "mesh_partition.h"

class Parallel_Mesh_Partitioner {

  public:
    Parallel_Mesh_Partitioner(Mesh* mesh, int nparts);
    virtual ~Parallel_Mesh_Partitioner();

    /** Tie mesh nodes in a range.  Only affects the element connectivity
     *  array; tied nodes are not removed from the mesh data structure.
     */
    void tie(double tol, int start=-1, int end=-1);

    /** Construct disjoint nodal partition of all nodes using METIS
     *  ptype = 0(Recursive), 1(Kway)
     */
    void partition_nodes(int ptype);

    /** Construct overlapping element partition. Elements can be duplicated
     *  on seperate partitions but are owned by only one parition. This is 
     *  to incorporate branch variables
     */
    void partition_elements();

    /** Send an receive information concerning nodes. Uses Non-blocking send
     *  and receives but waits at the end of each function.
     */
    void send_nodes_info(int sproc, int numpart);
    void receive_nodes_info(int rproc, int numpart);
    void sendreceive_nodes_info();
    void sendreceive_nodes_info_root();

    /** Send an receive information concerning elements. Uses Non-blocking send
     *  and receives but waits at the end of each function.
     */
    void send_elements_info(int sproc, int numpart);
    void receive_elements_info(int rproc, int numpart);
    void sendreceive_elements_info();
    void sendreceive_elements_info_root();

    /** Get each partition
     */
    int num_mymeshes() {return meshes.size();}
    Mesh_Partition* get_mesh(int i){return meshes[i];}

    /** Assigning IDG
     */
    void assign_ids();
    void assign_onp_idg();
    void assign_offp_idg();
    void allgather_idg_info();

  private:

    typedef std::vector<int> ivec;

    Mesh* mesh;               // mesh to partition
    Mesh_Partitioner* meshp;  // serial mesh partitioner

    int  process_id;          // process id
    int  numprocess;          // number of processes
    ivec PP;                  // processor partitioning
    int  nparts;              // number of total partitions

    std::vector<Mesh_Partition*> meshes; // partitions

    ivec IDP;                 // id partitioning
};

#endif /* MESH_PARTITIONERP_H */
