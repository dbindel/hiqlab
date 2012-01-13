#ifndef MESH_PARTITIONER_METIS_H 
#define MESH_PARTITIONER_METIS_H

#include <vector>
#include "mesh_partitioner.h"
#include "mesh.h"

/** Class for partitioning the mesh with METIS
*/

class Mesh_Partitioner_METIS : public Mesh_Partitioner {

  public:

    Mesh_Partitioner_METIS(int nparts, Mesh* mesh, int ptype=0);
    virtual ~Mesh_Partitioner_METIS();

    /** Tie mesh nodes in a range.  Only affects the element connectivity
     *  array; tied nodes are not removed from the mesh data structure.
     */
    void build_tie_map(double tol, int start=-1, int end=-1);

    /** Construct disjoint nodal partition of all nodes using METIS
     *  ptype = 0(Recursive), 1(Kway)
     */
    void partition_nodes();
    void set_metis(int ptype) {metis_type=ptype;}

    /** Construct overlapping element partition. Elements can be duplicated
     *  on seperate partitions but are owned by only one parition. This is 
     *  to incorporate branch variables
     */
    void partition_elements();

    /** Construct overlapping global variable partition. Global variables
     *  can be duplicated on seperate partitions but are owned by only one
     *  parition.
     *  FIXME: NOT IMPLEMENTED
     */
    void partition_globals() {}

    /** Assign ids to the nodes, such that nodes on a partition are numbered
     *  consecutively.
     */
    void assign_ids();
    int  get_startid(int i){return IDR[i];}
    int  get_end1id(int i) {return IDR[i+1];}

    /** For a given partition, extract,
     *       - global node number of nodes in partition
     *       - the global reduced node that it is tied to
     *       - the partition that each node is on
     */
    int  get_numnp(int n);
    void get_nodes         (int i, int numnpi, int* nodes);
    void get_tie_map       (int i, int numnpi, int* nodes, int* xmap);
    void get_node_partition(int i, int numnpi, int* nodes, int* parts);
    void get_node_partition(int* parts);
    void get_node_ids      (int n, int numnpn, int* nodes, int* nodeids);

    /** For a given partition, extract,
     *       - global element number of elements in partition
     *       - the partition that each element is on
     *  FIXME: get_branch_ids NOT IMPLEMENTED
     */
    int  get_numelt(int n);
    void get_elements         (int i, int numelti, int* elts);
    void get_element_partition(int i, int numelti, int* elts, int* parts);
    void get_branch_ids       (int i, int numelti, int* elts, int* branchids){}

    /** For a given partition, extract,
     *       - global global variable number of globals in partition
     *       - the partition that each global is on
     *  FIXME: NOT IMPLEMENTED
     */
    int  get_numglobal(int i){}
    void get_globals         (int i, int numglobali, int* globals) {}
    void get_global_partition(int i, int numglobali, int* globals, int* parts) {}
    void get_global_ids      (int i, int numglobali, int* globals, int* globalids) {}

    /** Get info about the underlying mesh
     */
    int get_ndm() {return mesh->get_ndm();}
    int get_ndf() {return mesh->get_ndf();}
    int get_nen() {return mesh->get_nen();}

    /** Get number of independent nodes that exist after tying
     */
    int   get_numrn()  {return numrn;}

    /** Get info about partition
     */
    int&  npart(int i) {return NPART[i];}
    Mesh* get_mesh()   {return mesh;}

  private:

    typedef std::vector<int>  ivec;

    ivec XMAP;                     // Node to Tied reduced Node Mapping
    ivec XMAPR;                    // Tied reduced node to Node Mapping
    ivec NXMAPR;                   // Offset of nodes into XMAPR

    ivec NPART;                    // Node partition index
    ivec EPART;                    // Element partition index
    ivec NEPART;                   // Offset of partitions in EPART
    ivec GPART;                    // Global  partition index
    ivec NGPART;                   // Offset of partitions in GPART

    ivec IDR;                      // ID offset of reduced ID into each partition
                                   // Partition[i] has IDs from IDR[i] to IDR[i]-1

    Mesh* mesh;
    int numrn;                     // Number of reduced nodes
    int metis_type;                // Method of partitioning used in METIS

    /** Build reverse map of the tying
     */
    void build_xmapr();

    /** Construct partition of nodes by constructing adjacency structure on
     *  the reduced nodes
     */
    void build_adj(std::vector<int>& jc, std::vector<int>& ir);
    void adj_sparsity_count(int* counts);
    void adj_build_index(int* jc, int* ir);
    void find_partitions(int* nodes, int n, std::vector<int>& part, int* npart);
};

#endif /* MESH_PARTITIONER_METIS_H */
