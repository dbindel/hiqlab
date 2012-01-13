#ifndef MESH_PARTITIONER_H 
#define MESH_PARTITIONER_H

/** Interface for partitioning the mesh
*/
class Mesh_Partitioner {

  public:

    Mesh_Partitioner(int nparts) : num_partitions(nparts) {}
    virtual ~Mesh_Partitioner() {}

    /** Tie mesh nodes in a range and build tie map.  Only affects the 
     *  element connectivity array; tied nodes are not removed from the 
     *  mesh data structure.
     */
    virtual void build_tie_map(double tol, int start=-1, int end=-1)=0;

    /** Construct nodal partition of all nodes
     */
    virtual void partition_nodes()=0;

    /** Construct overlapping element partition. Elements can be duplicated
     *  on seperate partitions but are owned by only one parition. This is 
     *  to incorporate branch variables
     */
    virtual void partition_elements()=0;

    /** Construct overlapping global variable partition. Global variables 
     *  can be duplicated on seperate partitions but are owned by only one 
     *  parition.
     */
    virtual void partition_globals()=0;

    /** Assign ids to the nodes, such that nodes on a partition are numbered
     *  consecutively.
     */
    virtual void assign_ids()=0;
    virtual int  get_startid(int i)=0;
    virtual int  get_end1id(int i)=0;

    /** For a given partition, extract,
     *       - global node number of nodes in partition
     *       - the global reduced node that it is tied to
     *       - the partition that each node is on
     */
    virtual int  get_numnp(int i)=0;
    virtual void get_nodes         (int i, int numnpi, int* nodes)=0;
    virtual void get_tie_map       (int i, int numnpi, int* nodes, int* xmap)=0;
    virtual void get_node_partition(int i, int numnpi, int* nodes, int* parts)=0;
    virtual void get_node_ids      (int i, int numnpi, int* nodes, int* nodeids)=0;

    /** For a given partition, extract,
     *       - global element number of elements in partition
     *       - the partition that each element is on
     */
    virtual int  get_numelt(int i)=0;
    virtual void get_elements         (int i, int numelti, int* elts)=0;
    virtual void get_element_partition(int i, int numelti, int* elts, int* parts)=0;
    virtual void get_branch_ids       (int i, int numelti, int* elts, int* branchids)=0;

    /** For a given partition, extract,
     *       - global global variable number of globals in partition
     *       - the partition that each global is on
     */
    virtual int  get_numglobal(int i)=0;
    virtual void get_globals         (int i, int numglobali, int* globals)=0;
    virtual void get_global_partition(int i, int numglobali, int* globals, int* parts)=0;
    virtual void get_global_ids      (int i, int numglobali, int* globals, int* globalids)=0;

    /** Get info about the underlying mesh
     */
    virtual int   get_ndm()=0;
    virtual int   get_ndf()=0;
    virtual int   get_nen()=0;
    int           get_num_partitions() {return num_partitions;};

  protected:

    int num_partitions;               // Number of partitions
};

#endif /* MESH_PARTITIONER_H */
