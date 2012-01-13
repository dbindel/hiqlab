/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef MESH_PARTITION_H
#define MESH_PARTITION_H

#include <vector>
#include <map>

#include "mesh.h"


/*@T -----------
 * \subsection{Accessor functions}
 *
 * The variables in our system can be real or complex, and they can be
 * associated with nodal variables, branch variables, or global variables.
 * Rather than write many tiny accessor functions by hand, we write a
 * few macros to produce these accessor functions.
 *
 *@c*/


#define defq_mesh_accessor(type, name, array) \
   type& name        (int i)        { return array[i];          }


/*@t ----------
 * \section{Partitioned Mesh object}
 *
 * A class derived from the Mesh object to include data structure for
 * a partition of a larger mesh. All data structure from the original Mesh
 * is inherited and enhanced with data required to interact with other partitions.
 * The mesh is constructed by going through exactly the same Lua mesh input file,
 * but through a filter, which selects the nodes, elements, and global variables
 * that should be added to the mesh. 
 *
 *
 * \subsection{Mesh_Partition data structure}
 *
 *
 * \subsubsection{Node mapping between local and global index, partitioning}
 *
 * The total number of nodes on the partition is given by the variable [pnumnp].
 * The mapping between the local partition node id and global node id is given by the NLGMAP.
 * The reverse mapping is given by NGLMAP, in the form of a look up table supplied by
 * the C++ STL "map". The mesh is usually tied, giving rise to redundant nodes. 
 * The mapping between global id to tied global id is given by NGTMAP, and 
 * the mapping between local  id to tied local  id is given by NLTMAP.
 *
 * The ownership of each node is defined by NPART. Nodes with NPART[i]==pid are the only nodes
 * the partition is responsible for keeping consistent.
 *
 *
 * \subsubsection{Element mapping between local and global index, partitioning}
 * 
 * The total number of elements on the partition is given by the variable [pnumelt].
 * The mapping between the local partition element id and global element id is given by the ELGMAP.
 *
 * The ownership of each element is defined by EPART. Elements with EPART[i]==pid are the only elements
 * the partition is responsible for keeping consistent.
 *
 *
 * \subsubsection{Globals mapping between local and global index, partitioning}
 * 
 * The total number of globalss on the partition is given by the variable [pnumglobals].
 * The mapping between the local partition globals id and global globals id is given by the GLGMAP.
 *
 * The ownership of each globals is defined by GPART. Globals with GPART[i]==pid are the only globals
 * the partition is responsible for keeping consistent.
 *
 *
 *@q*/


/*@t -----------
 * \subsection{Mesh_Partition declarations}
 *@c*/

class Mesh_Partition : public Mesh {

  public:
    Mesh_Partition(int pid, int ndm, int maxnen = 0, int maxndf = 0);
    virtual ~Mesh_Partition();

    /** Setting and getting partition info
     */
    void set_pid(int id) {pid = id;};
    int  get_pid()       {return pid;}

    /** Set the global number of the nodes(nps), elements(elts) that 
     *  are assigned to this partition. Also set the global node to
     *  tied node map(tmap) and global node and element partition map(pmap)
     *
     *  
     *  If [i] is the local node number,
     *
     *  nodes[i]: global node number
     *  tmap [i]: the global node that it is tied to
     *  pmap [i]: the partition number it is assigned to
     *
     *
     *  If [i] is the local element number
     *
     *  elts[i]: the global element number
     *  pmap[i]: the partition number it is assigned to
     *
     *
     *  If [i] is the local globals number
     *
     *  globals[i]: the global globals number
     *  pmap   [i]: the partition number it is assigned to
     *
     */
    void set_nodes    (int numnp,     int* nodes);
    void set_ntie_map (int numnp,     int* tmap);
    void set_npart_map(int numnp,     int* pmap);

    void set_elements (int numelt,    int* elts );
    void set_epart_map(int numelt,    int* pmap);

    void set_globals  (int numglobals,int* globals);
    void set_gpart_map(int numglobals,int* pmap);

    /** Add element(s) to the mesh directly.
     *  In each case, return the identifier of the first item created.
     */
    int add_node(double* x, int n=1);

    /** Add element(s) to the mesh directly.
     *  In each case, return the identifier of the first item created.
     *
     *  add_element      : e is given in global node indices
     *  add_element_local: e is given in local  node indices
     */
    int add_element      (int* e, Element* etype, int nen, int n=1);
    int add_element_local(int* e, Element* etype, int nen, int n=1);

    /** Add global to the mesh directly
     *  In each case, return the identifier of the first item created.
     *  //FIXME: NOT IMPLEMENTED YET
     */
    int add_global(int n) {}

    /** Construct a Cartesian block. */
    void add_block(double* x1, double* x2, int* m,
                   Element* etype, int order=1);
    void add_block(double x1, double x2, int m,
                   Element* etype, int order=1);
    void add_block(double x1, double y1,
                   double x2, double y2,
                   int nx, int ny,
                   Element* etype, int order=1);
    void add_block(double x1, double y1, double z1,
                   double x2, double y2, double z2,
                   int nx, int ny, int nz,
                   Element* etype, int order=1);

    /** Tie mesh nodes in a range.  Only affects the element connectivity
     *  array; tied nodes are not removed from the mesh data structure.
     */
    void tie(double tol, int start=-1, int end=-1);

    /** Allocate space for global arrays and assign initial index map. */
    void initialize();
    void initialize_minimal();

    /** Assigning of the ids is done in the following steps.
     *
     *  0. Call Mesh::assign_ids and set local reduced id array
     *  1. Allocate memory to hold global reduced id array in [IDG]
     *  2. Set the starting id(start) and one after the last id(end1)
     *  3. Set node ids
     *  //FIXME: NOT IMPLEMENTED YET
     *  4. Set branch ids
     *  5. Set globals ids
     *
     */
    void set_ids_range(int start, int end1) {pstartid=start;pend1id=end1;pnumid=end1-start;}
    void set_node_ids   (int numnp, int ndf, int* nodeids);

    /** Assemble local portion of the tangent matrix structure
     */
    void assemble_struct(QStructAssembler* K, int reduced=1);

    /** Assemble local portion of the stiffness and mass matrices
     */
    void assemble_dR(QAssembler* K, double cx, double cv, double ca, int reduced=1);

    /** Starting and ending reduced id info of this partition
     */
    int get_pstartid() {return pstartid;}
    int get_pend1id()  {return pend1id;}
    int get_pnumid()   {return pnumid;}

    /** Get and set the IDG array */
    defq_mesh_accessor( int, idg, IDG)

    int& nlgmap(int i) {return NLGMAP[i];}
    int& elgmap(int i) {return ELGMAP[i];}
    int& glgmap(int i) {return GLGMAP[i];}
    int& npart (int i) {return NPART[i];}
    int& epart (int i) {return EPART[i];}
    int& gpart (int i) {return GPART[i];}

    int  nglmap(int i) {return NGLMAP[i];}
    int  eglmap(int i) {return EGLMAP[i];}

  private:

    typedef std::vector<int> ivec;
    typedef std::map<int,int> imap;

    int  pid;          // Global partition id
    int  pnumnp;       // Partition number of nodal points
    int  pnumelt;      // Partition number of elements
    int  pnumglobals;  // Partition number of globals

    ivec NLGMAP;       // Local node to Global node mapping
    imap NGLMAP;       // Global node to Local node mapping
    ivec NGTMAP;       // Node global to reduced node tie mapping 
    ivec NLTMAP;       // Node local  to reduced node tie mapping 
    ivec ELGMAP;       // Local element to Global element mapping
    imap EGLMAP;       // Global element to Local element mapping
    ivec GLGMAP;       // Local globals to Global globals mapping

    ivec NPART;        // Nodal   partitioning
    ivec EPART;        // Element partitioning
    ivec GPART;        // Globals partitioning

    int  pnumid;       // Number of ids in this partition
    int  pstartid;     // Starting id of this partition
    int  pend1id;      // Last+1 id of this partition

    ivec IDG;          // ID array globaly

    int n_counter;     // Counter for global number of nodes
    int e_counter;     // Counter for global number of elements
    int g_counter;     // Counter for global number of globals

    void build_local_tie_map();
    void set_partition_neighbors();

};

/*@q -------------------
 * End macro definitions
 */
#undef defq_mesh_accessor


#endif /* MESH_PARTITION_H */
