/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef MESH_ADD_BLOCK_H
#define MESH_ADD_BLOCK_H

#include <vector>

#include "mesh.h"
#include "qassembly.h"

class Mesh_Partition;

/*@t ----------
 * \section{Mesh_Add_Block object}
 *
 * A class derived from the Mesh object to include data structure for 
 * recording the information about the blocks added by add_block. This
 * class is solely for the purpose of constructing geometrical prologators
 * in the multigrid solver. Nodes and elements can only be added through
 * the add_block generators. Otherwise a mapping between meshes cannot be
 * constructed.
 * 
 * FIXME: This class cannot handle branch and global variables.
 *@q*/

/*@t -----------
 * \subsection{Mesh_Add_Block declarations}
 *@c*/

class Mesh_Add_Block : public Mesh {

 public:

    Mesh_Add_Block(int ndm, int maxnen = 0, int maxndf = 0);
    virtual ~Mesh_Add_Block();

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

    /** Assemble mapping between 
     *      Nodes on THIS mesh and the element they belong to on FROM mesh
     *  [nemap] is a [numnp] size array.
     */
    void node_element_mapping(Mesh_Add_Block* from_mesh, int* nemap);

    int get_numblocks()      { return SNID.size()-1;}
    int get_snid(int i)      { return SNID[i];}
    int get_seid(int i)      { return SEID[i];}
    int get_order(int i)     { return BO[i];}
    int get_nx(int i, int j) { return NNPB[j*get_ndm()+i];}
    int map(int i)           { return XMAP[i];}

    /*@q ----------
     * Private data structures are documented in mesh_add_block.cc
     */
 private:
    typedef std::vector<int>           ivec;

    ivec SNID;         // Starting nodal id of the block
    ivec SEID;         // Starting element id of the block
    ivec NNPB;         // Number of nodes per block array(size dof*3 array)
    ivec BO;           // Order of the elements in the block
    ivec XMAP;         // Node to Tied reduced Node Mapping
};

/*@t ----------
 * \section{QAddBlockProlongator object}
 * 
 * Object to construct the prolongator mapping between two meshes,
 * the domain(from mesh) and the range(to mesh). The range is expected to
 * contain more dof than the domain.
 * 
 *@q*/

/*@t -----------
 * \subsection{QAddBlockProlongator declarations}
 *@c*/

class QAddBlockProlongator {

  public:
    QAddBlockProlongator(Mesh* from, Mesh* to, int* nemap=NULL);
    ~QAddBlockProlongator();

    /** Setting the nemap. */
    void set_nemap(int* nemap);

    /** Assemble prolongator matrix structure. */
    void assemble_struct(QStructAssembler* P, int reduced = 1);
    void assemble_struct_raw(QStructAssembler* P);

    /** Assemble global prolongator matrix. */
    void assemble_P(QAssembler* P, int reduced = 1);
    void assemble_P_raw(QAssembler* P);

    Mesh* get_from() {return from;}
    Mesh* get_to()   {return   to;}
    int*  get_nemap(){return &NEMAP[0];}

  private:
    Mesh* from;
    Mesh* to;
    std::vector<int> NEMAP;
};


/*@t ----------
 * \section{QAddBlockProlongator_Partition object}
 * 
 * Object to construct the prolongator mapping between two meshes,
 * the domain(from mesh) and the range(to mesh). The range is expected to
 * contain more dof than the domain. This differs from [QAddBlockProlongator]
 * in that this assume the Mesh to be of type [Mesh_Partition] 
 * 
 * The [nemap_partition] is given in terms of global node numbering and global 
 * element numbering.
 *
 * FIXME: This is a funny class since the only things that require derivation of
 *        this class from QAddBlockProlongator are the functions,
 *        assemble_struct_raw, assemble_P_raw.
 *@q*/


class QAddBlockProlongator_Partition : public QAddBlockProlongator {

  public:
    QAddBlockProlongator_Partition(Mesh_Partition* from, Mesh_Partition* to, int* nemap_global=NULL);
    ~QAddBlockProlongator_Partition();

    /** Setting the nemap. */
    void set_nemap_global(int* nemap_global);

    /** Assemble prolongator matrix structure. */
    void assemble_struct(QStructAssembler* P, int reduced = 1);

    /** Assemble global prolongator matrix. */
    void assemble_P(QAssembler* P, int reduced = 1);

    Mesh_Partition* get_from() {return from;}
    Mesh_Partition* get_to()   {return   to;}

  private:
    Mesh_Partition* from;
    Mesh_Partition* to;
};


#endif /* MESH_ADD_BLOCK_H */
