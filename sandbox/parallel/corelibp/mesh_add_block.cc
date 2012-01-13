/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <iostream>
#include <cassert>
#include <cstring>
#include "mesh_add_block.h"
#include "shapes.h"
#include "qassemblyp.h"
#include "mesh_partition.h"
#include "coordsorter.h"

using std::vector;
using std::fill;

/*@t ----------
 * \subsection{Mesh_Add_Block construction, destruction, and memory management}
 *
 *@c*/

#define ME Mesh_Add_Block


ME::ME(int ndm, int maxnen, int maxndf) : Mesh(ndm,maxnen,maxndf)
{
}


ME::~ME()
{
}


/*@t ----------
 * \subsection{Block mesh generation}
 * 
 * As each block is added, the id number of the first node,
 * the id number of the first element, the number of nodal points
 * per edge in each dimension, and the order of the element is 
 * recorded.
 *
 *@c*/


void ME::add_block(double* x1, double* x2, int* m,
                     Element* etype, int order)
{

    int ndm = get_ndm();
    int start_node = numnp();
    int start_elem = numelt();

    SNID.push_back(start_node);
    SEID.push_back(start_elem);
    BO.push_back(order);

    for (int i = 0; i < ndm; ++i)
        NNPB.push_back(m[i]);

    // -- Dispatch to each function in Mesh
    if (ndm == 1)
        Mesh::add_block(x1[0], x1[0], m[0], etype, order);

    else if (ndm == 2)
        Mesh::add_block(x1[0], x1[1],
                  x2[0], x2[1],
                  m[0],  m[1], etype, order);

    else if (ndm == 3)
        Mesh::add_block(x1[0], x1[1], x1[2],
                  x2[0], x2[1], x2[2],
                  m[0],  m[1],  m[2], etype, order);

}


void ME::add_block(double x1, double x2, int m,
                   Element* etype, int order)
{
    double x1a[1] = {x1};
    double x2a[1] = {x2};
    int    ma[1]  = {m};
    
    ME::add_block(x1a,x2a,ma,etype,order);
}


void ME::add_block(double x1, double y1,
                   double x2, double y2,
                   int nx, int ny,
                   Element* etype, int order)
{
    double x1a[2] = {x1, y1};
    double x2a[2] = {x2, y2};
    int    ma[2]  = {nx, ny};
    
    ME::add_block(x1a,x2a,ma,etype,order);
}


void ME::add_block(double x1, double y1, double z1,
                     double x2, double y2, double z2,
                     int nx, int ny, int nz,
                     Element* etype, int order)
{
    double x1a[3] = {x1, y1, z1};
    double x2a[3] = {x2, y2, z2};
    int    ma[3]  = {nx, ny, nz};
    
    ME::add_block(x1a,x2a,ma,etype,order);
}


/*@t ----------
 * \subsection{Tied meshes}
 *
 * Same as tie in Mesh.h. This additionally constructs the XMAP array.
 *
 **@c*/


void ME::tie(double tol, int start, int end)
{
    int np = numnp();
    if (start < 0) start = 0;
    if (end   < 0) end   = np;
    int nt = end-start;

    // Set scratch = start:end-1, then sort the indices according
    // to the order of the corresponding nodes.
    //
    CoordSorter sorter(this->get_ndm(), *this, tol);
    vector<int> scratch(nt);
    for (int i = 0; i < nt; ++i) 
        scratch[i] = start+i;
    sort(scratch.begin(), scratch.end(), sorter);

    // XMAP[i] := smallest index belonging to equivalence class of p[i]
    //
    XMAP.resize(np);
    for (int i = 0; i < np; ++i) 
        XMAP[i] = i;

    int inext;
    for (int i = 0; i < nt; i = inext) {
        int class_id = scratch[i];
        for (inext = i+1; inext < nt; ++inext) {
             if (sorter.compare(scratch[i], scratch[inext]) != 0)
                 break;
             if (scratch[inext] < scratch[i])
                class_id = scratch[inext];
        }
        for (int j = i; j < inext; ++j)
            XMAP[scratch[j]] = class_id;
    }

    // Apply map to element connectivity array
    // FIXME: Multiple call, since we cannot access the private 
    //        IX array.
//    Mesh::tie(tol,start,end);
    int nelt = numelt();
    for (int i = 0; i < nelt; ++i) {
        int nen = get_nen(i);
        for (int j = 0; j < nen; ++j)
            if (ix(j,i)>=0)
                ix(j,i) = XMAP[ix(j,i)];
    }
}


/*@T ----------
 * \subsection{Initialization}
 *
 * Mesh initialization occurs after the geometry of the mesh has been
 * established and before any analysis takes place.  This is where we
 * decide how big all the major arrays will be, and set up the index
 * structures that are used for all of our assembly operations.
 *
 * Set [SNID] and [SEID] so that they become offsets in the node numbers
 * into each block, and the difference become the number of node/elements.
 **@c*/


void ME::initialize()
{
    SNID.push_back(numnp());
    SEID.push_back(numelt());
    Mesh::initialize();
}


void ME::initialize_minimal()
{
    SNID.push_back(numnp());
    SEID.push_back(numelt());
    Mesh::initialize_minimal();
}


/*@t -----
 * \subsection{Construct the Node-Element mapping}
 *
 * The mapping between the nodes of this mesh and the elements of the
 * "from" mesh is constructed. This relation can be established from the
 * block generation structure. 
 *@c*/


void ME::node_element_mapping(Mesh_Add_Block* from, int* nemap)
{
    int ndm       = get_ndm();
    int numblocks = get_numblocks();

    assert(ndm==from->get_ndm() && numblocks==from->get_numblocks());

    for (int ib = 0; ib < numblocks; ++ib) {

        int to_order   =       get_order(ib);
        int from_order = from->get_order(ib);
        vector<int> to_nx(ndm+1);      
        vector<int> from_Nx(ndm);    

        // -- to_nx[j+1]: Number of node     per edge in direction j
        //    from_Nx[j]: Number of elements per edge in direction j
        for (int j = 0; j < ndm; ++j) {
            to_nx[j+1] =        get_nx(j,ib);
            from_Nx[j] = (from->get_nx(j,ib)-1)/from_order;
        }

        // -- to_nx[j]  : Offset into array containing info about nodes in direction j
        to_nx[0] = 0;
        for (int j = 0; j < ndm; ++j)
            to_nx[j+1] += to_nx[j];

        // -- Divide indexing of nodes [ix] into baskets corresponding to elements [Ix]
        vector<int> ix(to_nx[ndm]);
        for (int idm = 0; idm < ndm; ++idm) {
            int numn  = to_nx[idm+1] - to_nx[idm];
            double Dx = ((double)(numn-1))/((double)(from_Nx[idm]));

            ix[to_nx[idm]] = 0;
            for (int i = 1; i < numn-1; ++i)
                ix[to_nx[idm]+i] = (int)(ceil(   ((double)i)/Dx)) - 1;
            ix[to_nx[idm+1]-1] = from_Nx[idm]-1;
        }

        // -- Set indices into nemap
        int to_snid   =       get_snid(ib);    // starting node id (to)
        int from_seid = from->get_seid(ib);    // starting element id (from)

        if (ndm==1) {

            int ndivx = to_nx[1]-to_nx[0];
            for (int inx = 0; inx < ndivx; ++inx)
                nemap[to_snid + inx] = from_seid + ix[to_nx[0]+inx];

        } else if (ndm==2) {

            int ndivx = to_nx[1]-to_nx[0];
            int ndivy = to_nx[2]-to_nx[1];
            for (int inx = 0; inx < ndivx; ++inx)
                for (int iny = 0; iny < ndivy; ++iny)
                    nemap[to_snid + ndivy*inx + iny] = from_seid + from_Nx[1]*ix[to_nx[0]+inx] 
                                                                 +            ix[to_nx[1]+iny];

        } else if (ndm==3) {

            int ndivx = to_nx[1]-to_nx[0];
            int ndivy = to_nx[2]-to_nx[1];
            int ndivz = to_nx[3]-to_nx[2];
            for (int inx = 0; inx < ndivx; ++inx)
                for (int iny = 0; iny < ndivy; ++iny)
                    for (int inz = 0; inz < ndivz; ++inz)
                    nemap[to_snid + ndivz*ndivy*inx + ndivz*iny + inz] = from_seid 
                                                   + from_Nx[2]*from_Nx[1]*ix[to_nx[0]+inx] 
                                                   + from_Nx[2]           *ix[to_nx[1]+iny] 
                                                   +                       ix[to_nx[2]+inz];
        }
    }
}

#undef ME


/*@t ----------
 * \section{QAddBlockProlongator}
 *
 * \subsection{QAddBlockProlongator construction, destruction, and memory management}
 *
 * The node(to)-element(from) mapping is copied into a member of the class, NEMAP,
 * through the method [set_nemap].
 *@c*/


#define ME QAddBlockProlongator


ME::ME(Mesh* from, Mesh* to, int* nemap) :
    from(from), to(to)
{
    assert(from->get_ndm()==to->get_ndm() && from->get_ndf()==to->get_ndf());
    if (nemap!=NULL)
        set_nemap(nemap);
}


ME::~ME()
{
}


void ME::set_nemap(int* nemap)
{
    NEMAP.resize(to->numnp());
    memcpy(&NEMAP[0],nemap,to->numnp()*sizeof(int));
}


/*@t ----------
 * \subsection{Assembling the non-zero structure and matrix}
 *
 * The non-zero structure is formed by going through each node
 * in the TO mesh, finding the element which it belongs to in
 * the FROM mesh, and finding the nodes in the FROM to which this
 * element is connect to.
 *
 * To evaluate the prolongator, for each node in the TO mesh, it is
 * evaluated at the nodal shape functions on the FROM mesh. This is
 * conducted by evaluating the node in the TO mesh in the corresponding
 * element in the FROM mesh. A Newton-Raphson iteration is used to solve
 * for the corresponding parent coordinates of the node in the element,
 * then the shape functions are evaluated.
 * 
 *@c*/


void ME::assemble_struct_raw(QStructAssembler* P)
{
    int ndm      = to->get_ndm();
    int ndf      = to->get_ndf();
    int to_numnp = to->numnp();

    // -- Go through nodes in TO mesh
    for (int i = 0; i < to_numnp; ++i)
        if (NEMAP[i]>=0) {

            // -- Find element in FROM mesh
            int eltid    = NEMAP[i];
            Element* elt = from->etype(eltid);
            int from_nen = from->get_nen(eltid);
            vector<int> fid(from_nen);

            // -- Set local id of eltid
            for (int j = 0; j < from_nen; ++j)
                fid[j] = from->ix(j,eltid);

            // -- Add into assembler
            vector<int> from_id(from_nen);
            vector<int> to_id(1);
            for (int j = 0; j < ndf; ++j) {
                to_id[0] = to->inode(j,i);
                for (int k = 0; k < from_nen; ++k)
                    from_id[k] = from->inode(j,fid[k]); 
                P->add(&to_id[0], 1, &from_id[0], from_nen);
            }

        }

}


void ME::assemble_struct(QStructAssembler* P, int reduced)
{
    QGlobalPStructAssembler gassembler(this, *P);
    QReducePStructAssembler rassembler(this, gassembler);
    QStructAssembler* assembler = 
        (reduced ? 
         static_cast<QStructAssembler*>(&rassembler): 
         static_cast<QStructAssembler*>(&gassembler));

    assemble_struct_raw(assembler);
}


void ME::assemble_P_raw(QAssembler* P)
{
    int ndm      = to->get_ndm();
    int ndf      = to->get_ndf();
    int to_numnp = to->numnp();

    // -- Go through nodes in TO mesh
    for (int i = 0; i < to_numnp; ++i)
        if (NEMAP[i]>=0) {

            // -- Find element in FROM mesh
            int eltid    = NEMAP[i];
            Element* elt = from->etype(eltid);
            int from_nen = from->get_nen(eltid);
            vector<int> fid(from_nen);
            vector<int> from_id(from_nen);
            vector<int> to_id(1);
            vector<double> Pe(from_nen);
            vector<double> X(ndm);
            vector<double> x(ndm);

            // -- Set local id of eltid
            for (int j = 0; j < from_nen; ++j)
                fid[j] = from->ix(j,eltid);

            // -- Find X(box coordinate) corresponding to x(spatial coordinate)
            //    and evaluate shape function
            for (int j = 0; j < ndm; ++j) {
                X[j] = 0;
                x[j] = to->x(j,i); 
            }

            if (ndm==1) {

                Quad1d shape(from,eltid,from_nen);
                shape.inv(&x[0],&X[0],5,1e-32);
                shape.eval(&X[0]);
                for (int j = 0; j < from_nen; ++j)
                    Pe[j] = shape.N(j);

            } else if (ndm==2) {

                Quad2d shape(from,eltid,from_nen);
                shape.inv(&x[0],&X[0],5,1e-32);
                shape.eval(&X[0]);
                for (int j = 0; j < from_nen; ++j)
                    Pe[j] = shape.N(j);

            } else if (ndm==3) {

                Quad3d shape(from,eltid,from_nen);
                shape.inv(&x[0],&X[0],5,1e-32);
                shape.eval(&X[0]);
                for (int j = 0; j < from_nen; ++j)
                    Pe[j] = shape.N(j);

            }

            // -- Add into assembler
            for (int j = 0; j < ndf; ++j) {
                to_id[0] = to->inode(j,i);
                for (int k = 0; k < from_nen; ++k)
                    from_id[k] = from->inode(j,fid[k]); 
                P->add(&to_id[0], 1, &from_id[0], from_nen, &Pe[0]);
            }

        } //END:if (NEMAP[i]>=0)

}


void ME::assemble_P(QAssembler* P, int reduced)
{
    QGlobalPAssembler gassembler(this, *P);
    QReducePAssembler rassembler(this, gassembler);
    QAssembler* assembler = 
        (reduced ? 
         static_cast<QAssembler*>(&rassembler) : 
         static_cast<QAssembler*>(&gassembler));

    assemble_P_raw(assembler);
}
#undef ME


/*@t ----------
 * \section{QAddBlockProlongator_Partition}
 *
 * \subsection{QAddBlockProlongator_Partition construction, destruction, and memory management}
 *
 * The node(to)-element(from) mapping defined in terms of global indices, is mapped to local
 * partition indices and copied into NEMAP through the function [set_nemap_global].
 *
 *@c*/


#define ME QAddBlockProlongator_Partition


ME::ME(Mesh_Partition* from, Mesh_Partition* to, int* nemap_global)
 : QAddBlockProlongator(from,to,NULL), from(from), to(to)
{
    if (nemap_global!=NULL)
        set_nemap_global(nemap_global);
}


ME::~ME()
{
}


void ME::set_nemap_global(int* nemap_global)
{

    int pid    = to->get_pid();
    int numnpn = to->numnp();
    vector<int> nemap_local(numnpn);

    for (int i = 0; i < numnpn; ++i)
        if (to->npart(i)==pid)
            nemap_local[i] = from->eglmap(nemap_global[i]);
        else
            nemap_local[i] = -1;

    set_nemap(&nemap_local[0]);
}


/*@t ----------
 * \subsection{Assembling the non-zero structure and matrix}
 *
 * This part differs from the corresponding functions in the base class
 * QAddBlockProlongator only in the Assemblers that are constructed to
 * apply a mapping of the ids to assemble into the matrix. In this case
 * the local indices must be mapped to global indices.
 *
 *@c*/


void ME::assemble_struct(QStructAssembler* P, int reduced)
{
    assert(reduced==1);
    QReducePStructAssemblerP rassembler(this, *P);
    QStructAssembler* assembler = static_cast<QStructAssembler*>(&rassembler);

    assemble_struct_raw(assembler);
}


void ME::assemble_P(QAssembler* P, int reduced)
{
    assert(reduced==1);

    QReducePAssemblerP rassembler(this, *P);
    QAssembler* assembler = static_cast<QAssembler*>(&rassembler);

    assemble_P_raw(assembler);
}
