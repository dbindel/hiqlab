/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <iostream>
#include <cassert>

#include "mesh_partition.h"
#include "luasupport.h"
#include "tolua++.h"
#include "qassemblyp.h"

using std::vector;
using std::fill;
using std::copy;


/*@t ----------
 * \subsection{Mesh_Partition construction, destruction, and memory management}
 *
 *@c*/


#define ME Mesh_Partition


ME::ME(int pid, int ndm, int maxnen, int maxndf) : Mesh(ndm,maxnen,maxndf), pid(pid),
  n_counter(0), e_counter(0), g_counter(0), pnumnp(0), pnumelt(0), pnumglobals(0),
                                            pnumid(0), pstartid(0), pend1id(0)
{
}


ME::~ME()
{
}


void ME::set_nodes( int numnp,  int* nps)
{
    pnumnp = numnp;

    // -- Construct NLGMAP
    NLGMAP.resize(numnp);
    copy(nps, nps+numnp, NLGMAP.begin());

    // -- Construct NGLMAP
    for (int i = 0; i < numnp; ++i)
        NGLMAP.insert(std::pair<int,int>(NLGMAP[i],i));
}


void ME::set_ntie_map(int numnp, int* tmap)
{
    assert(pnumnp==numnp);
    NGTMAP.resize(numnp);
    copy(tmap, tmap+numnp, NGTMAP.begin());
}


void ME::build_local_tie_map()
{
    NLTMAP.resize(pnumnp);
    for (int i = 0; i < pnumnp; i++)
        NLTMAP[i] = NGLMAP[NGTMAP[i]];
}


void ME::set_npart_map(int numnp, int* pmap)
{
    assert(pnumnp==numnp);
    NPART.resize(numnp);
    copy(pmap, pmap+numnp, NPART.begin());
}


void ME::set_elements(int numelt, int* elts)
{
    pnumelt = numelt;

    // -- Construct ELGMAP
    ELGMAP.resize(numelt);
    copy(elts, elts+numelt, ELGMAP.begin());

    // -- Construct EGLMAP
    for (int i = 0; i < numelt; ++i)
        EGLMAP.insert(std::pair<int,int>(ELGMAP[i],i));
}


void ME::set_epart_map(int numelt, int* pmap)
{
    assert(pnumelt==numelt);
    EPART.resize(numelt);
    copy(pmap, pmap+numelt, EPART.begin());
}


void ME::set_globals(int numglobals, int* globals)
{
    pnumglobals = numglobals;
    GLGMAP.resize(numglobals);
    copy(globals, globals+numglobals, GLGMAP.begin());
}


void ME::set_gpart_map(int numglobals, int* pmap)
{
    assert(pnumglobals==numglobals);
    GPART.resize(numglobals);
    copy(pmap, pmap+numglobals, GPART.begin());
}


/*@t -----
 * \subsection{Basic mesh generation}
 *
 * Since it is assumed that the same Lua mesh input file will be run
 * through one must filter and accept only the nodes, elements, and globals
 * that are to be on this partition. The global number of the nodes, elements,
 * and globals are recorded in the counters,[n_counter,e_counter, g_counter].
 *
 *@c*/


int ME::add_node(double* x, int n)
{
    int ndm = get_ndm();
    int id  = numnp();
    int n_counter_l = id;
    vector<double> addn(n*ndm);

    int ind = 0;
    for (int i = 0; i < n; ++i) {

        // -- add node if the global number of this node
        //                [n_counter] 
        //    is equal to the global node number
        //                NLGMAP[n_counter_l]
        //    of the next local node number [n_counter_l]
        //    
        if(n_counter_l<pnumnp && n_counter==NLGMAP[n_counter_l]) {
            for (int j = 0; j < ndm; ++j)
                addn[ind++]=x[i*ndm+j];
            n_counter_l++;
        }

        // -- update number of global nodes added
        n_counter++;
    }

    Mesh::add_node(&addn[0],n_counter_l-id);

    return id;
}


int ME::add_element(int* e, Element* etype, int nen, int n)
{
    int id  = numelt();
    int e_counter_l = id;
    vector<int> adde(n*nen);

    int ind = 0;
    for (int i = 0; i < n; ++i) {

        if(e_counter_l<pnumelt && e_counter==ELGMAP[e_counter_l]) {
            for (int j = 0; j < nen; ++j)
                adde[ind++] = e[i*nen+j];
            e_counter_l++;
        }

        // -- update number of global elements added
        e_counter++;
    }

    // -- Must filter the e indices through GLMAP
    for (int i = 0; i < (e_counter_l-id); ++i)
        for (int j = 0; j < nen; ++j)
            adde[i*nen+j] = NGLMAP[adde[i*nen+j]];

    Mesh::add_element(&adde[0],etype, nen, e_counter_l-id);

    return id;
}


int ME::add_element_local(int* e, Element* etype, int nen, int n)
{
    int id = numelt();
    int e_counter_l = id;
    vector<int> adde(n*nen);

    int ind = 0;
    for (int i = 0; i < n; ++i) {

        if(e_counter_l<pnumelt && e_counter==ELGMAP[e_counter_l]) {
            for (int j = 0; j < nen; ++j)
                adde[ind++] = e[i*nen+j];
            e_counter_l++;
        }

        // -- update number of global elements added
        e_counter++;
    }

    Mesh::add_element(&adde[0],etype, nen, e_counter_l-id);

    return id;
}


/*@t ----------
 * \subsection{Block mesh generation}
 *
 * The block mesh generators produce rectilinear coordinate-aligned
 * block meshes, which can subsequently be turned into more
 * interesting shape by mapping and tying operations.  For each
 * dimension, we specify the coordinate range ($[x_1, x_2]$) and the
 * number of {\em nodes} (not elements) in that direction.  It is a
 * checked error to try to specify a number of nodes that doesn't
 * yield an integer number of elements.
 *
 * For this Mesh_Partition, the nodes and elements must be filtered.
 *
 **@c*/


void ME::add_block(double* x1, double* x2, int* m,
                     Element* etype, int order)
{
    int ndm = get_ndm();

    if (ndm == 1)
        ME::add_block(x1[0], x2[0], m[0], etype, order);

    else if (ndm == 2)
        ME::add_block(x1[0], x1[1],
                      x2[0], x2[1],
                       m[0],  m[1], etype, order);

    else if (ndm == 3)
        ME::add_block(x1[0], x1[1], x1[2],
                      x2[0], x2[1], x2[2],
                       m[0],  m[1],  m[2], etype, order);
}


void ME::add_block(double x1, double x2, int nx,
                     Element* etype, int order)
{
    assert( nx > 0 && (nx-1)%order == 0 );

    int start_node = numnp();

    // -- construct local mapping between full nodes and existing local nodes
    vector<int> lmap(nx);
    fill(lmap.begin(),lmap.end(),-1);
    int ncl = start_node;
    int nc  = n_counter;

    for (int i = 0; i < nx; ++i) {

        // -- add node if the global number of this node
        //                [n_counter] 
        //    is equal to the global node number
        //                NLGMAP[n_counter_l]
        //    of the next local node number [n_counter_l]
        //    
        if (ncl<pnumnp && nc==NLGMAP[ncl])
            lmap[i] = ncl++;

        // -- update number of global nodes added
        nc++;
    }

    // -- add nodes
    for (int ix = 0; ix < nx; ++ix) {
        double x[1];
        x[0] = ((nx-1-ix)*x1 + ix*x2) / (nx-1);
        add_node(x);
    }

    // -- add elements
    vector<int> adde(order+1);
    for (int ix = 0; ix < nx-order; ix += order) {

        for (int jx = 0; jx < order+1; ++jx)
            adde[jx] = lmap[(ix+jx)];

        add_element_local(&adde[0], etype, order+1, 1);
    }
}


void ME::add_block(double x1, double y1,
                     double x2, double y2,
                     int nx, int ny,
                     Element* etype, int order)
{
    assert( nx > 0 && (nx-1)%order == 0 );
    assert( ny > 0 && (ny-1)%order == 0 );

    int start_node = numnp();

    // -- construct local mapping between full nodes and existing local nodes
    vector<int> lmap(nx*ny);
    fill(lmap.begin(),lmap.end(),-1);
    int ncl = start_node;
    int nc  = n_counter;

    for (int i = 0; i < nx*ny; ++i) {

        // -- add node if the global number of this node
        //                [n_counter] 
        //    is equal to the global node number
        //                NLGMAP[n_counter_l]
        //    of the next local node number [n_counter_l]
        //    
        if (ncl<pnumnp && nc==NLGMAP[ncl])
            lmap[i] = ncl++;

        // -- update number of global nodes added
        nc++;
    }

    // -- add nodes
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            double x[2];
            x[0] = ((nx-1-ix)*x1 + ix*x2) / (nx-1);
            x[1] = ((ny-1-iy)*y1 + iy*y2) / (ny-1);
            add_node(x);
        }
    }

    // -- add elements
    vector<int> adde((order+1)*(order+1));
    for (int ix = 0; ix < nx-order; ix += order)
        for (int iy = 0; iy < ny-order; iy += order) {
            for (int jx = 0; jx < order+1; ++jx)
                for (int jy = 0; jy < order+1; ++jy)
                    adde[jx*(order+1)+jy] = lmap[(iy+jy) + ny*(ix+jx)];
            add_element_local(&adde[0], etype, (order+1)*(order+1), 1);
        }
}


void ME::add_block(double x1, double y1, double z1,
                     double x2, double y2, double z2,
                     int nx, int ny, int nz,
                     Element* etype, int order)
{
    assert( nx > 0 && (nx-1)%order == 0 );
    assert( ny > 0 && (ny-1)%order == 0 );
    assert( nz > 0 && (nz-1)%order == 0 );

    int start_node = numnp();

    // -- construct local mapping between full nodes and existing local nodes
    vector<int> lmap(nx*ny*nz);
    fill(lmap.begin(),lmap.end(),-1);
    int ncl = start_node;
    int nc  = n_counter;

    for (int i = 0; i < nx*ny*nz; ++i) {

        // -- add node if the global number of this node
        //                [n_counter] 
        //    is equal to the global node number
        //                NLGMAP[n_counter_l]
        //    of the next local node number [n_counter_l]
        //    
        if (ncl<pnumnp && nc==NLGMAP[ncl])
            lmap[i] = ncl++;

        // -- update number of global nodes added
        nc++;
    }

    // -- add nodes
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                double x[3];
                x[0] = ((nx-1-ix)*x1 + ix*x2) / (nx-1);
                x[1] = ((ny-1-iy)*y1 + iy*y2) / (ny-1);
                x[2] = ((nz-1-iz)*z1 + iz*z2) / (nz-1);
                add_node(x);
            }
        }
    }

    // -- add elements
    vector<int> adde((order+1)*(order+1)*(order+1));
    for (int ix = 0; ix < nx-order; ix += order)
        for (int iy = 0; iy < ny-order; iy += order)
            for (int iz = 0; iz < nz-order; iz += order) {
                for (int jx = 0; jx < order+1; ++jx)
                    for (int jy = 0; jy < order+1; ++jy)
                        for (int jz = 0; jz < order+1; ++jz)
                            adde[jx*(order+1)*(order+1)+jy*(order+1)+jz] = 
                              lmap[(iz+jz) + nz*((iy+jy) + ny*(ix+jx))];
                add_element_local(&adde[0], etype, (order+1)*(order+1)*(order+1), 1);
            }
}


/*@t -----
 * \subsection{Tie meshes}
 *
 * The mesh is tied according to the tie map given. First the
 * tie map between the local nodes and reduced nodes is constructed.
 * Then this is applied to the connectivity array.
 *@c*/


void ME::tie(double tol, int start, int end)
{
    // -- Build local tie map
    build_local_tie_map();

    // -- tie according to local tie map
    for (int i = 0; i < pnumelt; ++i) {
        int nen = get_nen(i);
        for (int j = 0; j < nen; ++j)
            ix(j,i) = NLTMAP[ix(j,i)];
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
 **@c*/


void ME::initialize()
{
    Mesh::initialize();

    // -- Initialize IDG array to -1
    int M = get_idsize();
    IDG.resize(M);
    fill(IDG.begin(),IDG.end(),-1);
}


void ME::initialize_minimal()
{
    Mesh::initialize_minimal();

    // -- Initialize IDG array to -1
    int M = get_idsize();
    IDG.resize(M);
    fill(IDG.begin(),IDG.end(),-1);
}


/*@t -----
 * \subsection{Assigining reduced indices}
 *
 * The reduced indices are set explicitly through arrays
 *
 */


void ME::set_node_ids(int numnp, int ndf, int* nodeids)
{
    assert(ndf == get_ndf() && numnp == pnumnp);

    for (int i = 0; i < numnp; ++i)
        for (int j = 0; j < ndf; ++j)
                IDG[inode(j,i)] = nodeids[inode(j,i)];
}


/*@t -----
 * \subsection{Assembly loops}
 *
 * Assembles the local portion of the structure and matrix
 *
 *@c*/


void ME::assemble_dR(QAssembler* K, double cx, double cv, double ca,
                       int reduced)
{
    assert(reduced==1);

    QReduceAssemblerP rassembler(this, *K);
    QAssembler* assembler = static_cast<QAssembler*>(&rassembler);

    Mesh::assemble_dR_raw(assembler,cx,cv,ca);
}


void ME::assemble_struct(QStructAssembler* K, int reduced)
{
    assert(reduced==1);

    QReduceStructAssemblerP rassembler(this, *K);
    QStructAssembler* assembler = static_cast<QStructAssembler*>(&rassembler);

    Mesh::assemble_struct_raw(assembler);
}

