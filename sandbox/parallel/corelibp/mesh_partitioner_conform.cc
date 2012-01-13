/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <iostream>
#include <algorithm>
#include <cassert>

#include "mesh_partitioner_conform.h"
#include "mesh.h"
#include "coordsorter.h"

using std::vector;

#define ME Mesh_Partitioner_Conform


/* -----
 * Define methods to evaluate order and equalities of PEpair
 */


static bool PEPair_eq(const ME::PEPair& x, const ME::PEPair& y)
{
    if (x.p==y.p && x.e==y.e) return true;
    else                      return false;
}


static bool PEPair_plt(const ME::PEPair& x, const ME::PEPair& y)
{
    if (x.p < y.p) return true;
    if (x.p > y.p) return false;
    if (x.e < y.e) return true;
    return false;
}


static bool PEPair_elt(const ME::PEPair& x, const ME::PEPair& y)
{
    if (x.e < y.e) return true;
    if (x.e > y.e) return false;
    if (x.p < y.p) return true;
    return false;
}


/*@t -----
 * \subsection{Mesh_Partitioner_Conform construction, destruction, and memory management}
 *c@*/


ME::ME(int nparts, Mesh* mesh) : Mesh_Partitioner(nparts), 
				 mesh(mesh), to_numnp(0), TO_NPART(0), NEMAP(0), own_to_data(0)
{
}


ME::~ME()
{
    if (own_to_data) {
        delete[] TO_NPART;
        delete[] NEMAP;
    }
}


/*@t --------
 * \subsection{Building the tie map}
 *
 * The array XMAP is constructed that gives the mapping between
 * the full node set and the node that it is tied to, i.e., the smallest
 * index belonging to the equivalence class of that node. This tie map
 * must coincide with the tie that is conducted in the mesh generation.
 * Thus [tol] must be the same as for the call to Mesh:tie(). Failure to
 * to do so will end in unexpected behavior.
 *
 * The function build_xmapr constructs the reverse map of XMAP into
 * XMAPR with NXMAPR as an offset into this array.
 *@c*/


void ME::build_tie_map(double tol, int start, int end)
{
    int np = mesh->numnp();
    if (start < 0) start = 0;
    if (end   < 0) end   = np;
    int nt = end-start;

    // Set scratch = start:end-1, then sort the indices according
    // to the order of the corresponding nodes.
    //
    CoordSorter sorter(mesh->get_ndm(), *mesh, tol);
    vector<int> scratch(nt);
    for (int i = 0; i < nt; ++i) 
        scratch[i] = start+i;
    sort(scratch.begin(), scratch.end(), sorter);

    // XMAP[i]  := smallest index belonging to equivalence class of p[i]
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

    // -- build reverse map 
    build_xmapr();
}


void ME::build_xmapr()
{
    int numnp = mesh->numnp();
    NXMAPR.clear();
    NXMAPR.resize(numnp+1);

    // -- NXMAPR[i] := number of nodes for reduced node i
    for (int i = 0; i < numnp; ++i)
        NXMAPR[XMAP[i]]++;

    // -- NXMAPR[i] := number of nodes for reduced node 0 .. i
    for (int i = 1; i < numnp+1; ++i)
        NXMAPR[i] += NXMAPR[i-1];

    // -- Build map via bucket sort
    XMAPR.clear();
    XMAPR.resize(numnp);
    for (int i = 0; i < numnp; ++i) {
        int node = XMAP[i];
        int slot = --NXMAPR[node];
        XMAPR[slot] = i;
    }
}


/*@t -----
 * \subsection{Partitioning nodes for conformity}
 *
 * The mesh is partitioned to conform with a given element partitioning.
 * Information is given through the two arrays [to_npart] and [nemap] each
 * of length [to_numnp] which is the number of nodes on the mesh that is the 
 * range of the mapping from this mesh. The entries in [to_npart] defines which
 * partition the node belongs to and the entries in [nemap] define the element 
 * on this mesh to which this node belongs.
 *
 * Node ownership between the partitions is defined in the following order. 
 * 1. Assign ownership of nodes connected to an element belonging to only one partition
 * 2. Assign ownership of nodes connected to an element belonging to more that on partition
 * 3. Assert that all nodes have been assigned ownership
 *
 *@c*/


void ME::build_pepairs(int to_np, int* to_npart, int* nemap, int own_data)
{
    // -- Store pointers
    to_numnp    = to_np;
    TO_NPART    = to_npart;
    NEMAP       = nemap;
    own_to_data = own_data;


    pepairs.resize(to_numnp);

    for (int i = 0; i < to_numnp; ++i) {
        pepairs[i].p = to_npart[i];
        pepairs[i].e = nemap[i];
    }

    // -- Sort and remove multiple copies of the same pair
    //    First sort by nodes and remove
    sort(pepairs.begin(),pepairs.end(),PEPair_plt);
    vector<PEPair>::iterator p = unique(pepairs.begin(),pepairs.end(),PEPair_eq);
    pepairs.erase(p,pepairs.end());

    //    Next sort by elements and remove
    sort(pepairs.begin(),pepairs.end(),PEPair_elt);
    p = unique(pepairs.begin(),pepairs.end(),PEPair_eq);
    pepairs.erase(p,pepairs.end());

}


void ME::partition_nodes()
{
    NPART.resize(mesh->numnp());
    fill(NPART.begin(), NPART.end(), -1);

    // -- npepairs[i] : Offset into pepairs
    int nesize = pepairs.size();
    npepairs.resize(mesh->numelt()+1);
    memset(&npepairs[0], 0, (mesh->numelt()+1)*sizeof(int));
    for (int i = 0; i < nesize; ++i)
        npepairs[pepairs[i].e+1]++;
    for (int i = 0; i < mesh->numelt(); ++i)
        npepairs[i+1]+=npepairs[i];

    // -- Assign partition to nodes connected to elements
    //    owned by one partition. For overlapping nodes, smallest
    //    number partition has ownership
    for (int i = 0; i < mesh->numelt(); ++i) {
        int npe = npepairs[i+1]-npepairs[i];
        if (npe==1) {
            int eltn = pepairs[npepairs[i]].e;
            int nen  = mesh->get_nen(eltn);
            for (int j = 0; j < nen; ++j)
                if (NPART[mesh->ix(j,eltn)]<0)
                    NPART[mesh->ix(j,eltn)] = pepairs[npepairs[i]].p;
        }
    }

    // -- Assign partition to nodes connected to elements
    //    owned by more than one partition. For overlapping nodes, smallest
    //    number partition has ownership
    for (int i = 0; i < mesh->numelt(); ++i) {
        int npe = npepairs[i+1]-npepairs[i];
        if (npe>1) {
            int eltn = pepairs[npepairs[i]].e;
            int nen  = mesh->get_nen(eltn);
            for (int j = 0; j < nen; ++j)
                if (NPART[mesh->ix(j,eltn)]<0)
                    NPART[mesh->ix(j,eltn)] = pepairs[npepairs[i]].p;
        }
    }

    // -- All nodes that are tied must also be partitioned
    //    Also assert if all nodes have been partitioned
    for (int i = 0; i < mesh->numnp(); ++i) {
        NPART[i] = NPART[XMAP[i]];
        assert(NPART[i]>=0);
    }
}


/*@t -----
 * \subsection{Partitioning elements for conformity}
 *
 * For constructing stiffness matrices from each partition
 * it is more convenient to have overlapping element partitions,
 * since this allows one to construct the slice of the stiffness
 * matrix locally without any message passing.
 *
 * This process is done by going through all elements and assigning
 * the partition id to each element. EPART keeps info on the partition
 * it is assigned to, and NEPART is an offset into this array.
 *
 * Other than the elements that are attached to the nodes that
 * are owned by the partition, the elements that are required for
 * prolongation are added in the overlapping element partition.
 * 
 *@c*/


void ME::partition_elements()
{
    int numelt = mesh->numelt();
    NEPART.clear();
    NEPART.push_back(0);

    for (int i = 0; i < numelt; ++i) {

        int nen = mesh->get_nen(i);
        int npart;
        vector<int> part(nen);
        find_partitions( &(mesh->ix(0,i)), nen, part, &npart);

        // -- Must add elements required for prolongation
        vector<int> partp(npepairs[i+1]-npepairs[i]);
        for (int j = 0; j < (npepairs[i+1]-npepairs[i]); ++j)
            partp[j] = pepairs[npepairs[i]+j].p;

        // -- Take set union of part and partp into partf
        vector<int> partf;
        set_union(part.begin(),part.end(),partp.begin(),partp.end(),back_inserter(partf));
        npart = partf.size();

        NEPART.push_back(npart);
        for (int j = 0; j < npart; ++j)
            EPART.push_back(partf[j]);
    }

    for (int i = 1; i < numelt+1; ++i)
        NEPART[i] += NEPART[i-1];

}


void ME::find_partitions(int* nodes, int n, vector<int>& part, int* npart)
{
    for (int i = 0; i < n; ++i)
        part[i] = NPART[nodes[i]];

    sort(part.begin(),part.end());
    vector<int>::iterator p = unique(part.begin(),part.end());
    part.erase(p,part.end());
    *npart = part.size();
}


/*@t -----
 * \subsection{Assigning ids consecutively on each partition}
 *
 * For construction of matrices it is convenient that the ids of the matrix
 * are consecutive and coincide with the partitioning of the mesh. Thus
 * a relabeling of the ids is conducted so that, partition [i] has all ids
 * from IDR[i] to IDR[i+1]-1. Within each partition, the ids are labeled
 * in order of nodeids, branchids, and globalids.
 *
 *@c*/


void ME::assign_ids()
{
    IDR.clear();
    IDR.resize(num_partitions+1);

    // -- For the nodeids count how many belong to each partition
    int numnp = mesh->numnp();
    int ndf   = mesh->get_ndf();
    for (int i = 0; i < numnp; ++i)
        for (int j = 0; j < ndf; ++j)
            if (mesh->id(j,i)>=0)
                IDR[NPART[i]+1]++;

    // -- For each element count how many branchids belong to each partition
    int numelt = mesh->numelt();
    for (int i = 0; i < numelt; ++i) {
        int nbranchid = mesh->nbranch_id(i);
        for (int j = 0; j < nbranchid; ++j)
            if (mesh->branchid(j,i)>=0)
                IDR[EPART[i]+1]++;
    }

    // -- For each element count how many globalids belong to each partition
    //FIXME: must implement

    // -- accumulate IDR so that it becomes offset
    vector<int> numidsp(num_partitions);
    for (int i = 0; i < num_partitions; ++i) {
        IDR[i+1] += IDR[i];
        numidsp[i] = IDR[i];
    }

    // -- Assign nodeids
    for (int i = 0; i < numnp; ++i)
        for (int j = 0; j < ndf; ++j)
            if (mesh->id(j,i)>=0)
                mesh->id(j,i) = numidsp[NPART[i]]++;

    // -- Assign branchids
    numelt = mesh->numelt();
    for (int i = 0; i < numelt; ++i) {
        int nbranchid = mesh->nbranch_id(i);
        for (int j = 0; j < nbranchid; ++j)
            if (mesh->branchid(j,i)>=0)
                mesh->branchid(j,i) = numidsp[EPART[i]]++;
    }

    // -- Assign gids
    //FIXME: must implement

    // -- Check to see if all ids have been reassigned
    for (int i = 0; i < num_partitions; ++i)
        assert(numidsp[i]==IDR[i+1]);
}


/*@t -----
 * \subsection{Extracting nodal information for each partition}
 *
 *
 *@c*/


int ME::get_numnp(int n)
{
    int numnp = mesh->numnp();
    int numnpn= 0;

    // -- Number of nodal points that are owned
    for (int i = 0; i < numnp; ++i)
        if (NPART[i]==n)
            numnpn++;

    // -- Number of nodal points that are not owned
    int numelt = mesh->numelt();
    vector<int> offpart;

    for (int i = 0; i < numelt; ++i)
        for (int j = NEPART[i]; j < NEPART[i+1]; ++j)
            if (EPART[j]==n) {
                int nen = mesh->get_nen(i);
                for (int k = 0; k < nen; ++k)
                    if(NPART[mesh->ix(k,i)]!=n)
                        offpart.push_back(mesh->ix(k,i)); 
            }

    // -- filter out duplicate not owned nodes
    sort(offpart.begin(),offpart.end());
    vector<int>::iterator p = unique(offpart.begin(),offpart.end());

    // -- count ones that are also tied
    for (vector<int>::iterator ip = offpart.begin(); ip!=p; ++ip)
        numnpn += (NXMAPR[*ip+1]-NXMAPR[*ip]);

   
    return numnpn;
}


void ME::get_nodes(int n, int numnpn, int* nodes)
{

    int numnp = mesh->numnp();
    vector<int> nodesp;

    // -- add nodes that are owned
    for (int i = 0; i < numnp; ++i)
        if (NPART[i]==n)
            nodesp.push_back(i);

    // -- Number of nodal points that are not owned
    int numelt = mesh->numelt();
    vector<int> offpart;

    for (int i = 0; i < numelt; ++i)
        for (int j = NEPART[i]; j < NEPART[i+1]; ++j)
            if (EPART[j]==n) {
                int nen = mesh->get_nen(i);
                for (int k = 0; k < nen; ++k)
                    if(NPART[mesh->ix(k,i)]!=n)
                        offpart.push_back(mesh->ix(k,i)); 
            }

    // -- filter out duplicate not owned nodes
    sort(offpart.begin(),offpart.end());
    vector<int>::iterator p = unique(offpart.begin(),offpart.end());

    // -- add ones that are not owned and also tied
    for (vector<int>::iterator ip = offpart.begin(); ip!=p; ++ip) {

        int offset = NXMAPR[*ip];
        for (int j = 0; j < NXMAPR[*ip+1]-NXMAPR[*ip]; ++j)
            nodesp.push_back(XMAPR[NXMAPR[*ip]+j]);
    }

    // -- sort all nodes
    sort(nodesp.begin(),nodesp.end());

    // -- assert that both are numnpn
    assert(numnpn==nodesp.size());

    int nt= 0;
    for (vector<int>::iterator ip = nodesp.begin(); ip!=nodesp.end(); ++ip)
        nodes[nt++] = *ip;
}


void  ME::get_tie_map(int n, int numnpn, int* nodes, int* xmap)
{
    for (int i = 0; i < numnpn; ++i)
        xmap[i] = XMAP[nodes[i]];
}


void  ME::get_node_partition(int n, int numnpn, int* nodes, int* parts)
{
    for (int i = 0; i < numnpn; ++i)
        parts[i] = NPART[nodes[i]];
}


void  ME::get_node_partition(int* parts)
{
    int numnp = mesh->numnp();
    for (int i = 0; i < numnp; ++i)
        parts[i] = NPART[i];
}


void  ME::get_node_ids(int n, int numnpn, int* nodes, int* nodeids)
{
    int ndf = mesh->get_ndf();
    for (int i = 0; i < numnpn; ++i)
        for (int j = 0; j < ndf; ++j)
            nodeids[i*ndf+j] = mesh->id(j,nodes[i]);
}


/*@t -----
 * \subsection{Extracting element information for each partition}
 *
 *
 *@c*/


int ME::get_numelt(int n)
{
    int numelt = mesh->numelt();
    int numeltn= 0;

    for (int i = 0; i < numelt; ++i)
        for (int j = NEPART[i]; j < NEPART[i+1]; ++j)
            if (EPART[j]==n)
                numeltn++;

    return numeltn;
}


void ME::get_elements(int n, int numeltn, int* elts)
{
    int numelt = mesh->numelt();
    int nt     = 0;

    for (int i = 0; i < numelt; ++i)
        for (int j = NEPART[i]; j < NEPART[i+1]; ++j)
            if (EPART[j]==n)
                elts[nt++] = i;

    // -- assert that both are numeltn
    assert(nt==numeltn);
}


void ME::get_element_partition(int n, int numeltn, int* elts, int* parts)
{
  for (int i = 0; i < numeltn; ++i)
      parts[i] = EPART[NEPART[elts[i]]];
}
