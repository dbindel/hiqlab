/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <vector>
#include <iostream>
#include <cassert>

#include "mesh_partitioner_metis.h"
#include "cscindexer.h"
#include "metispart.h"
#include "coordsorter.h"

using std::vector;

#define ME Mesh_Partitioner_METIS


/*@t --------
 * \subsection{Mesh_Partitioner_METIS construction, destruction, 
 *             and memory management}
 *
 *@c*/


ME::ME(int nparts, Mesh* mesh, int ptype) :
    Mesh_Partitioner(nparts), mesh(mesh), metis_type(ptype), numrn(0)
{
}


ME::~ME()
{
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
 * The NPART vector is initialized so that it contains the 
 * reduced ids that are given to METIS to partition. 
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
    // NPART[i] := 1 if it is class_id
    XMAP.resize(np);
    for (int i = 0; i < np; ++i) 
        XMAP[i] = i;

    NPART.clear();
    NPART.resize(np);
    numrn = 0;
    for (int inext, i = 0; i < nt; i = inext) {
        int class_id = scratch[i];
        for (inext = i+1; inext < nt; ++inext) {
             if (sorter.compare(scratch[i], scratch[inext]) != 0)
                 break;
             if (scratch[inext] < scratch[i])
                class_id = scratch[inext];
        }
        for (int j = i; j < inext; ++j)
            XMAP[scratch[j]] = class_id;
        numrn++;
        NPART[class_id] = 1;
    }

    // Construct temporary relabeling of the nodes for partitioning
    // NPART[i] := -1 if not class_id otherwise number
    int inum = 0;
    for (int i = 0; i < np; ++i)
        NPART[i] = ( (NPART[i] == 1) ? inum++ : -1 );

    assert(numrn==inum);

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
 * /subsection{Partitioning nodes with METIS}
 * 
 * A graph of the nodal connectivity must be formed
 * to allow METIS to construct a nonoverlapping nodal
 * partition.
 *
 *@c*/


void ME::build_adj(std::vector<int>& jc, std::vector<int>& ir)
{
    jc.resize(numrn+1);
    adj_sparsity_count(&jc[0]);
    int nnz = jc[numrn];
    ir.resize(nnz);
    adj_build_index(&jc[0], &ir[0]);
}


void ME::adj_sparsity_count(int* counts)
{
    CSCSparsityCounter counter(numrn);

    int numelt = mesh->numelt();
    for (int ie = 0; ie < numelt; ++ie) {

        int nen = mesh->get_nen(ie);

        // -- Add upper half
        for (int j = 0; j < nen; ++j)
            for (int i = j+1; i < nen; ++i) {
                int iie = NPART[mesh->ix(i,ie)];
                int jie = NPART[mesh->ix(j,ie)];
                counter.add(iie, jie);
                counter.add(jie, iie);
            }
    }

    counter.get_counts(counts);
    int nnz = 0;
    for (int j = 0; j < numrn; ++j) {
        int nnzj = counts[j];
        counts[j] = nnz;
        nnz += nnzj;
    }
    counts[numrn] = nnz;
}


void ME::adj_build_index(int* jc, int* ir)
{
    CSCIndexBuilder indexer(numrn, jc, ir);

    int numelt = mesh->numelt();
    for (int ie = 0; ie < numelt; ++ie) {

        int nen = mesh->get_nen(ie);

        // -- Add upper half
        for (int j = 0; j < nen; ++j)
            for (int i = j+1; i < nen; ++i) {
                int iie = NPART[mesh->ix(i,ie)];
                int jie = NPART[mesh->ix(j,ie)];
                indexer.add(iie, jie);
                indexer.add(jie, iie);
            }
    }
}


void ME::partition_nodes()
{
    vector<int> jc;
    vector<int> ir;
    vector<int> part;
    int         edgecut;

    // -- build adjacency structure
    build_adj(jc,ir);
    part.resize(numrn);

    // -- partition with METIS
    Metis partitioner(numrn,&jc[0],&ir[0]);
    if (metis_type==0)
        partitioner.PartGraphRecursive(num_partitions,&part[0],&edgecut);
    else if(metis_type==1)
        partitioner.PartGraphKway(num_partitions,&part[0],&edgecut);

    // -- copy data into NPART
    for (int i = 0; i < mesh->numnp(); ++i)
        if (NPART[i]>=0)
            NPART[i] = part[NPART[i]];
    for (int i = 0; i < mesh->numnp(); ++i)
        if (NPART[i] < 0)
            NPART[i] = NPART[XMAP[i]];

}

/*@t -----
 * /subsection{Constructing an overlapping element partition}
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
        find_partitions(&(mesh->ix(0,i)), nen, part, &npart);

        NEPART.push_back(npart);
        for (int j = 0; j < npart; ++j)
            EPART.push_back(part[j]);
    }

    for (int i = 1; i < numelt+1; ++i)
        NEPART[i] += NEPART[i-1];

}

void ME::find_partitions(int* nodes, int n, std::vector<int>& part, int* npart)
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
    std::vector<int> offpart;

    for (int i = 0; i < numelt; ++i)
        for (int j = NEPART[i]; j < NEPART[i+1]; ++j)
            if (EPART[j]==n) {
                int nen = mesh->get_nen(i);
                for (int k = 0; k < nen; ++k)
                    if (NPART[mesh->ix(k,i)] != n)
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
