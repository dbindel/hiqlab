#include <vector>
#include <iostream>
#include <cassert>

#include "mesh_manager.h"
#include "cscindexer.h"
#include "cscassembly.h"

using std::vector;
#define ME Mesh_Manager


/*@t -----
 * \section{Mesh_Manager}
 *
 * \subsection{Mesh_Manager construction, destruction, and memory management}
 *
 *@c*/


ME::ME()
{
}


ME::~ME()
{
}


/*@t -----
 * \subsection{Adding meshes and computing total numids}
 *
 *@c*/


int ME::add_mesh(Mesh_Partition* mesh)
{
    meshes.push_back(mesh);
}


int ME::compute_numid()
{
    int tnumid = 0;
    int nmesh = meshes.size();

    for (int i = 0; i < nmesh; ++i)
        tnumid += meshes[i]->get_pnumid();

    return tnumid;
}


/*@t -----
 * \subsection{Partitioning the mesh}
 *
 *@c*/


void ME::set_partition(Mesh_Partitioner* meshp)
{
    int nmesh = nummesh();
    assert(nmesh==meshp->get_num_partitions());

    for (int i = 0; i < nmesh; ++i) {

        int pid    = meshes[i]->get_pid();
        int numnpn = meshp->get_numnp (pid);
        int numeltn= meshp->get_numelt(pid);
        vector<int> nodes(numnpn);
        vector<int> xmap(numnpn);
        vector<int> nparts(numnpn);
        vector<int> elts(numeltn);
        vector<int> eparts(numeltn);

        meshp->get_nodes         (pid,numnpn,&nodes[0]);
        meshp->get_tie_map       (pid,numnpn,&nodes[0],&xmap[0]  );
        meshp->get_node_partition(pid,numnpn,&nodes[0],&nparts[0]);
        meshp->get_elements         (pid,numeltn,&elts[0]);
        meshp->get_element_partition(pid,numeltn,&elts[0],&eparts[0]);

        meshes[i]->set_nodes    (numnpn,&nodes[0]);
        meshes[i]->set_ntie_map (numnpn,&xmap[0]);
        meshes[i]->set_npart_map(numnpn,&nparts[0]);
        meshes[i]->set_elements (numeltn,&elts[0]);
        meshes[i]->set_epart_map(numeltn,&eparts[0]);

    }
}


/*@t -----
 * \subsection{Initializing the mesh}
 *
 *@c*/


void ME::initialize()
{
    int nmesh = nummesh();

    for (int i = 0; i < nmesh; ++i)
        meshes[i]->initialize();
} 


/*@t -----
 * \subsection{Assigning reduced ids}
 *
 *@c*/


void ME::set_ids(Mesh_Partitioner* meshp)
{
    int nmesh = nummesh();
    assert(nmesh==meshp->get_num_partitions());

    for (int i = 0; i < nmesh; ++i) {

        int pid    = meshes[i]->get_pid();
        int numnpn = meshp->get_numnp (pid);
        int ndf    = meshp->get_ndf();
        vector<int> nodes(numnpn);
        vector<int> nodeids(numnpn*ndf);

        meshp->get_nodes   (pid,numnpn,&nodes[0]);
        meshp->get_node_ids(pid,numnpn,&nodes[0],&nodeids[0]);

        meshes[i]->set_ids_range(meshp->get_startid(pid),meshp->get_end1id(pid));
        meshes[i]->set_node_ids(numnpn,ndf,&nodeids[0]);

    }
}


/*@t -----
 * \subsection{Assembling the residual}
 *
 *@c*/


void ME::assemble_R()
{
    int nmesh = nummesh();

    for (int i = 0; i < nmesh; ++i)
        meshes[i]->assemble_R();
}


/*@t ----------
 * \section{Building CSCMatrices with Mesh_Manager}
 *
 * \subsection{Functions for building sparsity structure}
 *@c*/


void csc_sparsity_count(Mesh_Manager* mm, int n, int* counts, int reduced)
{
    int nmesh = mm->nummesh();

    CSCSparsityCounter counter(n);
    for (int i = 0; i < nmesh; ++i)
        mm->get_mesh(i)->assemble_struct(&counter, reduced);
    counter.get_counts(counts);
    int nnz = 0;
    for (int j = 0; j < n; ++j) {
        int nnzj = counts[j];
        counts[j] = nnz;
        nnz += nnzj;
    }
    counts[n] = nnz;
}


void csc_build_index(Mesh_Manager* mm, int n, int* jc, int* ir, int reduced)
{
    int nmesh = mm->nummesh();

    CSCIndexBuilder indexer(n, jc, ir);
    for (int i = 0; i < nmesh; ++i)
        mm->get_mesh(i)->assemble_struct(&indexer, reduced);
}


void build_csc_matrix(Mesh_Manager* mm, int n,
                      vector<int>& jc, vector<int>& ir,
                      int reduced)
{
    jc.resize(n+1);
    csc_sparsity_count(mm, n, &jc[0], reduced);
    int nnz = jc[n];
    ir.resize(nnz);
    csc_build_index(mm, n, &jc[0], &ir[0], reduced);
}


CSCMatrix* build_csc_matrix(Mesh_Manager* mm, int reduced)
{
    int n = mm->compute_numid();
    CSCMatrix* result = new CSCMatrix(n, n, 0);
    int* jc = result->get_jc();
    csc_sparsity_count(mm, n, jc, reduced);
    int nnz = jc[n];
    result->set_nnz(nnz);
    int* ir = result->get_ir();
    csc_build_index(mm, n, jc, ir, reduced);
    return result;
}
#undef ME


/*@t ----------
 * \subsection{Assembling the stiffness and mass matrices}
 *
 * Matrices are constructed in a two pass algorithm. First the
 * non-zero structure is computed. Then the matrix is filled.
 *@c*/


CSCMatrix* assemble_dR(Mesh_Manager* mm, int cx, int cv, int ca, int reduced)
{
    int nmesh = mm->nummesh();

    // -- form nonzero structure
    CSCMatrix* K = build_csc_matrix(mm,reduced);

    // -- insert values into structure
    CSCAssembler assembler(K);
    for (int i = 0; i < nmesh; ++i)
        mm->get_mesh(i)->assemble_dR(&assembler, cx, cv, ca, reduced);

    return K;
}


/*@t ----------
 * \subsection{Obtaining the residual}
 *
 * To obtain the residual, assemble_R must be called on Mesh_Manager.
 * Then the information is extracted.
 *@c*/


void get_f(Mesh_Manager* mm, double* vr, double* vi)
{
    int nmesh = mm->nummesh();
    Mesh_Partition* mnp;

    for (int i = 0; i < nmesh; ++i) {

        mnp = mm->get_mesh(i);
        int M = mnp->get_idsize();

        for (int j = 0; j < M; ++j) {
            int ind = mnp->idg(j);
            if (ind >= mnp->get_pstartid() && ind < mnp->get_pend1id()) {
                vr[ind] = mnp->f(j);
                if (vi)
                    vi[ind] = mnp->fi(j);
            }
        }
    }
}


/*@t ----------
 * \section{Building CSCMatrix Prolongators with QAddBlockProlongator}
 *
 *@c*/


void csc_sparsity_count(QAddBlockProlongator* prolongator, int m, int n, int* counts, int reduced)
{
    CSCSparsityCounter counter(m,n);
    prolongator->assemble_struct(&counter, reduced);
    counter.get_counts(counts);
    int nnz = 0;
    for (int j = 0; j < n; ++j) {
        int nnzj = counts[j];
        counts[j] = nnz;
        nnz += nnzj;
    }
    counts[n] = nnz;
}


void csc_build_index(QAddBlockProlongator* prolongator, int m, int n, int* jc, int* ir, int reduced)
{
    CSCIndexBuilder indexer(m, n, jc, ir);
    prolongator->assemble_struct(&indexer, reduced);
}


CSCMatrix* build_csc_matrix(QAddBlockProlongator* prolongator, int m, int n, int reduced)
{
    CSCMatrix* result = new CSCMatrix(m, n, 0);
    int* jc = result->get_jc();
    csc_sparsity_count(prolongator, m, n, jc, reduced);

    int nnz = jc[n];
    result->set_nnz(nnz);
    int* ir = result->get_ir();
    csc_build_index(prolongator, m, n, jc, ir, reduced);
    return result;
}


CSCMatrix* assemble_P(Mesh* from, Mesh* to, int* nemap, int reduced)
{
    // -- Build QAddBlockProlongator structure
    QAddBlockProlongator prolongator(from,to,nemap);

    // -- form nonzero structure
    assert(reduced==1); //FIXME: So that reduced and not, as well as subassembly work
    int m = to->get_numid();
    int n = from->get_numid();
    CSCMatrix* P = build_csc_matrix(&prolongator,m,n,reduced);

    // -- insert values into structure
    CSCAssembler assembler(P);
    prolongator.assemble_P(&assembler,reduced); 

    return P;
}


/*@t ----------
 * \section{Assembling the Prolongator with Mesh_Manager}
 *
 *@c*/


#define ME Prolongator_Mesh_Manager


ME::ME(Mesh_Manager* from, Mesh_Manager* to) : from(from), to(to)
{
    assert(from->nummesh()==to->nummesh());
}


ME::~ME()
{
}


void ME::set_nemap(int numnp, int * nmp)
{
    NEMAP.resize(numnp);
    for (int i = 0; i < numnp; ++i)
        NEMAP[i] = nmp[i];
}


void ME::construct_partition_nemap(int n, int numnpn, int* nemap_global)
{
    Mesh_Partition* mpt = to->get_mesh(n);

    for (int i = 0; i < numnpn; ++i)
        nemap_global[i] = NEMAP[mpt->nlgmap(i)];
}


/*@t ----------
 * \section{Building CSCMatrices with Prolongator_Mesh_Manager}
 *
 * \subsection{Functions for building sparsity structure}
 *@c*/


void csc_sparsity_count(Prolongator_Mesh_Manager* pmm, int m, int n, int* counts, int reduced)
{
    int nmesh = pmm->get_to()->nummesh();

    CSCSparsityCounter counter(m,n);
    for (int i = 0; i < nmesh; ++i) {

        // -- Construct partition nemap
        int numnp_to = pmm->get_to()->get_mesh(i)->numnp();
        vector<int> nemap_global(numnp_to);
        pmm->construct_partition_nemap(i, numnp_to, &nemap_global[0]);
        QAddBlockProlongator_Partition abp(pmm->get_from()->get_mesh(i),
                                           pmm->get_to()->get_mesh(i),&nemap_global[0]);

        abp.assemble_struct(&counter, reduced);

    }
    counter.get_counts(counts);
    int nnz = 0;
    for (int j = 0; j < n; ++j) {
        int nnzj = counts[j];
        counts[j] = nnz;
        nnz += nnzj;
    }
    counts[n] = nnz;
}


void csc_build_index(Prolongator_Mesh_Manager* pmm, int m, int n, int* jc, int* ir, int reduced)
{
    int nmesh = pmm->get_to()->nummesh();

    CSCIndexBuilder indexer(m, n, jc, ir);
    for (int i = 0; i < nmesh; ++i) {

        // -- Construct partition nemap
        int numnp_to = pmm->get_to()->get_mesh(i)->numnp();
        vector<int> nemap_global(numnp_to);
        pmm->construct_partition_nemap(i, numnp_to, &nemap_global[0]);
        QAddBlockProlongator_Partition abp(pmm->get_from()->get_mesh(i),
                                           pmm->get_to()->get_mesh(i),&nemap_global[0]);

        abp.assemble_struct(&indexer, reduced);

    }
}


CSCMatrix* build_csc_matrix(Prolongator_Mesh_Manager* pmm, int reduced)
{
    assert(reduced==1);

    int m = pmm->get_to()->compute_numid();
    int n = pmm->get_from()->compute_numid();
    CSCMatrix* result = new CSCMatrix(m, n, 0);
    int* jc = result->get_jc();
    csc_sparsity_count(pmm, m, n, jc, reduced);
    int nnz = jc[n];
    result->set_nnz(nnz);
    int* ir = result->get_ir();
    csc_build_index(pmm, m, n, jc, ir, reduced);
    return result;
}


/*@t ----------
 * \subsection{Assembling the prolongator matrix}
 *
 * Matrices are constructed in a two pass algorithm. First the
 * non-zero structure is computed. Then the matrix is filled.
 *@c*/


CSCMatrix* assemble_P(Prolongator_Mesh_Manager* pmm, int reduced)
{
    int nmesh = pmm->get_to()->nummesh();

    // -- form nonzero structure
    CSCMatrix* P = build_csc_matrix(pmm,reduced);

    // -- insert values into structure
    CSCAssembler assembler(P);
    for (int i = 0; i < nmesh; ++i) {

        // -- Construct partition nemap
        int numnp_to = pmm->get_to()->get_mesh(i)->numnp();
        vector<int> nemap_global(numnp_to);
        pmm->construct_partition_nemap(i, numnp_to, &nemap_global[0]);
        QAddBlockProlongator_Partition abp(pmm->get_from()->get_mesh(i),
                                           pmm->get_to()->get_mesh(i),&nemap_global[0]);

        abp.assemble_P(&assembler, reduced);
    }

    return P;
}

