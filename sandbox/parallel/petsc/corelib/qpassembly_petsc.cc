/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include "mesh.h"
#include "mesh_add_block.h"
#include "qpassembly_petsc.h"
#include "qcomplex.h"
#include "csrindexer.h"

using std::vector;

/*@t ----------
 * \section{QPETScStructAssembler methods}
 *
 *@c*/
#define ME QPETScStructAssembler


ME::~ME()
{
}


void ME::add(int i, int j)
{
    if (is_reduced==QP_ASSEMBLE_REMOVE_DIRICHLET) {

        assembler.add(mesh->id(i)-Istart,mesh->id(j));

    } else if (is_reduced==QP_ASSEMBLE_NO_DIRICHLET) {

        assembler.add(i-Istart,j);

    } else if (is_reduced==QP_ASSEMBLE_IDENTITY_DIRICHLET) {

        assembler.add(i-Istart,j);

    }
}


void ME::add(int* eltid, int n)
{
    ME::add(eltid,n,eltid,n);
}


void ME::add(int* eltidm, int m, int* eltidn, int n)
{
    vector<int> mapped_eltidm(m);
    vector<int> mapped_eltidn(n);

    map_eltid(eltidm,m,eltidn,n,mapped_eltidm,mapped_eltidn);
    filter_eltid(&(mapped_eltidm[0]),m,Istart,Iend1,-Istart);

    assembler.add(&(mapped_eltidm[0]), m, &(mapped_eltidn[0]), n);
}


void ME::map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn)
{
    if (is_reduced==QP_ASSEMBLE_REMOVE_DIRICHLET) {
 
        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = mesh->id(eltidm[i]); 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = mesh->id(eltidn[i]);


    } else if (is_reduced==QP_ASSEMBLE_NO_DIRICHLET) {

        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = eltidm[i]; 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = eltidn[i];


    } else if (is_reduced==QP_ASSEMBLE_IDENTITY_DIRICHLET) {

        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = eltidm[i]; 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = eltidn[i];
    }
}


void ME::filter_eltid(int* eltid, int n, int sid, int eid1, int offset)
{
    for (int i = 0; i < n; ++i) 
        if (eltid[i] < sid || eltid[i] >=eid1)
            eltid[i] = -3;
        else
            eltid[i] += offset;
}


#undef ME


/*@t ----------
 * \section{QPETScAssembler methods}
 *
 *@c*/
#define ME QPETScAssembler


ME::~ME()
{
}


void ME::add(int* eltid, int n, dcomplex* Ke)
{
    ME::add(eltid, n, eltid, n, Ke);
}


void ME::add(int* eltid, int n, double* Ke)
{
    ME::add(eltid, n, eltid, n, Ke);
}


void ME::add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke)
{
#ifdef PETSC_USE_COMPLEX
    if (!is_assemble_real) {
        dcomplex ione(0,1);
        for (int i = 0; i < m*n; ++i)
            Ke[i] = Ke[i]*ione;
    }
    add<dcomplex>(eltidm, m, eltidn, n, Ke);
#endif
}


void ME::add(int* eltidm, int m, int* eltidn, int n, double* Ke)
{
#ifndef PETSC_USE_COMPLEX
    add<double>(eltidm, m, eltidn, n, Ke);
#endif
}


template<class T>
void ME::add(int* eltidm, int m, int* eltidn, int n, T* Ke)
{
    vector<int> mapped_eltidm(m);
    vector<int> mapped_eltidn(n);

    map_eltid<T>(eltidm,m,eltidn,n,mapped_eltidm,mapped_eltidn, Ke);
    filter_eltid(&(mapped_eltidm[0]),m,Istart,Iend1,0);

    MatSetValues(A, m, &(mapped_eltidm[0]), n, &(mapped_eltidn[0]),Ke,ADD_VALUES);
}


template<class T>
void ME::map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn, T* Ke)
{
    if (is_reduced==QP_ASSEMBLE_REMOVE_DIRICHLET) {
 
        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = mesh->id(eltidm[i]); 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = mesh->id(eltidn[i]);


    } else if (is_reduced==QP_ASSEMBLE_NO_DIRICHLET) {

        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = eltidm[i]; 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = eltidn[i];


    } else if (is_reduced==QP_ASSEMBLE_IDENTITY_DIRICHLET) {

        apply_identity_dirichlet(eltidm, m, Ke);
        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = eltidm[i]; 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = eltidn[i];
    }
}


void ME::filter_eltid(int* eltid, int n, int sid, int eid1, int offset)
{
    for (int i = 0; i < n; ++i) 
        if (eltid[i] < sid || eltid[i] >=eid1)
            eltid[i] = -3;
        else
            eltid[i] += offset;
}


template<class T>
void ME::apply_identity_dirichlet(int* eltid, int n, T* Ke)
{
    for (int i = 0; i < n; ++i) {
        if (mesh->id(eltid[i]) < 0) {
            for (int j = 0; j < n; ++j) {
#ifdef PETSC_USE_COMPLEX
                Ke[j*n+i] = dcomplex(0,0);
                Ke[i*n+j] = dcomplex(0,0);
#else
                Ke[j*n+i] = 0;
                Ke[i*n+j] = 0;
#endif
            }
        }
    }
}


template<class T>
void ME::add_identity_dirichlet(T diag)
{

    PetscInt Istart, Iend1;
    MatGetOwnershipRange(A, &Istart, &Iend1);
    for (int i = Istart; i < Iend1; ++i)
        if (mesh->id(i) < 0)
            MatSetValue(A, i, i, diag, ADD_VALUES);

}


#undef ME

/*@t ----------
 * \section{QPETScProlongatorStructAssembler methods}
 *
 *@c*/
#define ME QPETScProlongatorStructAssembler


ME::~ME()
{
}


void ME::add(int i, int j)
{
    if (is_reduced==QP_ASSEMBLE_REMOVE_DIRICHLET) {

        assembler.add(prolongator->get_to()->id(i)-Istart,prolongator->get_from()->id(j));

    } else if (is_reduced==QP_ASSEMBLE_NO_DIRICHLET) {

        assembler.add(i-Istart,j);

    } else if (is_reduced==QP_ASSEMBLE_IDENTITY_DIRICHLET) {

        assembler.add(i-Istart,j);

    }
}


void ME::add(int* eltid, int n)
{
    ME::add(eltid,n,eltid,n);
}


void ME::add(int* eltidm, int m, int* eltidn, int n)
{
    vector<int> mapped_eltidm(m);
    vector<int> mapped_eltidn(n);

    map_eltid(eltidm,m,eltidn,n,mapped_eltidm,mapped_eltidn);
    filter_eltid(&(mapped_eltidm[0]),m,Istart,Iend1,-Istart);

    assembler.add(&(mapped_eltidm[0]), m, &(mapped_eltidn[0]), n);
}


void ME::map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn)
{
    if (is_reduced==QP_ASSEMBLE_REMOVE_DIRICHLET) {
 
        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = prolongator->get_to()->id(eltidm[i]); 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = prolongator->get_from()->id(eltidn[i]);


    } else if (is_reduced==QP_ASSEMBLE_NO_DIRICHLET) {

        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = eltidm[i]; 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = eltidn[i];


    } else if (is_reduced==QP_ASSEMBLE_IDENTITY_DIRICHLET) {

        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = eltidm[i]; 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = eltidn[i];
    }
}


void ME::filter_eltid(int* eltid, int n, int sid, int eid1, int offset)
{
    for (int i = 0; i < n; ++i) 
        if (eltid[i] < sid || eltid[i] >=eid1)
            eltid[i] = -3;
        else
            eltid[i] += offset;
}


#undef ME


/*@t ----------
 * \section{QPETScProlongatorAssembler methods}
 *
 *@c*/
#define ME QPETScProlongatorAssembler


ME::~ME()
{
}


void ME::add(int* eltid, int n, dcomplex* Ke)
{
    ME::add(eltid, n, eltid, n, Ke);
}


void ME::add(int* eltid, int n, double* Ke)
{
    ME::add(eltid, n, eltid, n, Ke);
}


void ME::add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke)
{
    add<dcomplex>(eltidm, m, eltidn, n, Ke);
}


void ME::add(int* eltidm, int m, int* eltidn, int n, double* Ke)
{
#ifdef PETSC_USE_COMPLEX
    vector<dcomplex> Kez(m*n);
    for (int i = 0; i < m*n; ++i)
        Kez[i] = dcomplex(Ke[i],0);
    add<dcomplex>(eltidm, m, eltidn, n, &(Kez[0]));
#else
    add<double>(eltidm, m, eltidn, n, Ke);
#endif
}


template<class T>
void ME::add(int* eltidm, int m, int* eltidn, int n, T* Ke)
{
    vector<int> mapped_eltidm(m);
    vector<int> mapped_eltidn(n);

    map_eltid<T>(eltidm,m,eltidn,n,mapped_eltidm,mapped_eltidn, Ke);
    filter_eltid(&(mapped_eltidm[0]),m,Istart,Iend1,0);

    MatSetValues(A, m, &(mapped_eltidm[0]), n, &(mapped_eltidn[0]),Ke,ADD_VALUES);
}


template<class T>
void ME::map_eltid(int* eltidm, int m, int* eltidn, int n, vector<int>& mapped_eltidm, vector<int>& mapped_eltidn, T* Ke)
{
    if (is_reduced==QP_ASSEMBLE_REMOVE_DIRICHLET) {
 
        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = prolongator->get_to()->id(eltidm[i]); 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = prolongator->get_from()->id(eltidn[i]);


    } else if (is_reduced==QP_ASSEMBLE_NO_DIRICHLET) {

        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = eltidm[i]; 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = eltidn[i];


    } else if (is_reduced==QP_ASSEMBLE_IDENTITY_DIRICHLET) {

        for (int i = 0; i < m; ++i)
            mapped_eltidm[i] = eltidm[i]; 

        for (int i = 0; i < n; ++i)
            mapped_eltidn[i] = eltidn[i];
    }
}


void ME::filter_eltid(int* eltid, int n, int sid, int eid1, int offset)
{
    for (int i = 0; i < n; ++i) 
        if (eltid[i] < sid || eltid[i] >=eid1)
            eltid[i] = -3;
        else
            eltid[i] += offset;
}


/*@t ----------
 * \section{Assembling K with Petsc Mat methods}
 *
 *@c*/


static int num_loc_rows(int n, int ndf, MPI_Comm comm)
{
    int mpirank;
    int mpisize;
    int numnp = n/ndf;

    MPI_Comm_rank(comm, &mpirank);
    MPI_Comm_size(comm, &mpisize);

    int avep = (numnp/mpisize) * ndf;
    int lastp= (numnp%mpisize) * ndf + avep;

    if (mpirank < mpisize-1)
       return avep;
    else
       return lastp;
}


static void csr_sparsity_count(Mesh* mesh, int is_reduced, int num_global_rows, PetscInt Istart, PetscInt Iend1, int* counts)
{
    int num_my_rows = Iend1 - Istart;
    CSRSparsityCounter counter(Iend1-Istart, num_global_rows);
    QPETScStructAssembler assembler(mesh,Istart,Iend1,counter,is_reduced);
    mesh->assemble_struct_raw(&assembler);

    counter.get_counts(counts);
    int nnz = 0;
    for (int j = 0; j < num_my_rows; ++j) {
        int nnzj = counts[j];
        counts[j] = nnz;
        nnz += nnzj;
    }
    counts[num_my_rows] = nnz;
}


static void csr_build_index(Mesh* mesh, int is_reduced, int num_global_rows, PetscInt Istart, PetscInt Iend1, int* ir, int* jc)
{
    CSRIndexBuilder indexer(Iend1-Istart, num_global_rows, &ir[0], &jc[0]); 
    QPETScStructAssembler assembler(mesh,Istart,Iend1,indexer,is_reduced);
    mesh->assemble_struct_raw(&assembler);
}


static void build_nnz_structure(Mesh* mesh, int is_reduced, int num_global_rows,
                                PetscInt Istart, PetscInt Iend1,
                                int* d_nnz, int* o_nnz, int blocked)
{
    int num_my_rows = Iend1 - Istart;
    vector<int> ir;
    vector<int> jc;

    ir.resize(Iend1-Istart+1);

    // -- csr_sparsity_count
    csr_sparsity_count(mesh,is_reduced,num_global_rows,Istart,Iend1,&(ir[0]));

    // -- intermediate
    int nnz = ir[num_my_rows];
    jc.resize(nnz);

    // -- csr_build_index
    csr_build_index(mesh,is_reduced,num_global_rows,Istart,Iend1,&(ir[0]),&(jc[0]));


    // -- build d_nnz and o_nnz
    //    if blocked==1, then take every [ndf] rows of [ir] to construct 
    //    [d_nnz] and [o_nnz]
    int ndf = (blocked ? mesh->get_ndf() : 1);
    for (int i = 0; i < num_my_rows/ndf; ++i) {
        int rnnz = 0;
        for (int j = ir[i*ndf]; j < ir[i*ndf+1]; ++j)
            if (Istart <= jc[j] && jc[j] < Iend1)
                rnnz++;
        d_nnz[i] = rnnz/ndf;
        o_nnz[i] = ((ir[i*ndf+1]-ir[i*ndf]) - rnnz)/ndf;
    }

} 


static void mat_preallocate_memory(Mesh* mesh, Mat *A, int num_global_rows, int mat_type, int is_reduced)
{
    int blocked = 0;
    if (mat_type==0) {
        MatSetType(*A, MATSEQAIJ);
    } else if (mat_type==2) {
        MatSetType(*A, MATMPIAIJ);
    } else if (mat_type==1) {
        blocked = 1;
        MatSetType(*A, MATSEQBAIJ);
    } else if (mat_type==3) {
        blocked = 1;
        MatSetType(*A, MATMPIBAIJ);
    }
    int  ndf   = (blocked ? mesh->get_ndf() : 1);

    PetscInt Istart, Iend1;
    MatGetOwnershipRange(*A,&Istart,&Iend1);
    int num_my_rows = (Iend1 - Istart);
    vector<int> d_nnz((Iend1-Istart)/ndf);
    vector<int> o_nnz((Iend1-Istart)/ndf);

    build_nnz_structure(mesh, is_reduced, num_global_rows, Istart, Iend1, &(d_nnz[0]), &(o_nnz[0]), blocked);

    if (mat_type==0) {
//      MatSeqAIJSetPreallocation(A, 0, &(d_nnz[0]));
        MatDestroy(*A);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, num_my_rows, num_my_rows, 
                        0, &(d_nnz[0]), A);

    } else if (mat_type==2) {
//      MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz);
        MatDestroy(*A);
        MatCreateMPIAIJ(PETSC_COMM_WORLD, num_my_rows, num_my_rows, 
                   PETSC_DETERMINE, PETSC_DETERMINE, 0, &(d_nnz[0]), 0, &(o_nnz[0]), A);

    } else if (mat_type==1) {
//      MatSeqBAIJSetPreallocation(*A, ndf, 0, d_nnz);   //FIXME: It seems you can do both.
        MatDestroy(*A);
        MatCreateSeqBAIJ(PETSC_COMM_WORLD, ndf, num_my_rows, num_my_rows, 
                        0, &(d_nnz[0]), A);

    } else if (mat_type==3) {
//      MatMPIBAIJSetPreallocation(*A, ndf, 0, d_nnz, 0, o_nnz);
        MatDestroy(*A);
        MatCreateMPIBAIJ(PETSC_COMM_WORLD, ndf, num_my_rows, num_my_rows, 
                   PETSC_DETERMINE, PETSC_DETERMINE, 0, &(d_nnz[0]), 0, &(o_nnz[0]), A);
    }

}


static void get_mat_sizes(Mesh* mesh, int* n, int* nl, int is_reduced)
{
    if (is_reduced == QP_ASSEMBLE_NO_DIRICHLET) {

        *n  = mesh->numnp() * mesh->get_ndf();
        *nl = num_loc_rows(*n, mesh->get_ndf(),  PETSC_COMM_WORLD); 

    } else if (is_reduced == QP_ASSEMBLE_REMOVE_DIRICHLET) {

        *n  = mesh->get_numid();
        *nl = num_loc_rows(*n, 1, PETSC_COMM_WORLD); 

    } else if (is_reduced == QP_ASSEMBLE_IDENTITY_DIRICHLET) {

        *n  = mesh->numnp() * mesh->get_ndf();
        *nl = num_loc_rows(*n, mesh->get_ndf(),  PETSC_COMM_WORLD); 

    }
}


Mat Mesh_assemble_dR_petsc(Mesh* mesh, double cxr, double cxi, 
                                       double cvr, double cvi,
                                       double car, double cai,
                           int mat_type_code, int is_reduced)
{
    Mat A;
    int n, num_my_rows;
    PetscInt Istart, Iend1;

    // -- Determine size of global and local degrees of freedom
    get_mat_sizes(mesh,&n,&num_my_rows,is_reduced);

    // -- Create matrix
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, num_my_rows, num_my_rows, n, n);

    // -- Preallocate memory
    mat_preallocate_memory(mesh,&A,n,mat_type_code,is_reduced);

    MatSetFromOptions(A);

    // -- Assemble real part
    MatGetOwnershipRange(A,&Istart,&Iend1);
    QPETScAssembler assembler(mesh,A,Istart,Iend1,is_reduced);
    assembler.assemble_real();
    mesh->assemble_dR_raw(&assembler,cxr,cvr,car);

    // -- Assemble imaginary part
#ifdef PETSC_USE_COMPLEX
    if (cxi!=0 || cvi!=0 || cai!=0) {
        assembler.assemble_imag();
        mesh->assemble_dR_raw(&assembler, cxi, cvi, cai);
    }
#endif

    // -- Make diagonal component [diag] for Dirichlet BC nodes
    if (is_reduced == QP_ASSEMBLE_IDENTITY_DIRICHLET) {
#ifdef PETSC_USE_COMPLEX
        dcomplex diag(1.0,0);
#else
        double diag = 1.0;
#endif
        assembler.add_identity_dirichlet(diag);
    }

    // -- Finalize assembly
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    return A;
}

Mat Mesh_assemble_dR_petsc(Mesh* mesh, double cxr, double cvr, double car,
                                       int mat_type_code, int is_reduced)
{
    return Mesh_assemble_dR_petsc(mesh,cxr, 0, cvr, 0, car, 0, mat_type_code, is_reduced);
}


/*@t ----------
 * \section{Assembling R with Petsc Vec methods}
 *
 *@c*/


inline static void vec_set_value(Mesh* mesh, Vec x, int i, double* fr, double* fi, int is_reduced)
{
    int idt;

    if (is_reduced == QP_ASSEMBLE_REMOVE_DIRICHLET)
        idt = i;
    else if (is_reduced == QP_ASSEMBLE_IDENTITY_DIRICHLET)
        idt = mesh->id(i);

    if (idt < 0) {
#ifdef PETSC_USE_COMPLEX
        VecSetValue(x, i, dcomplex(0,0), INSERT_VALUES);
#else
        VecSetValue(x, i,             0, INSERT_VALUES);
#endif
    } else {
#ifdef PETSC_USE_COMPLEX
        VecSetValue(x, i, dcomplex(fr[idt],fi[idt]), INSERT_VALUES);
#else
        VecSetValue(x, i, fr[idt], INSERT_VALUES);
#endif
    }
}


Vec Mesh_assemble_R_petsc(Mesh* mesh, int is_reduced)
{
    Vec x;
    int n, num_my_rows;
    PetscInt Istart, Iend1;
    double* fr;
    double* fi;

    // -- Determine size of global and local degrees of freedom
    get_mat_sizes(mesh,&n,&num_my_rows,is_reduced);

    // -- Create matrix
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, num_my_rows, n);

    VecSetFromOptions(x);

    // -- Assemble
    VecGetOwnershipRange(x,&Istart,&Iend1);
    mesh->assemble_R();

    // -- Obtain residual(real) from mesh
    fr = new double[mesh->get_numid()];
    memset(fr, 0, mesh->get_numid()*sizeof(double));
    mesh->get_f(fr);

    // -- Obtain residual(imag) from mesh
#ifdef PETSC_USE_COMPLEX
    fi = new double[mesh->get_numid()];
    memset(fi, 0, mesh->get_numid()*sizeof(double));
    mesh->get_fi(fi);
#endif

    for (int i = Istart; i < Iend1; ++i)
        vec_set_value(mesh,x,i,fr,fi,is_reduced);

    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

    // -- clean up
    delete[] fr;
#ifdef PETSC_USE_COMPLEX
    delete[] fi;
#endif

    return x;
}


/*@t ----------
 * \section{Assembling Prolongators methods}
 *
 *@c*/


static void csr_sparsity_count(QAddBlockProlongator* prolongator, int is_reduced, int n_num_global_rows, PetscInt Istart, PetscInt Iend1, int* counts)
{
    int m_num_my_rows = Iend1 - Istart;
    CSRSparsityCounter counter(Iend1-Istart, n_num_global_rows);
    QPETScProlongatorStructAssembler assembler(prolongator,Istart,Iend1,counter,is_reduced);
    prolongator->assemble_struct_raw(&assembler);

    counter.get_counts(counts);
    int nnz = 0;
    for (int j = 0; j < m_num_my_rows; ++j) {
        int nnzj = counts[j];
        counts[j] = nnz;
        nnz += nnzj;
    }
    counts[m_num_my_rows] = nnz;
}


static void csr_build_index(QAddBlockProlongator* prolongator, int is_reduced, int n_num_global_rows, PetscInt Istart, PetscInt Iend1, int* ir, int* jc)
{
    CSRIndexBuilder indexer(Iend1-Istart, n_num_global_rows, &ir[0], &jc[0]); 
    QPETScProlongatorStructAssembler assembler(prolongator,Istart,Iend1,indexer,is_reduced);
    prolongator->assemble_struct_raw(&assembler);
}


static void build_nnz_structure(QAddBlockProlongator* prolongator, int is_reduced, int n_num_global_rows,
                                PetscInt Istart, PetscInt Iend1,
                                int* d_nnz, int* o_nnz, int blocked)
{
    int m_num_my_rows = Iend1 - Istart;
    vector<int> ir;
    vector<int> jc;

    ir.resize(Iend1-Istart+1);

    // -- csr_sparsity_count
    csr_sparsity_count(prolongator,is_reduced,n_num_global_rows,Istart,Iend1,&(ir[0]));

    // -- intermediate
    int nnz = ir[m_num_my_rows];
    jc.resize(nnz);

    // -- csr_build_index
    csr_build_index(prolongator,is_reduced,n_num_global_rows,Istart,Iend1,&(ir[0]),&(jc[0]));


    // -- build d_nnz and o_nnz
    //    if blocked==1, then take every [ndf] rows of [ir] to construct 
    //    [d_nnz] and [o_nnz]
    int ndf = (blocked ? prolongator->get_from()->get_ndf() : 1);
    for (int i = 0; i < m_num_my_rows/ndf; ++i) {
        int rnnz = 0;
        for (int j = ir[i*ndf]; j < ir[i*ndf+1]; ++j)
            if (Istart <= jc[j] && jc[j] < Iend1)
                rnnz++;
        d_nnz[i] = rnnz/ndf;
        o_nnz[i] = ((ir[i*ndf+1]-ir[i*ndf]) - rnnz)/ndf;
    }

} 


static void mat_preallocate_memory(QAddBlockProlongator* prolongator, Mat *P, 
                                  int m_num_global_rows, int n_num_global_rows,
                                  int m_num_my_rows,     int n_num_my_rows,    
                                  int mat_type, int is_reduced)
{
    int blocked = 0;
    if (mat_type==0) {
        MatSetType(*P, MATSEQAIJ);
    } else if (mat_type==2) {
        MatSetType(*P, MATMPIAIJ);
    } else if (mat_type==1) {
        blocked = 1;
        MatSetType(*P, MATSEQBAIJ);
    } else if (mat_type==3) {
        blocked = 1;
        MatSetType(*P, MATMPIBAIJ);
    }
    int  ndf   = (blocked ? prolongator->get_from()->get_ndf() : 1);

    PetscInt Istart, Iend1;
    MatGetOwnershipRange(*P,&Istart,&Iend1);
    int num_my_rows = (Iend1 - Istart);
    vector<int> d_nnz((Iend1-Istart)/ndf);
    vector<int> o_nnz((Iend1-Istart)/ndf);

    build_nnz_structure(prolongator, is_reduced, n_num_global_rows, Istart, Iend1, &(d_nnz[0]), &(o_nnz[0]), blocked);

    if (mat_type==0) {
//      MatSeqAIJSetPreallocation(*P, 0, &(d_nnz[0]));
        MatDestroy(*P);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, m_num_my_rows, n_num_my_rows, 
                        0, &(d_nnz[0]), P);

    } else if (mat_type==2) {
//      MatMPIAIJSetPreallocation(*P, 0, d_nnz, 0, o_nnz);
        MatDestroy(*P);
        MatCreateMPIAIJ(PETSC_COMM_WORLD, m_num_my_rows, n_num_my_rows, 
                   PETSC_DETERMINE, PETSC_DETERMINE, 0, &(d_nnz[0]), 0, &(o_nnz[0]), P);

    } else if (mat_type==1) {
//      MatSeqBAIJSetPreallocation(*P, ndf, 0, d_nnz);   //FIXME: It seems you can do both.
        MatDestroy(*P);
        MatCreateSeqBAIJ(PETSC_COMM_WORLD, ndf, m_num_my_rows, n_num_my_rows, 
                        0, &(d_nnz[0]), P);

    } else if (mat_type==3) {
//      MatMPIBAIJSetPreallocation(*P, ndf, 0, d_nnz, 0, o_nnz);
        MatDestroy(*P);
        MatCreateMPIBAIJ(PETSC_COMM_WORLD, ndf, m_num_my_rows, n_num_my_rows, 
                   PETSC_DETERMINE, PETSC_DETERMINE, 0, &(d_nnz[0]), 0, &(o_nnz[0]), P);
    }

}


Mat Mesh_assemble_P_petsc(Mesh* from, Mesh* to, int* nemap, int mat_type_code, int is_reduced)
{
    Mat P;
    int m, m_num_my_rows;
    int n, n_num_my_rows;
    PetscInt Istart, Iend1;

    // -- Construct prolongator object
    QAddBlockProlongator prolongator(from,to,nemap);

    // -- Determine size of global and local degrees of freedom
    get_mat_sizes(to  ,&m,&m_num_my_rows,is_reduced);
    get_mat_sizes(from,&n,&n_num_my_rows,is_reduced);

    // -- Create matrix
    MatCreate(PETSC_COMM_WORLD, &P);
    MatSetSizes(P, m_num_my_rows, n_num_my_rows, m, n);

    // -- Preallocate memory
    mat_preallocate_memory(&prolongator,&P,m,n,m_num_my_rows,n_num_my_rows,mat_type_code,is_reduced);

    MatSetFromOptions(P);

    // -- Assemble
    MatGetOwnershipRange(P,&Istart,&Iend1);
    QPETScProlongatorAssembler assembler(&prolongator,P,Istart,Iend1,is_reduced);
    prolongator.assemble_P_raw(&assembler);

    // -- Finalize assembly
    MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);

    return P;
}


/*@t ----------
 * \section{Remove redundat zeros from the Prolongators}
 *
 *
 *@c*/

Mat MatRemoveZeros(Mat A, double atol)
{
    PetscErrorCode ierr;
    PetscInt i,j;
    MatType mtype;
    PetscInt rmin, rmax1, cmin, cmax1;
    PetscInt mg,ng;

    PetscInt ncols;
    const PetscInt* icols;
    const PetscScalar* vals;
    PetscInt ncolsr;
    PetscInt* icolsr;
    PetscScalar* valsr;

    PetscInt* d_nnz;
    PetscInt* o_nnz;

    Mat Ar;
    Vec rvec, lvec;

    // -- construct new matrix
    ierr = MatGetType(A,&mtype);
    ierr = MatGetVecs(A,&rvec,&lvec);
    ierr = MatGetSize(A,&mg,&ng);
    ierr = VecGetOwnershipRange(lvec,&rmin,&rmax1);
    ierr = VecGetOwnershipRange(rvec,&cmin, &cmax1);
    ierr = VecDestroy(rvec);
    ierr = VecDestroy(lvec);

    ierr = MatCreate(PETSC_COMM_WORLD,&Ar);
    ierr = MatSetSizes(Ar,rmax1-rmin,cmax1-cmin,mg,ng);
    ierr = MatSetType(Ar,mtype);

    // -- allocate memory for nonzeros
    ierr = PetscMalloc((rmax1-rmin)*sizeof(PetscInt),&d_nnz);
    ierr = PetscMalloc((rmax1-rmin)*sizeof(PetscInt),&o_nnz);
    ierr = PetscMemzero(d_nnz,(rmax1-rmin)*sizeof(PetscInt));
    ierr = PetscMemzero(o_nnz,(rmax1-rmin)*sizeof(PetscInt));

    for (i = rmin; i < rmax1; ++i) {

        ierr = MatGetRow(A, i, &ncols, &icols, &vals);

        for (j = 0; j < ncols; ++j) {
            if (PetscAbsScalar(vals[j]) > atol) {
                if (icols[j] >= cmin && icols[j] < cmax1)
                    d_nnz[i-rmin]+=1;
                else
                    o_nnz[i-rmin]+=1;
            }
        }

        ierr = MatRestoreRow(A, i, &ncols, &icols, &vals);

    }
    ierr = MatSeqAIJSetPreallocation(Ar, 0, d_nnz);
    ierr = MatMPIAIJSetPreallocation(Ar, 0, d_nnz, 0, o_nnz);
    ierr = PetscFree(d_nnz);
    ierr = PetscFree(o_nnz);

    // -- set values into matrix
    for (i = rmin; i < rmax1; ++i) {

        ierr = MatGetRow(A, i, &ncols, &icols, &vals);
        
        ncolsr = 0;
        ierr = PetscMalloc(ncols*sizeof(PetscInt),&icolsr);
        ierr = PetscMalloc(ncols*sizeof(PetscScalar),&valsr);
        for (j = 0; j < ncols; ++j) {
            if (PetscAbsScalar(vals[j]) > atol) {
                icolsr[ncolsr] = icols[j];
                valsr[ncolsr]  = vals[j];
                ncolsr++;
            }
        }

        ierr = MatRestoreRow(A, i, &ncols, &icols, &vals);

        ierr = MatSetValues(Ar, 1, &i, ncolsr, icolsr, valsr, INSERT_VALUES);
        ierr = PetscFree(valsr);
        ierr = PetscFree(icolsr);

    }
    ierr = MatAssemblyBegin(Ar,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(Ar,MAT_FINAL_ASSEMBLY);

    return Ar;
}
