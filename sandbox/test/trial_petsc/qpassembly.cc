/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 * $Id: qpassembly.cc,v 1.5 2006/06/09 04:54:24 dbindel Exp $
 */
 
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>

#include "qassembly.h"
#include "qpassembly.h"
#include "mesh.h"
#include "petscksp.h"

#define MAXNEN 64

using namespace std;

QPETScAssembler::QPETScAssembler(Mat matrix, Mesh* mesh,
		                 int reduced) :
	matrix(matrix), mesh(mesh), reduced(reduced)
{
}

QPETScAssembler::~QPETScAssembler()
{
}

void QPETScAssembler::map_idvec(int* eltid, int n)
{
    // -- Apply filter to select on processor entries
    PetscInt Istart, Iend;
    MatGetOwnershipRange(matrix, &Istart, &Iend);

    myid.resize(n + mesh->numglobals());
    if (reduced) {
        for (int i = 0; i < n; ++i)
            myid[i] = mesh->id(eltid[i]);
        for (int i = 0; i < mesh->numglobals(); ++i)
            myid[n+i] = mesh->globalid(i);
    } else {
        for (int i = 0; i < n; ++i)
            myid[i] = eltid[i];
        for (int i = 0; i < mesh->numglobals(); ++i)
            myid[n+i] = mesh->iglobal(i);
    }

    myid_l.resize(n + mesh->numglobals());
    for (int i = 0; i < n+mesh->numglobals(); ++i) {
	if (Istart<=myid[i] && myid[i]<Iend)
	    myid_l[i] = myid[i];
        else
	    myid_l[i] = -3; //FIXME: since -1, -2 used for BC
    } 
}

int QPETScAssembler::get_Vlocal(QMatrix<double>& Vl, 
		                 int* eltid, int n)
{
    int ng = mesh->numglobals();
    int nzflag = 0;
    for (int j = 0; j < ng; ++j)
        for (int i = 0; i < n; ++i) {
             Vl(i,j) = mesh->shapeg(eltid[i],j);
             nzflag = nzflag || (Vl(i,j) != 0);
        }
    return nzflag;
}

void QPETScAssembler::add(int* eltid, int n, double* Ke)
{
/*
    map_idvec(eltid,n);
    int ng = mesh->numglobals();

    QMatrix<double> Ke1(n,n);
    Ke1 = Ke;
//    if (reduced == QP_ASSEMBLE_IDENTITY_DIRICHLET)
//        apply_identity_dirichlet(eltid, n, Ke1.data);

    if (ng == 0) {
        MatSetValues(matrix, n, &(myid_l[0]), n, &(myid[0]), 
                     Ke1.data, ADD_VALUES);
        return;
    }

    QMatrix<double>   Vl(n,ng);
    if (!get_Vlocal(Vl,eltid,n)) {
        MatSetValues(matrix, n, &(myid_l[0]), n, &(myid[0]), 
                     Ke1.data, ADD_VALUES);
        return;
    }

    QMatrix<double> Ke2(n+ng,n+ng);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            Ke2(i,j) = Ke1(i,j);
    
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < ng; ++j)
            for (int l = 0; l < n; ++l)
                Ke2(i,j+n) += Ke2(i,l) * Vl(l,j);

    for (int i = 0; i < ng; ++i)
        for (int j = 0; j < n; ++j)
            for (int l = 0; l < n; ++l)
                Ke2(i+n,j) += Vl(l,i) * Ke2(l,j);

    for (int i = 0; i < ng; ++i)
        for (int j = 0; j < ng; ++j)
            for (int l = 0; l < n; ++l)
                Ke2(i+n,j+n) += Ke2(i+n,l) * Vl(l,j);


    MatSetValues(matrix, n+ng, &(myid_l[0]), n+ng, &(myid[0]),
                 Ke2.data, ADD_VALUES);
*/
}
void QPETScAssembler::add(int* eltid, int n, dcomplex* Ke)
{
    map_idvec(eltid,n);
    int ng = mesh->numglobals();

    QMatrix<dcomplex> Ke1(n,n);
    Ke1 = Ke;
//    if (reduced == QP_ASSEMBLE_IDENTITY_DIRICHLET)
//        apply_identity_dirichlet(eltid, n, Ke1.data);
    if (ng == 0) {
        MatSetValues(matrix, n, &(myid_l[0]), n, &(myid[0]), 
                     Ke1.data, ADD_VALUES);
        return;
    }

    QMatrix<double>   Vl(n,ng);
    if (!get_Vlocal(Vl,eltid,n)) {
        MatSetValues(matrix, n, &(myid_l[0]), n, &(myid[0]), 
                     Ke1.data, ADD_VALUES);
        return;
    }

    QMatrix<dcomplex> Ke2(n+ng,n+ng);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            Ke2(i,j) = Ke1(i,j);
    
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < ng; ++j)
            for (int l = 0; l < n; ++l)
                Ke2(i,j+n) += Ke2(i,l) * Vl(l,j);

    for (int i = 0; i < ng; ++i)
        for (int j = 0; j < n; ++j)
            for (int l = 0; l < n; ++l)
                Ke2(i+n,j) += Vl(l,i) * Ke2(l,j);

    for (int i = 0; i < ng; ++i)
        for (int j = 0; j < ng; ++j)
            for (int l = 0; l < n; ++l)
                Ke2(i+n,j+n) += Ke2(i+n,l) * Vl(l,j);


    MatSetValues(matrix, n+ng, &(myid_l[0]), n+ng, &(myid[0]),
                 Ke2.data, ADD_VALUES);
}

void Mesh_petsc_my_nnz(Mesh* mesh, double cx, double cv, double ca, int is_reduced,
		       PetscInt Istart, PetscInt Iend, int* d_nnz, int* o_nnz)
{

    int* ir;
    int* jc;
    double* Ax;
    int NumGlobalElements = (is_reduced) ? mesh->get_numid() : 
	                                   mesh->get_ndf()*mesh->numnp();
    int NumMyElements     = Iend - Istart;
    int Base              = Istart;
    int MinMyGID          = Istart;
    int MaxMyGID          = Iend - 1;

    // -- Construct assembler for local process(range restricted)
    CoordMatrix assembler(NumGlobalElements);
    assembler.restrict_range(MinMyGID, MaxMyGID);
    mesh->assemble_dR(&assembler, cx, cv, ca, is_reduced);

    // -- Compute nnz structure of local process
    assembler.pack();
    int ncoord = assembler.get_ncoord();
    ir = new int[NumGlobalElements+1];
    jc = new int[ncoord];
    Ax = new double[ncoord];
    memset(ir, 0, (NumGlobalElements+1) * sizeof(int));
    memset(jc, 0,                ncoord * sizeof(int));
    memset(Ax, 0,                ncoord * sizeof(double));
    memset(d_nnz, 0,      NumMyElements * sizeof(int));
    memset(o_nnz, 0,      NumMyElements * sizeof(int));
    assembler.to_sparse_row(ir, jc, Ax, NULL);
    for (int i = 0; i < NumMyElements; ++i)
	for (int j = 0; j < ir[Base+i+1]-ir[Base+i]; ++j)
            if (MinMyGID <= jc[ir[Base+i]+j] && MaxMyGID >= jc[ir[Base+i]+j])
                d_nnz[i]++;
    for (int i = 0; i < NumMyElements; ++i)
        o_nnz[i] = ir[Base+i+1]-ir[Base+i] - d_nnz[i];

    // -- Clean up
    delete ir, jc, Ax;

}

Mat Mesh_assemble_dR_petsc(Mesh* mesh, double cx, double cv, double ca,
                           int mat_type_code, int is_reduced)
{
    Mat A;
    PetscInt Istart, Iend;
    int* d_nnz;
    int* o_nnz;
    int numid= mesh->get_numid();
    int ndf  = mesh->get_ndf();
    
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, numid, numid);
    MatSetType(A, MATAIJ);
    MatGetOwnershipRange(A, &Istart, &Iend);
    d_nnz = new int[Iend-Istart];
    o_nnz = new int[Iend-Istart];

    if (mat_type_code == 0) {
        Mesh_petsc_my_nnz(mesh,cx,cv,ca,is_reduced,Istart,Iend,d_nnz,o_nnz);
	MatSetUpPreallocation(A);
	MatSeqAIJSetPreallocation(A, 0, d_nnz);

    } else if (mat_type_code == 2) {
std::cout << "Creating MPIAIJ\n";
        Mesh_petsc_my_nnz(mesh,cx,cv,ca,is_reduced,Istart,Iend,d_nnz,o_nnz);
	MatSetUpPreallocation(A);
	MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz);
    }

    // -- Clean up
    delete[] d_nnz;
    delete[] o_nnz;

    MatSetFromOptions(A);
    QPETScAssembler assembler(A, mesh, is_reduced);
    mesh->assemble_dR(&assembler, cx, cv, ca);
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    return A;
}


Vec Mesh_assemble_R_petsc(Mesh* mesh)
{
    double* fx;
    double* fi;
	
    Vec x;
    PetscInt Istart, Iend;
    int numid = mesh->get_numid();
    mesh->assemble_R();
    fx = new double[numid];
    fi = new double[numid];
    memset(fx, 0, numid * sizeof(double));
    memset(fi, 0, numid * sizeof(double));
    mesh->get_fz(fx,fi); //FIXME: Should not create ANOTHER copy of f
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, numid);
    VecSetFromOptions(x);
    VecGetOwnershipRange(x, &Istart, &Iend);

    for (int i = Istart; i < Iend; ++i) {
        VecSetValue(x, i, dcomplex(fx[i],fi[i]), INSERT_VALUES);
    }

    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

VecView(x,0);

    return x;
}
