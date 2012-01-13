/* HiQLab
 * Copyright (c): Regents of the University of California
 *
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
    myid.resize(n + mesh->numglobals());
    if (reduced == QP_ASSEMBLE_REMOVE_DIRICHLET) {
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
}


int QPETScAssembler::get_Vlocal(QMatrix<double>& Vl, int* eltid, int n)
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


void QPETScAssembler::apply_identity_dirichlet(int* eltid, int n, double* Ke)
{
    for (int i = 0; i < n; ++i) {
        if (mesh->id(eltid[i]) < 0) {
            for (int j = 0; j < n; ++j) {
                Ke[j*n+i] = 0;
                Ke[i*n+j] = 0;
            }
        }
    }
}


void QPETScAssembler::add_identity_dirichlet()
{
    int n = mesh->numnp() * mesh->get_ndf();
    for (int i = 0; i < n; ++i)
        if (mesh->id(i) < 0)
            MatSetValue(matrix, i, i, 1, ADD_VALUES);
}


void QPETScAssembler::add(int* eltid, int n, dcomplex* Ke)
{
    fprintf(stderr, "Complex PETSc assembler not yet implemented\n");
}


void QPETScAssembler::add(int* eltid, int n, double* Ke)
{
    map_idvec(eltid,n);
    int ng = mesh->numglobals();

    QMatrix<double> Ke1(n,n);
    Ke1 = Ke;
    if (reduced == QP_ASSEMBLE_IDENTITY_DIRICHLET)
        apply_identity_dirichlet(eltid, n, Ke1.data);

    if (ng == 0) {
        MatSetValues(matrix, n, &(myid[0]), n, &(myid[0]), 
                     Ke1.data, ADD_VALUES);
        return;
    }

    QMatrix<double>   Vl(n,ng);
    if (!get_Vlocal(Vl,eltid,n)) {
        MatSetValues(matrix, n, &(myid[0]), n, &(myid[0]), 
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


    MatSetValues(matrix, n+ng, &(myid[0]), n+ng, &(myid[0]),
                 Ke2.data, ADD_VALUES);
}


int* Mesh_petsc_nnz(Mesh* mesh, int blocked)
{
    int  numnp = mesh->numnp();
    int  ndf   = blocked ? 1 : mesh->get_ndf();
    int* nnz = new int[mesh->numnp()*ndf];
    vector<int> adj;

    for (int i = 0; i < numnp; ++i) {

        // Accumulate adjacency list
        adj.clear();
        for (int j = 0; j < mesh->get_nne(i); ++j) 
            for (int k = 0; k < mesh->get_nen(j); ++k) 
                if (mesh->ix(k,j) != i)
                    adj.push_back(mesh->ix(k,j));
        sort(adj.begin(), adj.end());

        // Count unique elements
        int count = 1;
        for (int j = 1; j < adj.size(); ++j)
            if (adj[j] != adj[j-1])
                ++count;

        for (int j = 0; j < ndf; ++j)
            nnz[i*ndf+j] = count*ndf;
    }
    return nnz;
}


Mat Mesh_assemble_dR_petsc(Mesh* mesh, double cx, double cv, double ca,
                           int mat_type_code)
{
    Mat A;
    int ndf = mesh->get_ndf();
    int n = mesh->numnp()*ndf; // FIXME: No account of globals

    if (mat_type_code == 0) {
        int* nnz = Mesh_petsc_nnz(mesh, 0);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, n,n, 0,nnz, &A);
        delete[] nnz;

    } else if (mat_type_code == 1) {
        int* nnz = Mesh_petsc_nnz(mesh, 1);
        MatCreateSeqBAIJ(PETSC_COMM_WORLD, ndf, n,n, 0,nnz, &A);
        delete[] nnz;

    } else if (mat_type_code == 2) {
        int* nnz = Mesh_petsc_nnz(mesh, 0);
        MatCreateMPIAIJ(PETSC_COMM_WORLD,
                        n,n,n,n,
                        PETSC_DEFAULT, nnz,
                        PETSC_DEFAULT, PETSC_NULL,
                        &A);
        delete[] nnz;

    } else if (mat_type_code == 3) {
        int* nnz = Mesh_petsc_nnz(mesh, 1);
        MatCreateMPIBAIJ(PETSC_COMM_WORLD, ndf,
                         n,n,n,n,
                         PETSC_DEFAULT, nnz,
                         PETSC_DEFAULT, PETSC_NULL,
                         &A);
        delete[] nnz;

    } else if (mat_type_code == 4) {
        int* nnz = Mesh_petsc_nnz(mesh, 0);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, n,n, 0,nnz, &A);
        MatSetOption(A, MAT_SYMMETRIC);
        delete[] nnz;

    } else if (mat_type_code == 5) {
        int* nnz = Mesh_petsc_nnz(mesh, 1);
        MatCreateSeqSBAIJ(PETSC_COMM_WORLD, ndf, n,n, 0,nnz, &A);
        delete[] nnz;

    } else if (mat_type_code == 6) {
        int* nnz = Mesh_petsc_nnz(mesh, 0);
        MatCreateMPIAIJ(PETSC_COMM_WORLD,
                        n,n,n,n,
                        PETSC_DEFAULT, nnz,
                        PETSC_DEFAULT, PETSC_NULL,
                        &A);
        MatSetOption(A, MAT_SYMMETRIC);
        delete[] nnz;

    } else if (mat_type_code == 7) {
        int* nnz = Mesh_petsc_nnz(mesh, 1);
        MatCreateMPISBAIJ(PETSC_COMM_WORLD, ndf,
                         n,n,n,n,
                         PETSC_DEFAULT, nnz,
                         PETSC_DEFAULT, PETSC_NULL,
                         &A);
        delete[] nnz;
    }

    MatSetFromOptions(A);
    QPETScAssembler assembler(A, mesh, QP_ASSEMBLE_IDENTITY_DIRICHLET);
    mesh->assemble_dR(&assembler, cx, cv, ca);
    assembler.add_identity_dirichlet();
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    return A;
}


Vec Mesh_assemble_R_petsc(Mesh* mesh)
{
    Vec x;
    int ndf = mesh->get_ndf();
    int n = mesh->numnp() * ndf;
    mesh->assemble_R();
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, n);
    VecSetFromOptions(x);

    for (int i = 0; i < n; ++i) {
        if (mesh->id(i) >= 0)
            VecSetValue(x, i, mesh->f(i), INSERT_VALUES);
        else
            VecSetValue(x, i, 0, INSERT_VALUES);
    }

    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    return x;
}
