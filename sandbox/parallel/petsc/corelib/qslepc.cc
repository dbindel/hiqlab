/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>

#include "mesh.h"
#include "qpassembly_petsc.h"
#include "qcomplex.h"

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "slepceps.h"

#include "qslepc.h"
#include "qpetsc_pc.h"

int compute_eigs_slepc(Mat Kshift, Mat M, int nev,
                double* dr, double* di, Vec vz)
{
    return compute_eigs_slepc(NULL, Kshift, M, nev, dr, di, vz);
}

int compute_eigs_slepc(Mesh* mesh, Mat Kshift, Mat M, int nev,
                double* dr, double* di, Vec vz)
{
    EPS eps;
    EPSType eps_type;
    PetscScalar kz;
    int nconv;

    // -- Set EPS
    EPSCreate(PETSC_COMM_WORLD, &eps);
    EPSSetOperators(eps, M, Kshift);
    EPSSetProblemType(eps, EPS_GNHEP);
    EPSSetDimensions(eps, nev, PETSC_DECIDE);
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
    EPSSetType(eps,EPSKRYLOVSCHUR);
    EPSSetFromOptions(eps);

std::cout << "Don't know why it crashes?!\n";
    // -- Set coordinates if PC is Prometheus
        ST  st;
        KSP ksp;
        PC  pc;
        PC  pcc;
        PCType pctype;

        EPSGetST(eps,&st);
        STGetKSP(st,&ksp);
        KSPGetPC(ksp,&pc);
        PCGetType(pc,&pctype);
        KSPSetFromOptions(ksp);
        PCSetFromOptions(pc);

    if (mesh) {
        if (strcmp(pctype,PCPROMETHEUS)==0) {
std::cout << "Setting Prometheus coordinates\n";
           EPSSetUp(eps);
//           PCSetCoordinatesFromMesh(pc, mesh);
        }
    }
    KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

    // -- Solve EPS
    EPSSolve(eps);

    // -- Obtain results
    EPSGetConverged(eps,&nconv);
    nconv = (nconv <= nev) ? nconv : nev;
    if (nconv > 0) {
        for (int i = 0; i < nconv; ++i) {
            EPSGetEigenpair(eps, i, &kz, PETSC_NULL, PETSC_NULL, PETSC_NULL);
            dr[i] = PetscRealPart(kz);
            di[i] = PetscImaginaryPart(kz);
        }
    }

    // -- Clean up
    EPSDestroy(eps);

    return 0;
}

// -- Compute with mesh directly

int compute_eigs_slepc(Mesh* mesh, double w0,
                 int nev, double* dr, double* di,
                 Vec vz)
{

    dcomplex* d;
    d      = new dcomplex[nev];

    // -- Construct matrices
    int mat_type_code = 2;
    int is_reduced    = 1;
std::cout << "Assembling K\n";
    Mat Kshift = Mesh_assemble_dR_petsc(mesh, 1.0, 0.0, -w0*w0,
                           mat_type_code, is_reduced);
std::cout << "Assembling M\n";
    Mat M      = Mesh_assemble_dR_petsc(mesh, 0.0, 0.0, 1.0,
                           mat_type_code, is_reduced);
std::cout << "Starting computation\n";
    int status = compute_eigs_slepc(mesh, Kshift, M, nev, dr, di, vz);
    copy_complex(dr, di, d, nev);
std::cout << "Raw eigenvalues\n";
    for (int i = 0; i < nev; ++i) {
std::cout << d[i] << "\n";
        d[i] = sqrt(w0*w0 + 1.0/d[i]);
    }
    copy_complex(d, dr, di, nev);

    delete[] d;
    MatDestroy(M);
    MatDestroy(Kshift);

    return status;
}
