#ifndef _QPETSC_H
#define _QPETSC_H

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscpc.h"
#include "mesh.h"

inline PetscViewer get_PETSC_VIEWER_STDOUT_WORLD()
{
    return PETSC_VIEWER_STDOUT_WORLD;
}

inline MPI_Comm get_PETSC_COMM_WORLD()
{
    return PETSC_COMM_WORLD;
}

int get_MPI_Comm_size(MPI_Comm comm);
int get_MPI_Comm_rank(MPI_Comm comm);

const char* qPetscErrorMessage(int errnum);

PetscViewer qPetscViewerASCIIOpen(const char* name);
PetscViewer qPetscViewerBinaryOpen(const char* name, int type);
int qPetscViewerSetFormat(PetscViewer viewer, int format);

int qPetscOptionsPrint();
int PetscOptionsSetOption(const std::string& iname, const std::string& value);
//int PetscOptionsSetOption(const std::string& iname);

Vec qVecCreate();
double qVecNorm(Vec x, int type);
int qVecScale(Vec x, double sr, double si = 0);
int qVecSet(Vec x, double sr, double si = 0);
Vec qVecLoad(PetscViewer viewer, int outtype);
int qVecAXPBY(Vec y, double a, double b, Vec x);
int qVecAXPY (Vec y, double a,           Vec x);
int qVecSetRandom(Vec x);

Mat qMatLoad(PetscViewer viewer, int outtype);
int qMatScale(Mat m, double sr, double si);
Vec MatGetVecX(Mat m);
Vec MatGetVecY(Mat m);
Vec qMatGetRowSum(Mat m);
Vec qMatGetDiagonal(Mat m);
int qMatPrintInfo(Mat m, int itype);
Mat qMatConvert(Mat m, const char* mattype);

KSP qKSPCreate();
int qKSPSetOperators(KSP ksp, Mat A, Mat B, int sameflag);
PC  qKSPGetPC(KSP ksp);
int qKSPSetType(KSP ksp, const char* ksptype);
int qKSPTrueResidualNorm(KSP ksp, int type, double* rnorm);

PC qPCCreate();
int qPCSetOperators(PC pc, Mat A, Mat B, int sameflag);
int qPCSetType(PC pc, const char* pctype);

#endif /* _QPETSC_H */
