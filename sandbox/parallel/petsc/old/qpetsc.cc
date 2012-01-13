#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "mesh.h"

#include "qpetsc.h"


const char* qPetscErrorMessage(int errnum)
{
    const char* s = NULL;
    PetscErrorMessage(errnum, &s, PETSC_NULL);
    return s;
}


PetscViewer qPetscViewerASCIIOpen(const char* name)
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &viewer);
    return viewer;
}


PetscViewer qPetscViewerBinaryOpen(const char* name, int type)
{
    PetscFileMode file_mode[] = {
        FILE_MODE_WRITE,
        FILE_MODE_READ,
    };
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, file_mode[type], &viewer);
    return viewer;
}


int qPetscViewerSetFormat(PetscViewer viewer, int format)
{
    const PetscViewerFormat formats[] = {
        PETSC_VIEWER_ASCII_DEFAULT,
        PETSC_VIEWER_ASCII_MATLAB,
        PETSC_VIEWER_ASCII_DENSE,
        PETSC_VIEWER_ASCII_IMPL,
        PETSC_VIEWER_ASCII_INFO,
        PETSC_VIEWER_ASCII_INFO_DETAIL,
        PETSC_VIEWER_ASCII_COMMON,
        PETSC_VIEWER_ASCII_INDEX,
        PETSC_VIEWER_ASCII_SYMMODU
    };
    return PetscViewerSetFormat(viewer, formats[format]);
}


Vec qVecCreate()
{
    Vec x;
    VecCreate(PETSC_COMM_WORLD, &x);
    return x;
}


int qVecNorm(Vec x, int type, double* val)
{
    NormType norm_types[] = {
        NORM_INFINITY,
        NORM_1,
        NORM_2
    };
    return VecNorm(x, norm_types[type], val);
}


Vec qVecLoad(PetscViewer viewer, int outtype)
{
    VecType vec_types[] = { VECSEQ, VECMPI };
    Vec vec;
    VecLoad(viewer, vec_types[outtype], &vec);
    return vec;
}


Mat qMatLoad(PetscViewer viewer, int outtype)
{
    MatType mat_types[] = {
        MATSEQAIJ,   MATMPIAIJ,
        MATSEQBAIJ,  MATMPIBAIJ,
        MATSEQSBAIJ, MATMPISBAIJ,
        MATSEQBDIAG, MATMPIBDIAG,
        MATMPIROWBS,
        MATSEQDENSE, MATMPIDENSE
    };
    Mat mat;
    MatLoad(viewer, mat_types[outtype], &mat);
    return mat;
}


KSP qKSPCreate()
{
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    return ksp;
}


int qKSPSetOperators(KSP ksp, Mat A, Mat B, int sameflag)
{
    MatStructure mat_structures[] = {
        SAME_PRECONDITIONER,
        SAME_NONZERO_PATTERN,
        DIFFERENT_NONZERO_PATTERN
    };
    return KSPSetOperators(ksp, A, B, mat_structures[sameflag]);
}


int qKSPTrueResidualNorm(KSP ksp, int type, double* rnorm)
{
    NormType norm_types[] = {
        NORM_INFINITY,
        NORM_1,
        NORM_2
    };
    Vec residual;
    PetscErrorCode ierr;
    ierr = KSPBuildResidual(ksp, PETSC_NULL, 
                            PETSC_NULL, &residual); CHKERRQ(ierr);
    ierr = VecNorm(residual, norm_types[type], rnorm); CHKERRQ(ierr);
    ierr = VecDestroy(residual); CHKERRQ(ierr);
    return 0;
}


PC qKSPGetPC(KSP ksp)
{
    PC pc;
    KSPGetPC(ksp, &pc);
    return pc;
}


int PCSetCoordinatesFromMesh(PC pc, Mesh* mesh)
{
    // Currently, Prometheus doesn't copy correctly if ndf != 3
    if (mesh == NULL || mesh->get_ndf() != 3)
        return 0;

    double* coords = new double[mesh->numnp()*3];
    memset(coords, 0, mesh->numnp()*3*sizeof(double));
    for (int i = 0; i < mesh->numnp(); ++i)
        for (int j = 0; j < mesh->get_ndm(); ++j)
            coords[i*3+j] = mesh->x(j,i);
    PetscErrorCode ierr = PCSetCoordinates(pc, 3, coords); // FIXME
//    PetscErrorCode ierr;
    delete[] coords;
    return ierr;
}
