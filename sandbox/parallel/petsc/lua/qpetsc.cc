#include <iostream>
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "mesh.h"

#include "qpetsc.h"


int get_MPI_Comm_size(MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    return size;
}

int get_MPI_Comm_rank(MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

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

int qPetscOptionsPrint()
{
    return PetscOptionsPrint(stdout);
}

int PetscOptionsSetOption(const std::string& iname, const std::string& value)
{

    int narg = 4;
    char** args;
    args = (char **)malloc( narg*sizeof(char *) );

    // -- Allocate memory for the passing array
    int ilen = strlen(iname.c_str());
    int vlen = strlen(value.c_str());
    args[0] = (char *)malloc( (ilen+1)*sizeof(char) );
    args[1] = (char *)malloc( (vlen+1)*sizeof(char) );
    args[2] = (char *)malloc( (ilen+1)*sizeof(char) );
    args[3] = (char *)malloc( (vlen+1)*sizeof(char) );

    // -- Copy arguments
    iname.copy(args[0],ilen);
    value.copy(args[1],vlen);
    iname.copy(args[2],ilen);
    value.copy(args[3],vlen);

    // -- Append last character
    args[0][ilen] = '\0';
    args[1][vlen] = '\0';
    args[2][ilen] = '\0';
    args[3][vlen] = '\0';

    int ierr = PetscOptionsInsert(&narg,&args,(char *)0);

    // -- clean up
    for (int i = 0; i < narg; ++i)
        free(args[i]);
    free(args);
//    delete args;

    return ierr;
}

Vec qVecCreate()
{
    Vec x;
    VecCreate(PETSC_COMM_WORLD, &x);
    return x;
}


double qVecNorm(Vec x, int type)
{

    double val;
    NormType norm_types[] = {
        NORM_INFINITY,
        NORM_1,
        NORM_2
    };
    VecNorm(x, norm_types[type], &val);
    return val;
}

int qVecScale(Vec x, double sr, double si)
{
#ifdef PETSC_USE_COMPLEX
    dcomplex s(sr,si);
    return VecScale(x,s);
#else
    return VecScale(x,sr);
#endif
}

int qVecSet(Vec x, double sr, double si)
{
#ifdef PETSC_USE_COMPLEX
    dcomplex s(sr,si);
    return VecSet(x,s);
#else
    return VecSet(x,sr);
#endif
}

int qVecAXPBY(Vec y, double a, double b, Vec x)
{
    return VecAXPBY(y, a, b, x);
}

int qVecAXPY (Vec y, double a,           Vec x)
{
    return VecAXPY(y, a, x);
}

int qVecSetRandom(Vec x)
{
    return VecSetRandom(x,PETSC_NULL);
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

int qMatScale(Mat m, double sr, double si)
{
#ifdef PETSC_USE_COMPLEX
    dcomplex s(sr,si);
    return MatScale(m,s);
#else
    return MatScale(m,sr);
#endif
}

Vec MatGetVecX(Mat m)
{
    Vec vec;
    MatGetVecs(m, &vec, PETSC_NULL);
    return vec;
}

Vec MatGetVecY(Mat m)
{
    Vec vec;
    MatGetVecs(m, PETSC_NULL, &vec);
    return vec;
}

Vec qMatGetRowSum(Mat m)
{
    Vec vec;
    MatGetVecs(m, PETSC_NULL, &vec);
    MatGetRowSum(m,vec);
    return vec;
}

Vec qMatGetDiagonal(Mat m)
{
    Vec vec;
    MatGetVecs(m, PETSC_NULL, &vec);
    MatGetDiagonal(m,vec);
    return vec;
}

int qMatPrintInfo(Mat m, int itype)
{
    int mpicommrank;
    PetscErrorCode ierr;
    MatInfo matinfo;
    MatInfoType info_types[] = {
        MAT_LOCAL,
        MAT_GLOBAL_MAX,
        MAT_GLOBAL_SUM
    };
    ierr = MatGetInfo(m, info_types[itype], &matinfo);

    MPI_Comm_rank(PETSC_COMM_WORLD,&mpicommrank);
    if (itype!=0) {
        if (mpicommrank==0) {
            std::cout << "Global  (rows,columns):(" << matinfo.rows_global  << "," << matinfo.columns_global << ")\n";
            std::cout << "Local   (rows,columns):(" << matinfo.rows_local   << "," << matinfo.columns_local  << ")\n";
            std::cout << "Blocksize             : " << matinfo.block_size   << "\n";
            std::cout << "Nz_alloc (used,unneed): " << matinfo.nz_allocated << "(" << matinfo.nz_used << "," 
                                                                                    << matinfo.nz_unneeded << ")\n";
            std::cout << "Memory(byte)          : " << matinfo.memory       << "\n";
            std::cout << "Mallocs               : " << matinfo.mallocs      << "\n";
        }
    } else {
            std::cout << "[" << mpicommrank << "]Global  (rows,columns):(" << matinfo.rows_global  << "," << matinfo.columns_global << ")\n";
            std::cout << "[" << mpicommrank << "]Local   (rows,columns):(" << matinfo.rows_local   << "," << matinfo.columns_local  << ")\n";
            std::cout << "[" << mpicommrank << "]Blocksize             : " << matinfo.block_size   << "\n";
            std::cout << "[" << mpicommrank << "]Nz_alloc (used,unneed): " << matinfo.nz_allocated << "(" << matinfo.nz_used << "," 
                                                                                    << matinfo.nz_unneeded << ")\n";
            std::cout << "[" << mpicommrank << "]Memory(byte)          : " << matinfo.memory       << "\n";
            std::cout << "[" << mpicommrank << "]Mallocs               : " << matinfo.mallocs      << "\n";
    }

    return ierr;
}

Mat qMatConvert(Mat m, const char* mattype)
{
    Mat A;
    PetscErrorCode ierr;
    ierr = MatConvert(m,mattype,MAT_INITIAL_MATRIX,&A);
    return A;
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


int qKSPSetType(KSP ksp, const char* ksptype)
{
    return KSPSetType(ksp,ksptype);
}


PC qPCCreate()
{
    PC pc;
    PCCreate(PETSC_COMM_WORLD, &pc);
    return pc;
}

int qPCSetOperators(PC pc, Mat A, Mat B, int sameflag)
{
    MatStructure mat_structures[] = {
        SAME_PRECONDITIONER,
        SAME_NONZERO_PATTERN,
        DIFFERENT_NONZERO_PATTERN
    };
    return PCSetOperators(pc, A, B, mat_structures[sameflag]);
}

int qPCSetType(PC pc, const char* pctype)
{
    return PCSetType(pc,pctype);
}
