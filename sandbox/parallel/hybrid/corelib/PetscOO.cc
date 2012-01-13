#include <string>
#include "PetscOO.h"


PetscOO::PetscOO()
{
/*    int argc = 2;
    char* argv[2] = {"mpirun","hiqlab"};
    static char* help = "Test package";
    char** argvp;
    argvp = argv;
    PetscInitialize(&argc,&argvp,(char *)0,help);*/
    KSPCreate(PETSC_COMM_WORLD, &ksp);
}

PetscOO::~PetscOO()
{
    KSPDestroy(ksp);
/*    PetscFinalize();*/
}

int PetscOO::SetOperators(Mat A, Mat P, MatStructure flag)
{
    return KSPSetOperators(ksp,A,P,flag);   
}

int PetscOO::SetPCType(PCType type)
{
    PC pc;
    KSPGetPC(ksp,&pc);
    return PCSetType(pc,type);
}

int PetscOO::SetPCType(PCType type, Mesh* mesh)
{
    PC pc;
    KSPGetPC(ksp,&pc);
    PCSetType(pc,type);
    return PCSetCoordinatesFromMesh(pc, mesh);
}

void PetscOO::set_from_options()
{
    // -- Set PC options
    PC pc;
    KSPGetPC(ksp,&pc);
    PCSetFromOptions(pc);
 
    // -- Set KSP options
    KSPSetFromOptions(ksp);
}

int PetscOO::SetOption(const std::string& iname, const std::string& value)
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
        delete args[i];
    delete args;

    return ierr;
}

int PetscOO::Solve(Vec b, Vec x)
{
    return KSPSolve(ksp,b,x);
}

double PetscOO::GetResidualNorm()
{
    PetscReal rnorm;
    KSPGetResidualNorm(ksp,&rnorm);
    return rnorm;
}

int PetscOO::GetIterationNumber()
{
    PetscInt its;
    KSPGetIterationNumber(ksp,&its);
    return its;
}
