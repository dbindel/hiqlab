#include <iostream>
#include "Petsc_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_LocalMap.h"

#include "Epetra_Export.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_VectorOut.h"

#include "Epetra_MpiComm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Export.h"

#include "qpetsc_pc.h"

extern Epetra_MpiComm* HIQLAB_Comm;

using namespace std;

Petsc_Operator::Petsc_Operator(KSP ksp, Epetra_Map* dm, Epetra_Map* rm,
		    Epetra_CrsMatrix* B)
  : Label_("Petsc Operator"), ksp(ksp), M(B), is_petsc_same(1)
{
    dmap   = new Epetra_Map(*dm);
    rmap   = new Epetra_Map(*rm);
    set_petsc_map();
}

Petsc_Operator::Petsc_Operator(KSP ksp, Epetra_Map* dm,
		    Epetra_CrsMatrix* B)
  : Label_("Petsc Operator"), ksp(ksp), M(B), is_petsc_same(1)
{
    dmap   = new Epetra_Map(*dm);
    rmap   = new Epetra_Map(*dm);
    set_petsc_map();
}

Petsc_Operator::Petsc_Operator(KSP ksp, Epetra_CrsMatrix* B)
  : Label_("Petsc Operator"), ksp(ksp), M(B), is_petsc_same(1)
{
    dmap   = new Epetra_Map(B->OperatorDomainMap());
    rmap   = new Epetra_Map(B->OperatorRangeMap());
    set_petsc_map();
}

Petsc_Operator::~Petsc_Operator() 
{
    delete petsc_map;
    delete dmap;
    delete rmap;
    if (!is_petsc_same)
        delete petsc_map_r;
}

void Petsc_Operator::set_petsc_map()
{
    Mat Kshift;
    PetscInt Istart, Iend;
    int NumMyElements;
    int* MyGlobalElements;

    // -- Get range of indices on this process
    KSPGetOperators(ksp, &Kshift, PETSC_NULL, PETSC_NULL);
    MatGetOwnershipRange(Kshift, &Istart, &Iend);
    NumMyElements = Iend - Istart; 
    MyGlobalElements = new int[NumMyElements];
    for (int i = 0; i < NumMyElements; ++i)
        MyGlobalElements[i] = Istart + i;

    petsc_map   = new Epetra_Map(-1, NumMyElements, MyGlobalElements, 0, *HIQLAB_Comm);
    petsc_map_r = petsc_map;
    // -- Clean up
    delete[] MyGlobalElements;
}

void Petsc_Operator::SetMesh(Mesh* msh)
{
    mesh = msh;

    // -- Construct reduced map if matrix size and fdof are not equal
    int numid = mesh->get_numid();
    PetscInt m,n;
    Mat A;
    KSPGetOperators(ksp, &A, PETSC_NULL, PETSC_NULL);
    MatGetSize(A,&m,&n);

    if (numid!=m)
        SetReducedMap();

    // -- Call PCSetCoordinates for Prometheus if ndm==3
    if (mesh->get_ndm()==3) {
        PC pc;
        KSPGetPC(ksp, &pc);
        PCSetFromOptions(pc);
std::cout << "Setting PC Coordinates\n";
        PCSetCoordinatesFromMesh(pc, mesh);
    }
}

void Petsc_Operator::SetReducedMap()
{

    int ndf           = mesh->get_ndf();
    int NumMyElements = 0;
    int MinGID_p      = petsc_map->MinMyGID();
    int MaxGID_p      = petsc_map->MaxMyGID();
    int MinNode_p     = MinGID_p/ndf;
    int MaxNode_p     = MaxGID_p/ndf;
    int NumMyNodes_p  = petsc_map->NumMyElements()/ndf;
    int MyGlobalElements[petsc_map->NumMyElements()];

    // -- Look for nodes on this processor in PetscMap and 
    //    Identify DOFs which have ID>=0 assigned. Record this
    //    into the reduced mapping
    for (int i = 0; i < NumMyNodes_p; ++i)
        for (int j = 0; j < ndf; ++j)
            if (mesh->id(j, MinNode_p+i) >= 0) {
                MyGlobalElements[NumMyElements] = mesh->id(j, MinNode_p+i);
                NumMyElements++;
            }
      
    petsc_map_r = new Epetra_Map(-1, NumMyElements, MyGlobalElements, 0, *HIQLAB_Comm);
    is_petsc_same = 0;
}

void Petsc_Operator::SetFromOptions()
{
    PC pc;
    KSPGetPC(ksp,&pc);
    KSPSetFromOptions(ksp);
    PCSetFromOptions(pc);
}

void Petsc_Operator::reduced2full(Epetra_MultiVector* xtmp,
                                  Epetra_MultiVector* xtmp_f) const
{

    // -- If the mapping is the same then simply export
    if (is_petsc_same) {
        Epetra_Export Exporter(xtmp->Map(), *petsc_map);
        xtmp_f->Export(*xtmp, Exporter, Add);
    
        return;
    }

    // -- Put GIDs on correct process
    Epetra_MultiVector xtmp_r(*petsc_map_r,xtmp->NumVectors());
    Epetra_Export Exporter(xtmp->Map(), *petsc_map_r);
    xtmp_r.Export(*xtmp, Exporter, Add);

    // -- Copy data from xtmp_r to xtmp_f
    double* vals;
    int lda;
    int ndf = mesh->get_ndf();
    xtmp_r.ExtractView(&vals,&lda);
    int NumMyNodes    = petsc_map->NumMyElements()/ndf;
    int MinNode       = petsc_map->MinMyGID()     /ndf;
    for (int i = 0; i < NumMyNodes; ++i)
        for (int j = 0; j < ndf; ++j) {
            int ind = mesh->id(j, MinNode+i);
            if (ind >= 0) {
                for (int k = 0; k < xtmp->NumVectors(); ++k)
                    xtmp_f->ReplaceMyValue(ndf*i+j,k,
                               vals[lda*k+petsc_map_r->LID(ind)] );
            }
        }

}

void Petsc_Operator::full2reduced(Epetra_MultiVector* xtmp_f,
                                  Epetra_MultiVector* xtmp) const
{
    // -- If the mapping is the same then simply export
    if (is_petsc_same) {
        Epetra_Export Exporter(*petsc_map, xtmp->Map());
        xtmp->Export(*xtmp_f, Exporter, Add);
    
        return;
    }

    // -- Put GIDs on correct process
    Epetra_MultiVector xtmp_r(*petsc_map_r,xtmp_f->NumVectors());

    // -- Copy data from xtmp_f to xtmp_r
    double* vals;
    int lda;
    int ndf = mesh->get_ndf();
    xtmp_f->ExtractView(&vals,&lda);
    int NumMyNodes    = petsc_map->NumMyElements()/ndf;
    int MinNode       = petsc_map->MinMyGID()     /ndf;
    for (int i = 0; i < NumMyNodes; ++i)
        for (int j = 0; j < ndf; ++j) {
            int ind = mesh->id(j, MinNode+i);
            if (ind >= 0) {
                for (int k = 0; k < xtmp->NumVectors(); ++k)
                    xtmp_r.ReplaceGlobalValue(ind,k,
                             vals[lda*k+ndf*i+j] );
            }
        }
    // -- Export to xtmp
    Epetra_Export Exporter(*petsc_map_r, xtmp->Map());
    xtmp->Export(xtmp_r, Exporter, Add);

}

int Petsc_Operator::Apply(
		 const Epetra_MultiVector& X, 
		       Epetra_MultiVector& Y) const {
//    if (!X.Map().SameAs(OperatorDomainMap())) EPETRA_CHK_ERR(-1);
//    if (!Y.Map().SameAs(OperatorRangeMap())) EPETRA_CHK_ERR(-2);
    if (Y.NumVectors()!=X.NumVectors()) EPETRA_CHK_ERR(-3);

double* normx;
double* normy;
normx = new double[X.NumVectors()];
normy = new double[X.NumVectors()];
X.Norm2(normx);
for (int i = 0; i < X.NumVectors(); ++i){
std::cout << "x[" << i << "]:" << normx[i] << "\n";
}

    // -- Construct temporary vector
    int nvecs = X.NumVectors();
    Epetra_MultiVector* xtmp;
    Epetra_MultiVector* xtmp_f;
    Epetra_MultiVector* ytmp_f;
    Mat Kshift;
    Vec rhs[nvecs];
    Vec lhs[nvecs];

    // -- Initialize Y
    Y.PutScalar(0.0); 

    // -- Multiply by M if exists and set RHS
    if (M) {
       xtmp = new Epetra_MultiVector(X.Map(),X.NumVectors());
       M->Multiply(false, X, *xtmp);
    } else {
       xtmp = new Epetra_MultiVector(X);
    }
xtmp->Norm2(normx);
for (int i = 0; i < X.NumVectors(); ++i){
std::cout << "xtmp[" << i << "]:" << normx[i] << "\n";
}
    // -- Construct rhs petsc vector from xtmp;
    xtmp_f = new Epetra_MultiVector(*petsc_map,X.NumVectors());
    reduced2full(xtmp, xtmp_f);

    KSPGetOperators(ksp, &Kshift, PETSC_NULL, PETSC_NULL);
    for (int i = 0; i < nvecs; ++i) 
        MatGetVecs(Kshift, &(rhs[i]), PETSC_NULL);
    Trilinos2PetscMultiVector(xtmp_f, rhs); //FIXME: assumed same distribution

    // -- Construct lhs petsc vector from matrix
    for (int i = 0; i < nvecs; ++i) 
        MatGetVecs(Kshift, &(lhs[i]), PETSC_NULL);

    // -- Solve for each rhs
    int ierr = 0;
    for (int i = 0; i < nvecs; ++i)
        KSPSolve(ksp,rhs[i],lhs[i]);

    // -- Copy into Y
    ytmp_f = new Epetra_MultiVector(*petsc_map,X.NumVectors());
    Petsc2TrilinosMultiVector(lhs, ytmp_f); //FIXME: assumed same distribution
    full2reduced(ytmp_f, &Y);

    // -- Show results
    if ( X.Comm().MyPID() == 0 ) {
        PetscReal norm;
        PetscInt  its;
        KSPGetResidualNorm(ksp,&norm);
        KSPGetIterationNumber(ksp,&its);
        std::cout << "Solver performed " << its
             << " iterations.\n";
        std::cout << "Norm of the true residual = " << norm << std::endl;
    }

Y.Norm2(normy);
for (int i = 0; i < X.NumVectors(); ++i){
std::cout << "y[" << i << "]:" << normy[i] << "\n";
}

    // -- Delete temporaries
delete normx;
delete normy;
    delete xtmp;
    delete xtmp_f;
    delete ytmp_f;
    for (int i = 0; i < nvecs; ++i) {
        VecDestroy(lhs[i]);
        VecDestroy(rhs[i]);
    }

    return(ierr);
}

void Trilinos2PetscMultiVector(Epetra_MultiVector* xtmp_f, Vec* rhs)
{
    //FIXME: assumed same distribution
    int NumVecs = xtmp_f->NumVectors();
    int NumMyElements = xtmp_f->MyLength();
    int MaxMyGID= xtmp_f->Map().MaxMyGID();
    int MinMyGID= xtmp_f->Map().MinMyGID();

    // -- Extract values and copy
    double* vals;
    int lda;
    xtmp_f->ExtractView(&vals,&lda);
    for (int i = 0; i < NumVecs; ++i) 
        for (int j = 0; j < NumMyElements; ++j) 
            VecSetValue(rhs[i],MinMyGID+j,vals[lda*i+j],INSERT_VALUES);

}

void Petsc2TrilinosMultiVector(Vec* lhs, Epetra_MultiVector* xtmp_f)
{
    //FIXME: assumed same distribution
    int NumVecs = xtmp_f->NumVectors();
    int NumMyElements = xtmp_f->MyLength();
    int MaxMyGID= xtmp_f->Map().MaxMyGID();
    int MinMyGID= xtmp_f->Map().MinMyGID();

    // -- Extract values and copy
    int ix[NumMyElements];
    PetscScalar vals[NumMyElements];
    for (int i = 0; i < NumMyElements; ++i) ix[i] = MinMyGID+i;
    for (int i = 0; i < NumVecs; ++i) {
        VecGetValues(lhs[i],NumMyElements,ix,vals);
        for (int j = 0; j < NumMyElements; ++j)
            xtmp_f->ReplaceMyValue(j, i, vals[j]);
    }

}
