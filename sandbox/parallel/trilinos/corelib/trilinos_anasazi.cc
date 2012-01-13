/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include <cstdio>
#include <cstring>
#include <cmath>

#include "mesh.h"
#include "qpassembly_trilinos.h"
#include "qcomplex.h"

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziTypes.hpp"

#include "trilinos_anasazi.h"
#include "Amesos_Operator.h"

#include "AnasaziEpetraMultiVectorComplexAdapter.hpp"
#include "AnasaziEpetraOperatorComplexAdapter.hpp"
#include "Amesos_Operator_Complex.h"
#include "EpetraExt_CrsMatrixIn.h"

template<class ST, class MV, class OP>
void create_eigenproblem(
    Teuchos::RefCountPtr< Anasazi::BasicEigenproblem<ST,MV,OP> >& eigenproblem,
    Teuchos::RefCountPtr<OP>& invKshiftM_r,
    Teuchos::RefCountPtr<MV>& ivec_r, int nev, Teuchos::ParameterList* pl)
{

    // -- Extract parameters from list
    bool is_sym   = (pl->get("IsSymmetric",0)==1) ? true : false;

    // -- Setup the eigenproblem
    eigenproblem->setOperator(invKshiftM_r);
    eigenproblem->setInitVec(ivec_r);
    eigenproblem->setNEV(nev);
    eigenproblem->setHermitian(is_sym);

    // -- Signal set up
    bool ierr = eigenproblem->setProblem();
    assert(ierr == true);
}

template<class ST, class MV, class OP>
int compute_eigs_anasazi_general(OP* invKshiftM,
                                 MV* ivec,
                 int nev, Teuchos::ParameterList* pl,
                 std::vector<Anasazi::Value<ST> >& d, MV* v)
{
    int status;

    // -- Construct eigensolver
    Teuchos::RefCountPtr<Anasazi::SolverManager<ST,MV,OP> > eigensolver;

    // -- Construct RCP of Kshift^{-1}*M, M, ivec
    Teuchos::RefCountPtr<OP> invKshiftM_r = Teuchos::rcp(invKshiftM,false);
    Teuchos::RefCountPtr<MV> ivec_r = Teuchos::rcp(ivec,false);

    // -- Create the eigenproblem
    Teuchos::RefCountPtr< Anasazi::BasicEigenproblem<ST,MV,OP> > eigenproblem =
        Teuchos::rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>() );
    create_eigenproblem<ST,MV,OP>(eigenproblem, invKshiftM_r, ivec_r, nev, pl);

    // -- Create the output manager
    int vb    = pl->template get<int>("Verbosity");
    vb = 0x1 + 0x2 + 0x4 + 0x8 + 0x10 + 0x20;

    // -- Create the sort manager
    std::string which(pl->template get<std::string>("Sigma"));
    pl->set( "Which", which);

    // -- Create the eigen solver
    switch ( pl->template get<int>("Solver") ) {

      case 0:
        eigensolver = Teuchos::rcp(new Anasazi::BlockKrylovSchurSolMgr<ST,MV,OP>(eigenproblem, *pl ));
        break;

      case 1:
        eigensolver = Teuchos::rcp(new Anasazi::BlockDavidsonSolMgr<ST,MV,OP>(eigenproblem, *pl ));
        break;

      case 2:
        eigensolver = Teuchos::rcp(new Anasazi::LOBPCGSolMgr<ST,MV,OP>(eigenproblem, *pl ));
        break;

    }

    // -- Solve the eigenvalue problem, and save the return code
    Anasazi::ReturnType solverreturn = eigensolver->solve();

    // -- Check return code of the solver: Unconverged, Failed, or OK
    switch (solverreturn) {

      case Anasazi::Converged:
        status =  0;
        break;

      case Anasazi::Unconverged:
        status = -1;
        break;

    }

    // -- Get eigenvalues and eigenvectors
    if (status==0) {
        Anasazi::Eigensolution<ST, MV> sol = eigenproblem->getSolution();
        std::vector<int> index = sol.index;
        int numev = sol.numVecs;

        // -- Get eigenvectors
        Teuchos::RefCountPtr<MV>  evecs = sol.Evecs;
//std::cout << "numvecs:" << sol->numVecs() << "\n";
std::cout << "NUMVECS:" << evecs->GetNumberVecs() << "\n";
std::cout << "NUMVECSv:" << v->GetNumberVecs() << "\n";
        v->MvAddMv( Teuchos::ScalarTraits<ST>::one(), *(evecs.get()),
                    Teuchos::ScalarTraits<ST>::zero(),*(evecs.get()) );

        // -- Get eigenvalues
        d = sol.Evals;
    }

    return status;
}

// -- Anasazi for Real problem
int compute_eigs_anasazi(Epetra_Operator* Op,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector* v)
{
    typedef double ST;

    Anasazi::EpetraMultiVec* vz_EMV;
    Anasazi::EpetraMultiVec* ivec_EMV;
    Anasazi::EpetraOp* Op_EOP;

    // -- Extract parameters from list
    int blocksize = pl->get<int>("Block Size");
    int is_sym    = pl->get("IsSymmetric",0);
    int is_real   = pl->get("IsReal",1);
    int vec_size  = v->NumVectors();
//    int vec_size  = (is_sym & is_real) ? nev : 1*nev;

    std::vector<Anasazi::Value<ST> > dz(vec_size);
    std::vector<int> vind(vec_size); for(int i = 0; i < vec_size; ++i) vind[i]=i;
    const Epetra_Map* Map = &(Op->OperatorDomainMap());

    // -- Create Epetra MultiVector to store eigenvectors
    vz_EMV = new Anasazi::EpetraMultiVec(View, *v, vind);

    // -- Create operator for Kshift^{-1}*M, M
    Op_EOP = new Anasazi::EpetraOp(Teuchos::rcp(Op, false));

    // -- Create initial Epetra MultiVector
    ivec_EMV = new Anasazi::EpetraMultiVec(*Map,blocksize);
    ivec_EMV->Random();

    int status =
      compute_eigs_anasazi_general<ST, Anasazi::MultiVec<ST>, Anasazi::Operator<ST> >
                                  (Op_EOP, ivec_EMV, nev, pl, dz, vz_EMV);

    // -- Copy real part
    for (int j=0; j < nev; ++j){
        dr[j] = dz[j].realpart;
        di[j] = dz[j].imagpart;
    }

    // -- Free allocated stuff
    delete ivec_EMV;
    delete Op_EOP;
    delete vz_EMV;

    return status;
}

int compute_eigs_anasazi(Epetra_CrsMatrix* K,
                         Epetra_CrsMatrix* M,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector* v)
{
//    PCKSP_Operator_Complex invKM(K, M);
    // -- Default solver is MUMPS
    Amesos_Operator invKM(K, M, 6);

    // -- Compute eigs
    return compute_eigs_anasazi(&invKM, nev, pl, dr, di, v);
}

int compute_eigs_anasazi(Epetra_CrsMatrix* Kshift,
                         Epetra_CrsMatrix* M, double w0, int form,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector* v)
{
    int ierr = compute_eigs_anasazi(Kshift, M, nev, pl, dr, di, v);

    // -- Reverse spectral transformation depending on form
    undo_spectral_trans(nev, w0, 0, form, dr, di);

    return ierr;
}

// -- Anasazi for Complex problem
int compute_eigs_anasazi(Epetra_Operator_Complex* Op,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector_Complex* v)
{
    typedef dcomplex ST;

    Anasazi::EpetraMultiVecComplex* vz_EMV;
    Anasazi::EpetraMultiVecComplex* ivec_EMV;
    Anasazi::EpetraOpComplex* Op_EOP;

    // -- Extract parameters from list
    int blocksize = pl->get<int>("Block Size");
    int is_sym    = pl->get("IsSymmetric",0);
    int is_real   = pl->get("IsReal",1);
    int vec_size  = (is_sym & is_real) ? nev : 1*nev;

    std::vector<Anasazi::Value<ST> > dz(vec_size);
    std::vector<int> vind(vec_size); for(int i = 0; i < vec_size; ++i) vind[i]=i;
    const Epetra_Map* Map = &(Op->OperatorDomainMap());

    // -- Create Epetra MultiVector to store eigenvectors
    vz_EMV = new Anasazi::EpetraMultiVecComplex(View, *v, vind);

    // -- Create operator for Kshift^{-1}*M, M
    Op_EOP = new Anasazi::EpetraOpComplex(Teuchos::rcp(Op, false));

    // -- Create initial Epetra MultiVector
    ivec_EMV = new Anasazi::EpetraMultiVecComplex(*Map,blocksize);
    ivec_EMV->Random();

    int status =
      compute_eigs_anasazi_general<ST, Anasazi::MultiVec<ST>, Anasazi::Operator<ST> >
                                  (Op_EOP, ivec_EMV, nev, pl, dz, vz_EMV);

    // -- Copy real part
    for (int j=0; j < nev; ++j){
        dr[j] = dz[j].realpart;
        di[j] = dz[j].imagpart;
    }

    // -- Free allocated stuff
    delete ivec_EMV;
    delete Op_EOP;
    delete vz_EMV;

    return status;
}

int compute_eigs_anasazi(Epetra_CrsMatrix_Complex* K,
                         Epetra_CrsMatrix_Complex* M,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector_Complex* v)
{
//    PCKSP_Operator_Complex invKM(K, M);
    // -- Default solver is MUMPS
    Amesos_Operator_Complex invKM(K, M, 6);

    // -- Compute eigs
    return compute_eigs_anasazi(&invKM, nev, pl, dr, di, v);
}

int compute_eigs_anasazi(Epetra_CrsMatrix_Complex* Kshift,
                         Epetra_CrsMatrix_Complex* M, double w0, int form,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector_Complex* v)
{
    int ierr = compute_eigs_anasazi(Kshift, M, nev, pl, dr, di, v);

    // -- Reverse spectral transformation depending on form
    undo_spectral_trans(nev, w0, 0, form, dr, di);

    return ierr;
}

void undo_spectral_trans(int nev, double sr, double si, int form, double* dr, double*di)
{

    // -- Reverse spectral transformation depending on form
    dcomplex d[nev];
    dcomplex s(sr,si);
    copy_complex(dr, di, d, nev);
    for (int i = 0; i < nev; ++i) {
        if (form==0)
            d[i] = s + 1.0/d[i];
        else if (form==1)
            d[i] = sqrt(s*s + 1.0/d[i]);
    }
    copy_complex(d, dr, di, nev);

}

// -- Compute with mesh directly

int compute_eigs_anasazi(Mesh* mesh, double w0,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector_Complex* vz)
{

    Epetra_CrsMatrix_Complex* Kshift;
    Epetra_CrsMatrix_Complex* M;

    Kshift = Mesh_assemble_dRz_trilinos(mesh, 1.0, 0.0, -w0*w0, 1);
    M      = Mesh_assemble_dRz_trilinos(mesh, 0.0, 0.0,    1.0, 1);
    int status = compute_eigs_anasazi(Kshift, M, nev, pl, dr, di, vz);

    // -- Reverse spectral transformation depending on form
    int form = 1; // Squared form
    undo_spectral_trans(nev, w0, 0, form, dr, di);

    delete M;
    delete Kshift;

    return status;
}

int compute_eigs_anasazi(Mesh* mesh, double w0,
                 int nev, Teuchos::ParameterList* pl, double* dr, double* di,
                 Epetra_MultiVector* v)
{

    Epetra_CrsMatrix* Kshift;
    Epetra_CrsMatrix* M;

    Kshift = Mesh_assemble_dR_trilinos(mesh, 1.0, 0.0, -w0*w0, 1);
    M      = Mesh_assemble_dR_trilinos(mesh, 0.0, 0.0,    1.0, 1);
    int status = compute_eigs_anasazi(Kshift, M, nev, pl, dr, di, v);

    // -- Reverse spectral transformation depending on form
    int form = 1; // Squared form
    undo_spectral_trans(nev, w0, 0, form, dr, di);

    delete M;
    delete Kshift;

    return status;
}

