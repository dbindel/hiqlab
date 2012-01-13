#include <iostream>
#include "petsc_jdqz_structure.h"
#include "petsc_jdqz_util.h"
#include "qlapack_util.h"
#include "qpetscmg.h"

using std::complex;
using std::real;
using std::imag;

/*@t ----------
 * \section{JDQZ_State class}
 *@c*/

#define ME JDQZ_State


ME::ME(Mat A, Mat B, Vec v, int neigs, int minq, int maxq) 
 : A(A), B(B), neigs(neigs), minq(minq), maxq(maxq), ceigs(0), sizeq(0), orthotype(1)
{

    // -- Construct vectors
    Q  = new Vec[neigs];
    Z  = new Vec[neigs];
    KiZ= new Vec[neigs];
    V  = new Vec[maxq ];
    W  = new Vec[maxq ];
    Vt = new Vec[maxq ];

    for (int i = 0; i < neigs; ++i) {
        VecDuplicate(v, &Q[i]);
        VecDuplicate(v, &Z[i]);
        VecDuplicate(v, &KiZ[i]);
    }
    for (int i = 0; i < maxq; ++i) {
        VecDuplicate(v, &V[i]);
        VecDuplicate(v, &W[i]);
        VecDuplicate(v, &Vt[i]);
    }
    VecDuplicate(v,&tempv);

    // -- construct reduced matrices
    S  = new complex<double>[neigs*neigs];
    T  = new complex<double>[neigs*neigs];
    MA = new complex<double>[maxq*maxq];
    MB = new complex<double>[maxq*maxq];
    SA = new complex<double>[maxq*maxq];
    TB = new complex<double>[maxq*maxq];
    UR = new complex<double>[maxq*maxq];
    UL = new complex<double>[maxq*maxq];

    // -- zero out arrays
    memset(MA, 0, maxq*maxq*sizeof(complex<double>));
    memset(MB, 0, maxq*maxq*sizeof(complex<double>));

}


ME::~ME()
{
    // -- delete vectors
    for (int i = 0; i < neigs; ++i) {
        VecDestroy(Q[i]);
        VecDestroy(Z[i]);
        VecDestroy(KiZ[i]);
    }
    for (int i = 0; i < maxq; ++i) {
        VecDestroy(V[i]);
        VecDestroy(W[i]);
        VecDestroy(Vt[i]);
    }
    VecDestroy(tempv);

    delete[] Q;
    delete[] Z;
    delete[] KiZ;
    delete[] V;
    delete[] W;
    delete[] Vt;

    // -- delete arrays
    delete[] S;
    delete[] T;
    delete[] MA;
    delete[] MB;
    delete[] SA;
    delete[] TB;
    delete[] UR;
    delete[] UL;
}


void ME::initialize_arnoldi(int initq, complex<double> alpha, complex<double> beta, Vec v0)
{
    if (sizeq < initq) {
        complex<double> sigma[] = {alpha, beta};

        InitArnoldi(V, v0, A, B, sigma, initq);

        // -- Orthogonalize V with Q
        if (ceigs==0)
            ProjectedGS(V, initq, NULL, 0, NULL, NULL, 0, tempv, orthotype, 1, 3);
        else
            ProjectedGS(V, initq,Q, ceigs, NULL, NULL, 0, tempv, orthotype, 1, 3);

        // -- W0 = sigma(1)*A*V0 - sigma(0)*B*V0
        for (int i = 0; i < initq; ++i)
            HarmonicVector(A,B,sigma,V[i],W[i],tempv);

        // -- Orthogonalize W with Z
        if (ceigs==0)
            ProjectedGS(W, initq, NULL, 0, NULL, NULL, 0, tempv, orthotype, 1, 3);
        else
            ProjectedGS(W, initq,Z, ceigs, NULL, NULL, 0, tempv, orthotype, 1, 3);

        sizeq = initq;
    }
}


void ME::compute_projections()
{
    ProjectedOP(A,V,sizeq,W,sizeq,tempv,MA,maxq,0,0);
    ProjectedOP(B,V,sizeq,W,sizeq,tempv,MB,maxq,0,0);
}


void ME::compute_projected_eigenvalues(complex<double> alpha, complex<double> beta)
{
    assert(sizeq>0);
    // -- Solve projected problem
    //    UL'*MA*UR = SA   UL'*MB*UR = TB
    //std::cout << "Solving project system:"<< sizeq <<"\n";
    SortQZ(sizeq,UR,maxq,UL,maxq,SA,maxq,TB,maxq,MA,maxq,MB,maxq,
           alpha,beta,(sizeq>=maxq ? minq : 1) );

    // -- Compute approximate eigenpair and residual
    //std::cout << "Compute approximate eigenpair and residual\n";
    //    q = V(:,1:sizeq)*UR(:,1)
    //    z = W(:,1:sizeq)*UL(:,1)
    VecZeroEntries(Q[ceigs]);
    VecZeroEntries(Z[ceigs]);
    VecMAXPY(Q[ceigs],sizeq,UR,V);
    VecMAXPY(Z[ceigs],sizeq,UL,W);
}


void ME::expand_subspace(Vec tt, complex<double> alpha, complex<double> beta)
{
    Vec* ospace;

    // orthogonalize against [Q,V]
    ospace = new Vec[ceigs+sizeq];
    for (int i = 0; i < ceigs; ++i)
        ospace[i] = Q[i];
    for (int i = 0; i < sizeq; ++i)
        ospace[ceigs+i] = V[i];

    ProjectedGS(&tt, 1, ospace, ceigs+sizeq, &V[sizeq], NULL, 0, tempv, orthotype, 1, 3);

    // -- Harmonic vector w = sigma[1]*A*v - sigma[0]*B*v
    complex<double> sigma[] = {alpha, beta};
    HarmonicVector(A,B,sigma,V[sizeq],W[sizeq],tempv);

    // orthogonalize against [Q,V]
    for (int i = 0; i < ceigs; ++i)
        ospace[i] = Z[i];
    for (int i = 0; i < sizeq; ++i)
        ospace[ceigs+i] = W[i];

    ProjectedGS(&W[sizeq], 1, ospace, ceigs+sizeq, NULL, NULL, 0, tempv, orthotype, 1, 3);

    sizeq++;

    // -- Update projected matrix  MA = W'*A*V MB = W'*B*V
    ProjectedOP(A, V         ,sizeq-1,&W[sizeq-1],      1,tempv,MA,maxq,sizeq-1,      0);
    ProjectedOP(A,&V[sizeq-1],      1, W         ,  sizeq,tempv,MA,maxq,      0,sizeq-1);
    ProjectedOP(B, V         ,sizeq-1,&W[sizeq-1],      1,tempv,MB,maxq,sizeq-1,      0);
    ProjectedOP(B,&V[sizeq-1],      1, W         ,  sizeq,tempv,MB,maxq,      0,sizeq-1);

    // -- clean up
    delete[] ospace;
}


void ME::compute_residual_vector(Vec rt)
{
    complex<double> theta[] = {get_SA(0,0), get_TB(0,0)}; 
    HarmonicVector(A,B,theta,Q[ceigs],rt,tempv);
    if (ceigs > 0)
        ProjectedGS(&rt, 1, Z, ceigs, NULL, NULL, 0, tempv, orthotype, 0, 3);
}


void ME::expand_partial_schur()
{
    ceigs++;
    ProjectedOP(A, &(Q[ceigs-1]), 1, Z, ceigs, tempv, S, neigs, 0, ceigs-1);
    ProjectedOP(B, &(Q[ceigs-1]), 1, Z, ceigs, tempv, T, neigs, 0, ceigs-1);
}


void ME::remove_converged_vector()
{
    // -- V(:,1:sizeq-1) = V(:,1:sizeq)*Ur(:,2:sizeq);
    for (int i = 1; i < sizeq; ++i) {
        VecZeroEntries(Vt[i-1]);
        VecMAXPY(Vt[i-1],sizeq,&UR[maxq*i],V);
    }
    for (int i = 1; i < sizeq; ++i)
        VecCopy(Vt[i-1],V[i-1]);
    VecZeroEntries(V[sizeq-1]);

    // -- W(:,1:sizeq-1) = W(:,1:sizeq)*Ul(:,2:sizeq);
    for (int i = 1; i < sizeq; ++i) {
        VecZeroEntries(Vt[i-1]);
        VecMAXPY(Vt[i-1],sizeq,&UL[maxq*i],W);
    }
    for (int i = 1; i < sizeq; ++i)
        VecCopy(Vt[i-1],W[i-1]);
    VecZeroEntries(W[sizeq-1]);

    // -- Update MA, MB
    for (int j = 0; j < sizeq-1; ++j)
        for (int i = 0; i < sizeq-1; ++i) {
            MA[maxq*j+i] = SA[maxq*(j+1)+(i+1)];
            MB[maxq*j+i] = TB[maxq*(j+1)+(i+1)];
        }

    sizeq--;
}


void ME::restart()
{
    assert(sizeq==maxq);

    // -- V(:,1:minq) = V(:,1:sizeq)*Ur(:,1:minq);
    for (int i = 0; i < minq; ++i) {
        VecZeroEntries(Vt[i]);
        VecMAXPY(Vt[i],minq,&UR[maxq*i],V);
    }
    for (int i = 0; i < minq; ++i)
        VecCopy(Vt[i],V[i]);

    // -- W(:,1:minq) = W(:,1:sizeq)*Ul(:,1:minq);
    for (int i = 0; i < minq; ++i) {
        VecZeroEntries(Vt[i]);
        VecMAXPY(Vt[i],minq,&UL[maxq*i],W);
    }
    for (int i = 0; i < minq; ++i)
        VecCopy(Vt[i],W[i]);

    // -- Update MA, MB
    for (int j = 0; j < minq; ++j)
        for (int i = 0; i < minq; ++i) {
            MA[maxq*j+i] = SA[maxq*j+i];
            MB[maxq*j+i] = TB[maxq*j+i];
        }

    sizeq = minq;
}

#undef ME


/*@t ----------
 * \section{JDQZ_Solve class}
 *@c*/


#define ME JDQZ_Solve


ME::ME(Mat A, Mat B, PC pc)
 : A(A), B(B), pc(pc), Q(0), Z(0), KiZ(0), neigs(0)
{
    // -- Set up shell operator
    PetscInt ln,lm, msize, nsize;
    MatGetSize(A,&msize,&nsize);
    MatGetOwnershipRange(A,&lm,&ln);
    ln = ln - lm;
    lm = ln;
    MatShellABCreate(&AsB);
    MatShellABSetOperators(AsB,A,B);
    MatCreateShell(PETSC_COMM_WORLD,lm,ln,msize,nsize,(void *)AsB,&AsBmat);
    MatShellSetOperation(AsBmat,MATOP_MULT,(void (*)(void))MatShellABApply);
    MatShellSetOperation(AsBmat,MATOP_MULT_TRANSPOSE,(void (*)(void))MatShellABApplyTranspose);

    // -- Register and set up Projected PC preconditioner in KSP
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPGetPC(ksp,&ppc);

    PCRegisterProjectedPC();
    PCSetType(ppc,"projectedpc");
    PCSetOperators(ppc,AsBmat,AsBmat,DIFFERENT_NONZERO_PATTERN);
    PCProjectedPCSetPC(ppc,pc);

    KSPSetTolerances(ksp,1.e-8,PETSC_DEFAULT,PETSC_DEFAULT,60);
    KSPSetFromOptions(ksp);
}


ME::~ME()
{
    // -- Destroy KSP
    KSPDestroy(ksp);

    // -- Destroy shell matrix
    MatDestroy(AsBmat);

    // -- Destroy shell operator
    MatShellABDestroy(AsB);
}


void ME::set_projection_spaces(Vec* Q, Vec* Z, Vec* KiZ, int neigst)
{
    neigs = neigst;
    PCProjectedPCSetSpaces(ppc,Q,Z,KiZ,neigs);
    PCProjectedPCSetCurrentSize(ppc,1);
    PCProjectedPCSetComputedSize(ppc,0);
}


void ME::apply(Vec x, Vec y)
{
    KSPSolve(ksp,x,y);
}


void ME::set_computed_size(int n)
{
    PCProjectedPCSetComputedSize(ppc,n);
}


void ME::set_current_size(int n)
{
    PCProjectedPCSetCurrentSize(ppc,n);
}


void ME::set_shift(complex<double> talpha, complex<double> tbeta)
{
    alpha = talpha;
    beta  = tbeta;

    MatShellABSetShift(AsB,alpha,beta);
}


void ME::set_tolerances(double rtol, int maxits)
{
    KSPSetTolerances(ksp,rtol,PETSC_DEFAULT,PETSC_DEFAULT,maxits);
}


void ME::set_ksp_type(string ksptype)
{
    KSPSetType(ksp,ksptype.c_str());

}


int ME::get_computed_size()
{
    int n;
    PCProjectedPCGetComputedSize(ppc,&n);
    return n;
}


int ME::get_current_size()
{
    int n;
    PCProjectedPCGetCurrentSize(ppc,&n);
    return n;
}

#undef ME


/*@t ------------
 * \section{JDQZ_Tracker}
 *
 * Class to hold the JDQZ tracking
 *@c*/


#define ME JDQZ_Tracker

void ME::set_target(complex<double>* target, complex<double>* aim, 
                    complex<double>* approx, double nr)
{
    double crit = std::real(approx[0]/approx[1]);
    double crita= std::real(aim[0]/aim[1]);
    double critd= std::abs(crit-crita)/std::abs(crita);

    if (nr < resid_atol && critd < eigdist_rtol) {
        ScaleEig(&approx[0],&approx[1],1);
        target[0] = approx[0];
        target[1] = approx[1];
    } else {
        target[0] = aim[0];
        target[1] = aim[1];
    }
}


void ME::set_ksp_tolerance(JDQZ_Solve* jsolve)
{
    double tol = 1/pow(10,4*citer);
    tol = (tol > ksp_rtol) ? tol : ksp_rtol;

    jsolve->set_tolerances(tol,ksp_maxits);
}

#undef ME

/*@t ------------
 * \section{JDQZ_Data}
 *
 * Class to hold the JDQZ data obtained from computation.
 *@c*/


#define ME JDQZ_Data


ME::ME(int n, Vec v0) : nvalues(n)
{
    if (nvalues>0) {
        revecs = new Vec[n];
        levecs = new Vec[n];

        for (int i = 0; i < n; ++i) {
            VecDuplicate(v0,&(revecs[i]));
            VecDuplicate(v0,&(levecs[i]));
        }

        alpha  = new complex<double>[n];
        beta   = new complex<double>[n];
    }
}


ME::~ME()
{
    if (nvalues>0) {
        for (int i = 0; i < nvalues; ++i) {
            VecDestroy(revecs[i]);
            VecDestroy(levecs[i]);
        }
        delete[] revecs;
        delete[] levecs;

        delete[] alpha;
        delete[] beta;
    }
}

void ME::compute_evecs(double* SAr, double* SAi, int ldSA,
                       double* TBr, double* TBi, int ldTB,
                       Vec* Q, Vec* Z)
{
    int n = nvalues;

    if (n > 0) {
        complex<double>* SA = new complex<double>[n*n];
        complex<double>* TB = new complex<double>[n*n];
        complex<double>* VL = new complex<double>[n*n];
        complex<double>* VR = new complex<double>[n*n];

        // -- copy data into complex format
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
                SA[n*j+i] = complex<double>(SAr[ldSA*j+i],SAi[ldSA*j+i]);
                TB[n*j+i] = complex<double>(TBr[ldTB*j+i],TBi[ldTB*j+i]);
            } 

        // -- compute eigenvalues
        zggev(n,SA,n,TB,n,alpha,beta,VL,n,VR,n);

        // -- compute eigenvectors
        for (int i = 0; i < n; ++i) {
            VecZeroEntries(revecs[i]);
            VecZeroEntries(levecs[i]);

            VecMAXPY(revecs[i],n,&(VR[n*i]),Q);
            VecMAXPY(levecs[i],n,&(VL[n*i]),Z);
        }

        // -- clean up
        delete[] SA;
        delete[] TB;
        delete[] VL;
        delete[] VR;
    }
}


void ME::compute_evecs(complex<double>* SA, int ldSA, complex<double>* TB, int ldTB,
                       Vec* Q, Vec* Z)
{
    int n = nvalues;

    if (n>0) {
        complex<double>* VL = new complex<double>[n*n];
        complex<double>* VR = new complex<double>[n*n];

        // -- compute eigenvalues
        zggev(n,SA,ldSA,TB,ldTB,alpha,beta,VL,n,VR,n);

        // -- compute eigenvectors
        for (int i = 0; i < n; ++i) {
            VecZeroEntries(revecs[i]);
            VecZeroEntries(levecs[i]);

            VecMAXPY(revecs[i],n,&(VR[n*i]),Q);
            VecMAXPY(levecs[i],n,&(VL[n*i]),Z);
        }

        // -- clean up
        delete[] VL;
        delete[] VR;
    }
}


#undef ME


