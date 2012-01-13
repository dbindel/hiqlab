#include <complex>
#include <iostream>
#include "petsc_jdqz_util.h"
#include "petsc_jdqz_structure.h"
#include "petsc_jdqz.h"

using std::complex;

JDQZ_Data* JDQZ(Mat A, Mat B, PC pc, Vec v0, 
          complex<double>* shift, int neigs, JDQZ_ParameterList* jpl)
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // -- Check and extract parameters
    jpl->check_parameters();
    int maxiter = jpl->maxiter;
    double atol = jpl->resid_atol;
    int initq   = jpl->initq;
    int minq    = jpl->minq;
    int maxq    = jpl->maxq;

    // -- Assert correctness of parameters
    int Am,An;
    MatGetSize(A,&Am,&An);
    assert(Am   >=neigs);

    // -- Set up data structure
    JDQZ_State jdstate(A, B, v0, neigs, minq, maxq);
    
    // -- Set up JD solver
    JDQZ_Solve jdsolve(A,B,pc);
    jdsolve.set_projection_spaces(jdstate.get_Q(),jdstate.get_Z(),jdstate.get_KiZ(),neigs);
    jdsolve.set_ksp_type(jpl->ksp_type);

    // -- Set up JD Tracker
    JDQZ_Tracker jtracker(maxiter);
    jtracker.set_resid_atol(jpl->tresid_atol);
    jtracker.set_eigdist_rtol(jpl->teigdist_rtol);
    jtracker.set_ksp_rtol(jpl->ksp_rtol);
    jtracker.set_ksp_maxits(jpl->ksp_maxits);

    // -- Allocate memory
    Vec  rt, tt;
    VecDuplicate(v0,&rt);
    VecDuplicate(v0,&tt);

    // -- Other variables
    complex<double> target[2], theta[2];
    complex<double> Sigma[2], sigma[2];
    Sigma[0] = shift[0];
    Sigma[1] = shift[1];
    ScaleEig(&(Sigma[0]),&(Sigma[1]),1);

    // -- Control flags
    int INITIATE = 1;
    int rKNOWN   = 0;
    int DETECTED = 0;

    // -- Analysis parameters
    double nr;

    if (rank==0) {
        std::cout << "neigs   :" << neigs   << "\n";
        std::cout << "maxiter :" << maxiter << "\n";
        std::cout << "atol    :" << atol    << "\n";
        std::cout << "initq   :" << initq   << "\n";
        std::cout << "Sigma[0]:" << Sigma[0]<< "\n";
        std::cout << "Sigma[1]:" << Sigma[1]<< "\n";
    }

    // -- JD iteration
    int niter = 0;
    while ( jdstate.get_num_converged() < neigs && niter < maxiter ) {

        niter++;
        jtracker.increment();
        // std::cout << "Iteration[" << niter << "]\n";

        if (INITIATE) {
            //std::cout << "Initiating\n";
            sigma[0] = Sigma[0];    
            sigma[1] = Sigma[1];    

            // -- Construct initial subspace
            jdstate.initialize_arnoldi(initq,sigma[0],sigma[1],v0);

            // -- Reset variables
            target[0] = sigma[0];
            target[1] = sigma[1];

            // -- Set projected matrices
            jdstate.compute_projections();

            INITIATE = 0;
            rKNOWN   = 0;
            DETECTED = 0;
        }

        // -- Solve correction equation and expand
        if (rKNOWN) {
            //std::cout << "Solving correction:" << target[0]/target[1] << "\n";
            // -- normalize right hand side
            VecScale(rt,1/nr);

            jdsolve.set_shift(target[0],target[1]);
            jtracker.set_ksp_tolerance(&jdsolve);
            jdsolve.apply(rt,tt);
            jdsolve.set_computed_size(jdsolve.get_computed_size()-1);

            rKNOWN = 0;

            // -- Expand space
            //std::cout << "Expanding subspace\n";
            jdstate.expand_subspace(tt,sigma[0],sigma[1]);
        }

        //std::cout << "Compute approximate eigenpair and residual\n";
        // -- Solve projected problem
        //    UL'*MA*UR = SA   UL'*MB*UR = TB
        //    Compute approximate eigenpair and residual
        //    q = V(:,1:sizeq)*UR(:,1)
        //    z = W(:,1:sizeq)*UL(:,1)
        jdstate.compute_projected_eigenvalues(target[0],target[1]);
        jdstate.compute_residual_vector(rt);
        rKNOWN = 1;

        theta[0] = jdstate.get_SA(0,0);
        theta[1] = jdstate.get_TB(0,0);
        VecNorm(rt,NORM_2,&nr);        

        if (rank==0){
            printf("Iter[%d]: Theta(%16.16f,%16.16f), Target(%16.16f,%16.16f): Resid:%16.16e\n",
                   niter, std::real(theta[0]/theta[1]),   std::imag(theta[0]/theta[1]),
                          std::real(target[0]/target[1]), std::imag(target[0]/target[1]),nr);
        }

        // -- If residual is small enough track
        jtracker.set_target(target,sigma,theta,nr);

        // -- Check convergence
        if (nr < atol) {
            //std::cout << "Converged\n";
            // -- Expand Schur form
            jdstate.expand_partial_schur();

            // -- Apply KiZ = Prec^{1}*Z
            jdsolve.set_current_size(jdstate.get_num_converged()+1);

            // -- Break if enough
            if (jdstate.get_num_converged()>=neigs)  break;

            // -- Remove eigenvector from space and update 
            jdstate.remove_converged_vector();

            DETECTED = 1;
            rKNOWN   = 0;

            // -- reset current eigenvalue counter since converged
            jtracker.set_citer(0);

        } else if (DETECTED)
            INITIATE = 1;

        // -- Restart if dim(V) >= jmax
        if (jdstate.get_size()==jdstate.get_max_size()) {
            //std::cout << "Restarting\n";
            jdstate.restart();
        }

    }

    // -- compute eigenvalues and eigenvectors
    JDQZ_Data* jdata = new JDQZ_Data(jdstate.get_num_converged(),v0);
    jdata->compute_evecs(jdstate.get_S(),jdstate.get_num_eigs(),
                         jdstate.get_T(),jdstate.get_num_eigs(),
                         jdstate.get_Q(),jdstate.get_Z());

    // -- Clean up
    VecDestroy(rt);
    VecDestroy(tt);

    return jdata;
}


JDQZ_Data* JDQZ(Mat A, Mat B, PC pc, Vec v0, 
          complex<double>* shift, int neigs, int maxiter, double atol, int initq, int minq, int maxq)
{
    JDQZ_ParameterList jpl;
    jpl.set_maxiter(maxiter);
    jpl.set_resid_atol(atol);
    jpl.set_initq(initq);
    jpl.set_minq(minq);
    jpl.set_maxq(maxq);
    JDQZ_Data* jdata = JDQZ(A,B,pc,v0,shift,neigs,&jpl);

    return jdata;
}


JDQZ_Data* JDQZ(Mat A, Mat B, PC pc, Vec v0, 
          double shift, int neigs, JDQZ_ParameterList* jpl)
{
    complex<double> shifts[2] = {complex<double>(shift), complex<double>(1)};

    JDQZ_Data* jdata =  JDQZ(A,B,pc,v0,shifts,neigs,jpl);
    return jdata;
}


JDQZ_Data* JDQZ(Mat A, Mat B, PC pc, Vec v0, 
          double shift, int neigs, int maxiter, double atol, int initq, int minq, int maxq)
{
    complex<double> shifts[2] = {complex<double>(shift), complex<double>(1)};

    JDQZ_Data* jdata =  JDQZ(A,B,pc,v0,shifts,neigs,maxiter,atol,initq,minq,maxq);
    return jdata;
}
