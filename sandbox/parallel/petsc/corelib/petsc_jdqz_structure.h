/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PETSC_JDQZ_STRUCTURE_H
#define PETSC_JDQZ_STRUCTURE_H

#include <complex>
#include <string>
#include "petscksp.h"
#include "matshellab.h"

using std::complex;
using std::string;

/*@t ------------
 * \section{JDQZ related objects}
 * 
 * This file includes declarations for the classes that are required to
 * run the Jacobi-Davidson QZ algorithm. The components are
 *
 * - JDQZ_State
 * - JDQZ_Solve
 * - JDQZ_Tracker
 * - JDQZ_Data
 * - JDQZ_ParameterList
 *
 * The purpose of each class is described in the following subsections.
 *
 *@q*/


/*@t ----------
 * \subsection{JDQZ_State}
 *
 * This class manages and owns the vectors necessary in the JDQZ iteration.
 * It also includes the methods that act on these vectors to compute the 
 * approximate eigenvalues and expand the subspace with a given new direction.
 *
 *@c*/


class JDQZ_State {

  public:
    JDQZ_State(Mat A, Mat B, Vec v, int neigs, int minq, int maxq);
    virtual ~JDQZ_State();

    /** Create a size [initq] subspace with a starting vector v0 with the complex shift 
     *  pair (alpha,beta), where [alpha/beta] is the eigenvalue.
     */
    void initialize_arnoldi(int initq, complex<double> alpha, complex<double> beta, Vec v0);

    /* project matrices onto subspace
     */
    void compute_projections();

    /* Compute approximate eigenpair
     * The eigenpairs are sorted according to distance from (alpha,beta)
     */
    void compute_projected_eigenvalues(complex<double> alpha, complex<double> beta);

    /* Expand subspace
     * Construct harmonic vector from pair (alpha,beta)
     */
    void expand_subspace(Vec tt, complex<double> alpha, complex<double> beta);

    /* Compute residual
     */
    void compute_residual_vector(Vec rt);

    /* Expand partial Schur form
     */
    void expand_partial_schur();

    /* Remove converged vector
     */
    void remove_converged_vector();

    /* Conduct restart of vectors when sizeq==maxq
     */
    void restart();

    void set_size(int n){sizeq = n;}

    Vec  get_Q  (int i) {return Q[i];  }
    Vec* get_Q  ()      {return Q;     }
    Vec  get_Z  (int i) {return Z[i];  }
    Vec* get_Z  ()      {return Z;     }
    Vec  get_KiZ(int i) {return KiZ[i];}
    Vec* get_KiZ()      {return KiZ;   }
    Vec  get_V  (int i) {return V[i];  }
    Vec  get_W  (int i) {return W[i];  }
    int  get_size()     {return sizeq; }
    int  get_max_size() {return maxq; }
    int  get_num_converged(){return ceigs; }
    int  get_num_eigs()     {return neigs; }

    complex<double>* get_MA() {return MA;}
    complex<double>* get_MB() {return MB;}
    complex<double>* get_SA() {return SA;}
    complex<double>* get_TB() {return TB;}
    complex<double>* get_UR() {return UR;}
    complex<double>* get_UL() {return UL;}
    complex<double>* get_S()  {return S ;}
    complex<double>* get_T()  {return T ;}

    complex<double>  get_SA(int i, int j) {return SA[maxq*j+i];}
    complex<double>  get_TB(int i, int j) {return TB[maxq*j+i];}

  private:

    // -- Size parameters
    int  minq;
    int  maxq;
    int  neigs;

    // -- Counting parameters
    int  ceigs;
    int  sizeq;

    int  orthotype;

    Mat  A;
    Mat  B;
    Vec* Q;
    Vec* Z;
    Vec* KiZ;
    Vec* V;
    Vec* W;
    Vec  tempv;
    Vec* Vt;

    complex<double>* MA;
    complex<double>* MB;
    complex<double>* SA;
    complex<double>* TB;
    complex<double>* UR;
    complex<double>* UL;
    complex<double>* S;
    complex<double>* T;
};


/*@t ------------
 * \section{JDQZ_Solve}
 *
 * Class to hold the JDQZ solving method
 *@c*/


class JDQZ_Solve {

  public:
    JDQZ_Solve(Mat A, Mat B, PC pc);
    virtual ~JDQZ_Solve();

    void apply(Vec x, Vec y);
    void set_computed_size(int n);
    void set_current_size(int n);
    void set_shift(complex<double> alpha, complex<double> beta);
    void set_rtol(double rtol);
    void set_tolerances(double rtol, int maxits);
    void set_projection_spaces(Vec* Q, Vec* Z, Vec* KiZ, int n);
    void set_ksp_type(string ksptype);

    int  get_computed_size();
    int  get_current_size();

  private:
    MatShellAB* AsB;
    Mat AsBmat;
    KSP ksp;
    PC  ppc;
    PC  pc;
    Mat A,B;

    Vec* Q;
    Vec* Z;
    Vec* KiZ;

    int neigs;
    complex<double> alpha;
    complex<double> beta;
};


/*@t ------------
 * \section{JDQZ_Tracker}
 *
 * Class to hold the JDQZ tracking
 *@c*/


class JDQZ_Tracker {

  public:
    JDQZ_Tracker(int maxiter) : titer(0), citer(0), maxiter(maxiter), 
                     resid_atol(1e-4), eigdist_rtol(1e-1),
                     ksp_rtol(1e-10), ksp_maxits(60) {}
    virtual ~JDQZ_Tracker() {}

    void set_citer(int i)         {citer = i;}
    void set_niter(int i)         {titer = i;}
    void set_maxiter(int i)       {maxiter = i;}
    void set_resid_atol  (double tol) {resid_atol   = tol;}
    void set_eigdist_rtol(double tol) {eigdist_rtol = tol;}
    void set_ksp_rtol    (double tol) {ksp_rtol     = tol;}
    void set_ksp_maxits  (int i)      {ksp_maxits   = i;}
    void set_target(complex<double>* target, complex<double>* aim, 
                    complex<double>* approx,double nr);
    void set_ksp_tolerance(JDQZ_Solve* jsolve);
    void increment() {titer++;citer++;}

    int niter() {return titer;}
    int nciter(){return citer;}

  private:
    int    titer;
    int    citer;
    int    maxiter;
    double resid_atol;
    double eigdist_rtol;
    double ksp_rtol;
    int    ksp_maxits;
};


/*@t ------------
 * \section{JDQZ_Data}
 *
 * Class to hold the JDQZ data obtained from computation.
 *@c*/


class JDQZ_Data {

  public:
    JDQZ_Data(int n, Vec v0);
    virtual ~JDQZ_Data();

    void compute_evecs(double* SAr, double* SAi, int ldSA,
                       double* TBr, double* TBi, int ldTB,
                       Vec* Q, Vec* Z);

    void compute_evecs(complex<double>* SA, int ldSA, complex<double>* TB, int ldTB,
                       Vec* Q, Vec* Z);

    int  get_n() {return nvalues;}

    Vec  get_revec(int i) {return revecs[i];}
    Vec  get_levec(int i) {return levecs[i];}

    complex<double> get_alpha (int i){return           alpha[i];}
    double          get_alphar(int i){return std::real(alpha[i]);}
    double          get_alphai(int i){return std::imag(alpha[i]);}
    complex<double> get_beta  (int i){return           beta[i];}
    double          get_betar (int i){return std::real(beta[i]);}
    double          get_betai (int i){return std::imag(beta[i]);}

    complex<double> get_eval  (int i){return           alpha[i]/beta[i];}
    double          get_evalr (int i){return std::real(alpha[i]/beta[i]);}
    double          get_evali (int i){return std::imag(alpha[i]/beta[i]);}

  private:
    int  nvalues;
    Vec* revecs;
    Vec* levecs;
    complex<double>* alpha;
    complex<double>* beta;
};


/*@t ------------
 * \section{JDQZ_ParameterList}
 *
 * Class to hold the JDQZ parameter list.
 *@c*/


class JDQZ_ParameterList {

  public:
    JDQZ_ParameterList() :
      maxiter(20), initq(1), minq(50), maxq(100), 
      resid_atol(1e-10), tresid_atol(1e-4), teigdist_rtol(1e-1),
      ksp_rtol(1e-9), ksp_maxits(60), ksp_type("gmres") {}
    virtual ~JDQZ_ParameterList() {}

    void set_maxiter(int v) {maxiter=v;}
    void set_initq  (int v) {initq  =v;}
    void set_minq   (int v) {minq   =v;}
    void set_maxq   (int v) {maxq   =v;}

    void set_resid_atol   (double v) {resid_atol   =v;}
    void set_tresid_atol  (double v) {tresid_atol  =v;}
    void set_teigdist_rtol(double v) {teigdist_rtol=v;}

    void set_ksp_rtol    (double v) {ksp_rtol    =v;}
    void set_ksp_maxits  (int    v) {ksp_maxits  =v;}
    void set_ksp_type    (string v) {ksp_type    =v;}

    void check_parameters() {assert(minq<=maxq);assert(initq<=maxq);}

    // -- variables
    int          maxiter;
    int          initq;
    int          minq;
    int          maxq;

    double       resid_atol;
    double       tresid_atol;
    double       teigdist_rtol;

    double       ksp_rtol;
    int          ksp_maxits;
    string       ksp_type;
};

#endif /* PETSC_JDQZ_STRUCTURE_H */
