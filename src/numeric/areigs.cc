/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstdio>
#include <vector>
#include <algorithm>

#include "qcomplex.h"
#include "areigs.h"

using std::vector;
using std::copy;
using std::fill;


extern "C" {

int dsaupd_(int*        IDO,
            const char* BMAT,
            int*        N,
            const char* WHICH,
            int*        NEV,
            double*     TOL,
            double*     RESID,
            int*        NCV,
            double*     V,
            int*        LDV,
            int*        IPARAM,
            int*        IPNTR,
            double*     WORKD,
            double*     WORKL,
            int*        LWORKL,
            int*        INFO );

int dseupd_(int*        RVEC,
            const char* HOWMNY,
            int*        SELECT,
            double*     D,
            double*     Z,
            int*        LDZ,
            double*     SIGMA,
            const char* BMAT,
            int*        N,
            const char* WHICH,
            int*        NEV,
            double*     TOL,
            double*     RESID,
            int*        NCV,
            double*     V,
            int*        LDV,
            int*        IPARAM,
            int*        IPNTR,
            double*     WORKD,
            double*     WORKL,
            int*        LWORKL,
            int*        INFO );

int dnaupd_(int*        IDO,
            const char* BMAT,
            int*        N,
            const char* WHICH,
            int*        NEV,
            double*     TOL,
            double*     RESID,
            int*        NCV,
            double*     V,
            int*        LDV,
            int*        IPARAM,
            int*        IPNTR,
            double*     WORKD,
            double*     WORKL,
            int*        LWORKL,
            int*        INFO );

int dneupd_(int*        RVEC,
            const char* HOWMNY,
            int*        SELECT,
            double*     DR,
            double*     DI,
            double*     Z,
            int*        LDZ,
            double*     SIGMAR,
            double*     SIGMAI,
            double*     WORKEV,
            const char* BMAT,
            int*        N,
            const char* WHICH,
            int*        NEV,
            double*     TOL,
            double*     RESID,
            int*        NCV,
            double*     V,
            int*        LDV,
            int*        IPARAM,
            int*        IPNTR,
            double*     WORKD,
            double*     WORKL,
            int*        LWORKL,
            int*        INFO );

int znaupd_(int*        IDO,
            const char* BMAT,
            int*        N,
            const char* WHICH,
            int*        NEV,
            double*     TOL,
            dcomplex*   RESID,
            int*        NCV,
            dcomplex*   V,
            int*        LDV,
            int*        IPARAM,
            int*        IPNTR,
            dcomplex*   WORKD,
            dcomplex*   WORKL,
            int*        LWORKL,
            double*     RWORK,
            int*        INFO );

int zneupd_(int*        RVEC,
            const char* HOWMNY,
            int*        SELECT,
            dcomplex*   D,
            dcomplex*   Z,
            int*        LDZ,
            dcomplex*   SIGMA,
            dcomplex*   WORKEV,
            const char* BMAT,
            int*        N,
            const char* WHICH,
            int*        NEV,
            double*     TOL,
            dcomplex*   RESID,
            int*        NCV,
            dcomplex*   V,
            int*        LDV,
            int*        IPARAM,
            int*        IPNTR,
            dcomplex*   WORKD,
            dcomplex*   WORKL,
            int*        LWORKL,
            double*     RWORK,
            int*        INFO );

}


#include "arinfo.cc"


// ---- General settings ----


Arpack::Arpack()
{
    n = 0;
    nev = 6;
    ncv = 20;
    rvec = 1;
    mode = 1;
    maxitr = 30;

    bmat[0] = 'I';
    bmat[1] = '\0';

    which[0] = 'L';
    which[1] = 'M';
    which[2] = '\0';

    tol   = 0.0;
}


Arpack::~Arpack()
{
}


void Arpack::set_n(int n_)
{
    n = n_;
}


void Arpack::set_nev(int nev_)
{
    nev = nev_;
    ncv = 2*nev_;
    if (ncv < 20)
        ncv = 20;
}


void Arpack::set_ncv(int ncv_)
{
    ncv = ncv_;
    if (ncv < nev)
        ncv = nev;
}


void Arpack::set_rvec(int rvec_)
{
    rvec = rvec_;
}


void Arpack::set_maxitr(int maxitr_)
{
    maxitr = maxitr_;
}


void Arpack::set_which(char* which_)
{
    which[0] = which_[0];
    which[1] = which_[1];
}


void Arpack::set_tol(double tol_)
{
    tol = tol_;
}


// ---- DS settings ----


ArpackDS::ArpackDS()
{
    set_mode(1);
    set_shift(0);
}


ArpackDS::~ArpackDS()
{
}


void ArpackDS::set_mode(int mode_)
{
    mode = mode_;
    if (mode == 1) {
        bmat[0] = 'I';
    } else {
        bmat[0] = 'G';
    }
}


void ArpackDS::set_shift(double sigma_)
{
    sigma = sigma_;
}


int ArpackDS::compute_eigs(double* d, double* v)
{
    vector<double> resid(n);
    vector<double> workd(3*n);
    vector<double> workl((ncv+8)*ncv);
    vector<int>    select(ncv);

    int lworkl = (ncv+8)*ncv;
    int iparam[11];
    int ipntr[14];
    int info   = 0;
    int ierr   = 0;
    int status = 0;
    int iter   = 0;
    int ido    = 0;

    // -- Clear temporaries
    fill(iparam, iparam+11, 0);
    fill(ipntr,  ipntr+14,  0);

    // -- Set eigensolver heuristics
    iparam[0] = 1;       // Exact shift strategy
    iparam[2] = maxitr;  // Maximum number of iterations
    iparam[6] = mode;    // Which mode?

    // -- Reverse communication Arnoldi loop
    ido = 0;
    while (1) {

        dsaupd_(&ido, bmat, &n, which, &nev, &tol, &resid[0],
                &ncv, v, &n, iparam, ipntr, &workd[0], 
                &workl[0], &lworkl, &info);

        if (ido == -1) {
            times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else if (ido == 1) {
            if (mode == 3 || mode == 4 || mode == 5) {
                times_OP2(&workd[ipntr[0]-1], &workd[ipntr[1]-1],
                          &workd[ipntr[2]-1]);
            } else {
                times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
            }
        } else if (ido == 2) {
            times_M(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else {
            break;
        }

        ++iter;
    }

    if (info != 0) {

        // -- Check for error from Arnoldi
        printf("Error in dsaupd_ = %d\n", info);
        dsaupd_info_string(info);
        status = -1;

    } else {

        // -- Post processing
        dseupd_(&rvec, "A", &select[0], d, v, &n, &sigma,
                bmat, &n, which, &nev, &tol, &resid[0], &ncv,
                v, &n, iparam, ipntr, &workd[0], &workl[0], &lworkl, &ierr);

        if (ierr != 0) {
            printf("Error in dseupd_ = %d\n", ierr);
            dseupd_info_string(ierr);
            status = -2;
        }
    }

    return status;
}


void ArpackDS::times_OP1(double* x, double* opx)
{
}


void ArpackDS::times_OP2(double* x, double* opx, double* Mx)
{
    times_OP1(x, opx);
}


void ArpackDS::times_M(double* x, double* Mx)
{
    copy(Mx, Mx+n, x);
}


// ---- DN settings ----


ArpackDN::ArpackDN()
{
    set_mode(1);
    set_shift(0);
}


ArpackDN::~ArpackDN()
{
}


void ArpackDN::set_mode(int mode_)
{
    mode = mode_;
    if (mode == 1) {
        bmat[0] = 'I';
    } else if (mode == 2) {
        bmat[0] = 'G';
    } else if (mode == 3) {
        bmat[0] = 'G';
    } else if (mode == 4) {
        bmat[0] = 'G';
    }
}


void ArpackDN::set_shift(double sigmar_, double sigmai_)
{
    sigmar = sigmar_;
    sigmai = sigmai_;
}


int ArpackDN::compute_eigs(double* dr, double* di, double* v)
{
    vector<double> resid(n);
    vector<double> workd(3*n);
    vector<double> workl((3*ncv+6)*ncv);
    vector<double> workev(3*ncv);
    vector<int>    select(ncv);

    int lworkl = (3*ncv+6)*ncv;
    int iparam[11];
    int ipntr[14];
    int info   = 0;
    int ierr   = 0;
    int status = 0;
    int iter   = 0;
    int ido    = 0;

    // -- Clear temporaries
    fill(iparam, iparam+11, 0);
    fill(ipntr,  ipntr+14,  0);

    // -- Set eigensolver heuristics
    iparam[0] = 1;       // Exact shift strategy
    iparam[2] = maxitr;  // Maximum number of iterations
    iparam[6] = mode;    // Which mode?

    // -- Reverse communication Arnoldi loop
    ido = 0;
    while (1) {

        dnaupd_(&ido, bmat, &n, which, &nev, &tol, &resid[0],
                &ncv, v, &n, iparam, ipntr, &workd[0], &workl[0], &lworkl,
                &info);

        if (ido == -1) {
            times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else if (ido == 1) {
            if (mode == 3 || mode == 4)
                times_OP2(&workd[ipntr[0]-1], &workd[ipntr[1]-1],
                          &workd[ipntr[2]-1]);
            else
                times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else if (ido == 2) {
            times_M(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else {
            break;
        }

        ++iter;
    }

    if (info != 0) {

        // -- Check for error from Arnoldi
        printf("Error in dnaupd_ = %d\n", info);
        dnaupd_info_string(info);
        status = -1;

    } else {

        // -- Post processing
        dneupd_(&rvec, "A", &select[0], dr, di, v, &n, &sigmar, &sigmai,
                &workev[0], bmat, &n, which, &nev, &tol, &resid[0], &ncv,
                v, &n, iparam, ipntr, &workd[0], &workl[0], &lworkl, &ierr);

        if (ierr != 0) {
            printf("Error in dneupd_ = %d\n", ierr);
            dneupd_info_string(ierr);
            status = -2;
        }
    }

    return status;
}


void ArpackDN::times_OP1(double* x, double* opx)
{
}


void ArpackDN::times_OP2(double* x, double* opx, double* Mx)
{
    times_OP1(x, opx);
}


void ArpackDN::times_M(double* x, double* Mx)
{
    copy(Mx, Mx+n, x);
}


// ---- ZN settings ----


ArpackZN::ArpackZN()
{
    set_mode(1);
    set_shift(0);
}


ArpackZN::~ArpackZN()
{
}


void ArpackZN::set_mode(int mode_)
{
    mode = mode_;
    if (mode == 1) {
        bmat[0] = 'I';
    } else if (mode == 2) {
        bmat[0] = 'G';
    } else if (mode == 3) {
        bmat[0] = 'G';
    }
}


void ArpackZN::set_shift(dcomplex sigma_)
{
    sigma = sigma_;
}


int ArpackZN::compute_eigs(dcomplex* d, dcomplex* v)
{
    vector<dcomplex> resid(n);
    vector<dcomplex> workd(3*n);
    vector<dcomplex> workl((3*ncv+5)*ncv);
    vector<dcomplex> workev(3*ncv);
    vector<double>   rwork(ncv);
    vector<int>      select(ncv);

    int lworkl = (3*ncv+5)*ncv;
    int iparam[11];
    int ipntr[14];
    int info   = 0;
    int ierr   = 0;
    int status = 0;
    int iter   = 0;
    int ido    = 0;

    // -- Clear temporaries
    fill(iparam, iparam+11, 0);
    fill(ipntr,  ipntr+14,  0);

    // -- Set eigensolver heuristics
    iparam[0] = 1;       // Exact shift strategy
    iparam[2] = maxitr;  // Maximum number of iterations
    iparam[6] = mode;    // Which mode?

    // -- Reverse communication Arnoldi loop
    ido = 0;
    while (1) {

        znaupd_(&ido, bmat, &n, which, &nev, &tol, &resid[0],
                &ncv, v, &n, iparam, ipntr, &workd[0], &workl[0], &lworkl,
                &rwork[0], &info);

        if (ido == -1) {
            times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else if (ido == 1) {
            if (mode == 3)
                times_OP2(&workd[ipntr[0]-1], &workd[ipntr[1]-1],
                          &workd[ipntr[2]-1]);
            else
                times_OP1(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else if (ido == 2) {
            times_M(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
        } else {
            break;
        }

        ++iter;
    }

    if (info != 0) {

        // -- Check for error from Arnoldi
        printf("Error in znaupd_ = %d\n", info);
        znaupd_info_string(info);
        status = -1;

    } else {

        // -- Post processing
        zneupd_(&rvec, "A", &select[0], d, v, &n, &sigma, &workev[0], bmat, &n,
                which, &nev, &tol, &resid[0], &ncv, v, &n, iparam, ipntr,
                &workd[0], &workl[0], &lworkl, &rwork[0], &ierr);
        if (ierr != 0) {
            printf("Error in zneupd_ = %d\n", ierr);
            zneupd_info_string(ierr);
            status = -2;
        }
    }

    return status;
}


void ArpackZN::times_OP1(dcomplex* x, dcomplex* opx)
{
}


void ArpackZN::times_OP2(dcomplex* x, dcomplex* opx, dcomplex* Mx)
{
    times_OP1(x, opx);
}


void ArpackZN::times_M(dcomplex* x, dcomplex* Mx)
{
    copy(Mx, Mx+n, x);
}
