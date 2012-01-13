/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>

#include "qcomplex.h"
#include "pareigs.h"

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


// ---- DS settings ----


PArpackDS::PArpackDS()
{
    set_mode(1);
    set_shift(0);
}


PArpackDS::~PArpackDS()
{
}


void PArpackDS::set_mode(int mode_)
{
    mode = mode_;
    if (mode == 1) {
        bmat[0] = 'I';
    } else {
        bmat[0] = 'G';
    }
}


void PArpackDS::set_shift(double sigma_)
{
    sigma = sigma_;
}


int PArpackDS::compute_eigs(double* d, double* v)
{
    double* resid;
    double* workd;
    double* workl;

    if (mypid() == 0) {
        resid   = new double[n];
        workd   = new double[3*n];
        workl   = new double[(ncv+8)*ncv];
    }

    int  lworkl = (ncv+8)*ncv;
    int  iparam[11];
    int  ipntr[14];
    int  info = 0;
    int  ierr = 0;
    int  status = 0;
    int* select = new int[ncv];
    int  iter = 0;
    int  ido = 0;

    // -- Clear temporaries
    memset(iparam, 0, 11*sizeof(int));
    memset(ipntr,  0, 14*sizeof(int));

    // -- Set eigensolver heuristics
    iparam[0] = 1;       // Exact shift strategy
    iparam[2] = maxitr;  // Maximum number of iterations
    iparam[6] = mode;    // Which mode?

    // -- Reverse communication Arnoldi loop
    barrier(); // -- 1
    ido = 0;
    while (1) {

        barrier();
        // -- Calling dnaupd on root processor
        if (mypid() == 0) {

            dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
                    &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl,
                    &info);
        }

        // -- Broadcast parameter
        barrier(); // -- 2
        broadcast( &ido, 1, 0);
        broadcast(&info, 1, 0);

        barrier(); // -- 3
        if (ido == -1) {
            times_OP1(workd+ipntr[0]-1, workd+ipntr[1]-1);
        } else if (ido == 1) {
            if (mode == 3 || mode == 4 || mode == 5) {
                times_OP2(workd+ipntr[0]-1, workd+ipntr[1]-1,
                          workd+ipntr[2]-1);
            } else {
                times_OP1(workd+ipntr[0]-1, workd+ipntr[1]-1);
            }
        } else if (ido == 2) {
            times_M(workd+ipntr[0]-1, workd+ipntr[1]-1);
        } else {
            break;
        }

        ++iter;
    }

    barrier(); // -- 4
    if (info != 0) {

        // -- Check for error from Arnoldi
        printf("Error in dsaupd_ = %d\n", info);
        dsaupd_info_string(info);
        status = -1;

    } else {

        // -- Post processing
        // -- Calling dseupd on root processor
        if (mypid() == 0) {
            dseupd_(&rvec, "A", select, d, v, &n, &sigma,
                    bmat, &n, which, &nev, &tol, resid, &ncv,
                    v, &n, iparam, ipntr, workd, workl, &lworkl, &ierr);
        }

        // -- Broadcast parameter
        barrier(); // -- 5
        broadcast(&ierr, 1, 0);

        if (ierr != 0) {
            printf("Error in dseupd_ = %d\n", ierr);
            dseupd_info_string(ierr);
            status = -2;
        }
    }

    // -- Free allocated stuff
    barrier(); // --6
    if (mypid() == 0) {
        delete[] resid;
        delete[] workd;
        delete[] workl;
    }
    delete[] select;

    barrier();
    return status;
}

void PArpackDS::times_OP1(double* x, double* opx)
{
}

void PArpackDS::times_OP2(double* x, double* opx, double* Mx)
{
    times_OP1(x, opx);
}

void PArpackDS::times_M(double* x, double* Mx)
{
    memcpy(Mx, x, n*sizeof(double));
}

// ---- DN settings ----

PArpackDN::PArpackDN()
{
    set_mode(1);
    set_shift(0);
}


PArpackDN::~PArpackDN()
{
}


void PArpackDN::set_mode(int mode_)
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


void PArpackDN::set_shift(double sigmar_, double sigmai_)
{
    sigmar = sigmar_;
    sigmai = sigmai_;
}


int PArpackDN::compute_eigs(double* dr, double* di, double* v)
{
    double* resid;
    double* workd;
    double* workl;
    double* workev;

    if (mypid() == 0) {
        resid   = new double[n];
        workd   = new double[3*n];
        workl   = new double[(3*ncv+6)*ncv];
        workev  = new double[3*ncv];
    }

    int  lworkl = (3*ncv+6)*ncv;
    int  iparam[11];
    int  ipntr[14];
    int  info = 0;
    int  ierr = 0;
    int  status = 0;
    int* select = new int[ncv];
    int  iter = 0;
    int  ido = 0;

    // -- Clear temporaries
    memset(iparam, 0, 11*sizeof(int));
    memset(ipntr,  0, 14*sizeof(int));

    // -- Set eigensolver heuristics
    iparam[0] = 1;       // Exact shift strategy
    iparam[2] = maxitr;  // Maximum number of iterations
    iparam[6] = mode;    // Which mode?

    // -- Reverse communication Arnoldi loop
    barrier(); // -- 1
    ido = 0;
    while (1) {

        barrier();
        // -- Calling dnaupd on root processor
        if (mypid() == 0) {
            dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
                    &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl,
                    &info);
        }

        // -- Broadcast parameter
        barrier(); // -- 2
        broadcast( &ido, 1, 0);
        broadcast(&info, 1, 0);

        barrier(); // -- 3
        if (ido == -1) {
            times_OP1(workd+ipntr[0]-1, workd+ipntr[1]-1);
        } else if (ido == 1) {
            if (mode == 3 || mode == 4)
                times_OP2(workd+ipntr[0]-1, workd+ipntr[1]-1,
                          workd+ipntr[2]-1);
            else
                times_OP1(workd+ipntr[0]-1, workd+ipntr[1]-1);
        } else if (ido == 2) {
            times_M(workd+ipntr[0]-1, workd+ipntr[1]-1);
        } else {
            break;
        }

        ++iter;
    }

    barrier(); // -- 4
    if (info != 0) {

        // -- Check for error from Arnoldi
        printf("Error in dnaupd_ = %d\n", info);
        dnaupd_info_string(info);
        status = -1;

    } else {

        // -- Post processing
        // -- Calling dneupd on root processor
        if (mypid() == 0) {
            dneupd_(&rvec, "A", select, dr, di, v, &n, &sigmar, &sigmai,
                    workev, bmat, &n, which, &nev, &tol, resid, &ncv,
                    v, &n, iparam, ipntr, workd, workl, &lworkl, &ierr);
        }

        // -- Broadcast parameter
        barrier(); // -- 5
        broadcast(&ierr, 1, 0);

        if (ierr != 0) {
            printf("Error in dneupd_ = %d\n", ierr);
            dneupd_info_string(ierr);
            status = -2;
        }
    }

    // -- Free allocated stuff
    barrier(); // --6
    if (mypid() == 0) {
        delete[] resid;
        delete[] workd;
        delete[] workl;
        delete[] workev;
    }
    delete[] select;

    barrier();
    return status;
}

void PArpackDN::times_OP1(double* x, double* opx)
{
}

void PArpackDN::times_OP2(double* x, double* opx, double* Mx)
{
    times_OP1(x, opx);
}

void PArpackDN::times_M(double* x, double* Mx)
{
    memcpy(Mx, x, n*sizeof(double));
}

// ---- ZN settings ----

PArpackZN::PArpackZN()
{
    set_mode(1);
    set_shift(0);
}


PArpackZN::~PArpackZN()
{
}


void PArpackZN::set_mode(int mode_)
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


void PArpackZN::set_shift(dcomplex sigma_)
{
    sigma = sigma_;
}


int PArpackZN::compute_eigs(dcomplex* d, dcomplex* v)
{
    dcomplex* resid;
    dcomplex* workd;
    dcomplex* workl;
    dcomplex* workev;
    double*   rwork;

    if (mypid() == 0) {
        resid   = new dcomplex[n];
        workd   = new dcomplex[3*n];
        workl   = new dcomplex[(3*ncv+5)*ncv];
        workev  = new dcomplex[3*ncv];
        rwork   = new double  [ncv];
    }

    int  lworkl = (3*ncv+5)*ncv;
    int  iparam[11];
    int  ipntr[14];
    int  info = 0;
    int  ierr = 0;
    int  status = 0;
    int* select = new int[ncv];
    int  iter = 0;
    int  ido = 0;

    // -- Clear temporaries
    memset(iparam, 0, 11*sizeof(int));
    memset(ipntr,  0, 14*sizeof(int));

    // -- Set eigensolver heuristics
    iparam[0] = 1;       // Exact shift strategy
    iparam[2] = maxitr;  // Maximum number of iterations
    iparam[6] = mode;    // Which mode?

    // -- Reverse communication Arnoldi loop
    barrier(); // -- 1
    ido = 0;
    while (1) {

        barrier();
        // -- Calling znaupd on root processor
        if (mypid() == 0) {
std::cout << "Starting znaupd:" << iter << "\n";
            znaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
                    &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl,
                    rwork, &info);
std::cout << "Ending znaupd:" << iter << "\n";
        }

        // -- Broadcast parameter
        barrier(); // -- 2
        broadcast( &ido, 1, 0);
        broadcast(&info, 1, 0);

        barrier(); // -- 3
        if (ido == -1) {
            times_OP1(workd+ipntr[0]-1, workd+ipntr[1]-1);
        } else if (ido == 1) {
            if (mode == 3)
                times_OP2(workd+ipntr[0]-1, workd+ipntr[1]-1,
                          workd+ipntr[2]-1);
            else
                times_OP1(workd+ipntr[0]-1, workd+ipntr[1]-1);
        } else if (ido == 2) {
            times_M(workd+ipntr[0]-1, workd+ipntr[1]-1);
        } else {
            break;
        }

        ++iter;
    }

    barrier(); // -- 4
    if (info != 0) {

        // -- Check for error from Arnoldi
        printf("Error in znaupd_ = %d\n", info);
        znaupd_info_string(info);
        status = -1;

    } else {

        // -- Post processing
        // -- Calling dneupd on root processor
        if (mypid() == 0) {
            zneupd_(&rvec, "A", select, d, v, &n, &sigma, workev, bmat, &n,
                    which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr,
                    workd, workl, &lworkl, rwork, &ierr);
        }

        // -- Broadcast parameter
        barrier(); // -- 5
        broadcast(&ierr, 1, 0);

        if (ierr != 0) {
            printf("Error in zneupd_ = %d\n", ierr);
            zneupd_info_string(ierr);
            status = -2;
        }
    }

    // -- Free allocated stuff
    barrier(); // --6
    if (mypid() == 0) {
        delete[] resid;
        delete[] workd;
        delete[] workl;
        delete[] workev;
        delete[] rwork;
    }
    delete[] select;

    barrier();
    return status;
}

void PArpackZN::times_OP1(dcomplex* x, dcomplex* opx)
{
}

void PArpackZN::times_OP2(dcomplex* x, dcomplex* opx, dcomplex* Mx)
{
    times_OP1(x, opx);
}

void PArpackZN::times_M(dcomplex* x, dcomplex* Mx)
{
    memcpy(Mx, x, n*sizeof(dcomplex));
}
