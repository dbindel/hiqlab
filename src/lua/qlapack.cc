/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include "qlapack.h"
#include <cassert>
#include <cstring>
#include <cstdio>
#include <cmath>

extern "C"
int dgemm_(const char& transA,
           const char& transB,
           const int& M, const int& N, const int& K,
           const double& alpha,
           const double* A, const int& ldA,
           const double* B, const int& ldB,
           const double& beta,
           double* C, const int& ldC);


void dgemm(char* trans, double alpha, QArray* A, QArray* B,
           double beta, QArray* C)
{
    assert(A->type() == 0 &&
           B->type() == 0 &&
           C->type() == 0);

    char transA = 'N';
    char transB = 'N';
    if (trans) {
        if (trans[0])
            transA = trans[0];
        if (trans[0] && trans[1])
            transB = trans[1];
    }

    assert(transA == 'N' || transA == 'T');
    assert(transB == 'N' || transB == 'T');

    int M = C->m();
    int N = C->n();
    int K;

    if (transA == 'N') {
        assert(A->m() == M);
        K = A->n();
    } else {
        assert(A->n() == M);
        K = A->m();
    }

    if (transB == 'N') {
        assert(B->m() == K);
        assert(B->n() == N);
    } else {
        assert(B->m() == N);
        assert(B->n() == K);
    }

    dgemm_(transA, transB, M, N, K, alpha,
           A->data_r(), A->lda(),
           B->data_r(), B->lda(),
           beta,
           C->data_r(), C->lda());
}


extern "C"
int dgetrf_(const int& M, const int& N, double* A,
            const int& ldA, int* ipiv, int& info);

int dgetrf(QArray* A, QIArray* ipiv)
{
    assert(A->type() == 0);
    assert(ipiv->m() >= A->m() &&
           ipiv->m() >= A->n());

    int info = 0;
    dgetrf_(A->m(), A->n(), A->data_r(), A->lda(),
            ipiv->data(), info);
    return info;
}


extern "C"
int dgetrs_(const char& trans, const int& N, const int& NRHS,
            double* A, const int& ldA, int* ipiv,
            double* B, const int& ldB, int& info);

int dgetrs(char* trans, QArray* A, QIArray* ipiv, QArray* B)
{
    assert(A->type() == 0);
    assert(A->m() == A->n());
    assert(ipiv->m() >= A->m());
    assert(B->m() == A->m());

    char t = (trans ? *trans : 'N');
    assert(t == 'N' || t == 'T');

    int info = 0;
    dgetrs_(t, A->m(), B->n(), A->data_r(), A->lda(), ipiv->data(),
            B->data_r(), B->lda(), info);
    return info;
}


extern "C"
int dgeev_(const char& jobVL, const char& jobVR,
           const int& N, double* A, const int& ldA,
           double* wr, double* wi,
           double* VL, const int& ldVL,
           double* VR, const int& ldVR,
           double* work, const int& lwork, int& info);

int dgeev(QArray* A, QArray* w, QArray* vl, QArray* vr)
{
    char jobVL = 'N';
    char jobVR = 'N';
    double* vlx = NULL;
    double* vrx = NULL;
    double* wrx = NULL;
    double* wix = NULL;

    int ldVR = 1;
    int ldVL = 1;

    assert(A->type() == 0);
    assert(A->m() == A->n());
    assert(w->type() != 1);

    if (vl) {
        jobVL = 'V';
        vlx = vl->data_r();
        ldVL = vl->lda();
        assert(vl->type() == 0);
        assert(vl->m() >= A->m());
        assert(vl->n() >= A->m());
    }

    if (vr) {
        jobVR = 'V';
        vrx = vr->data_r();
        ldVR = vr->lda();
        assert(vr->type() == 0);
        assert(vr->m() >= A->m());
        assert(vr->n() >= A->m());
    }

    if (w->type() == 0) {
        assert(w->m() >= A->m());
        assert(w->n() == 2);
        wrx = w->data_r();
        wix = w->data_r() + A->lda();
    } else if (w->type() == 2) {
        assert(w->m() >= A->m());
        assert(w->n() == 1);
        wrx = w->data_r();
        wix = w->data_i();
    }

    int info = 0;
    int lwork = 0;
    double* work = NULL;
    double  work1 = 0;

    // Workspace query call
    dgeev_(jobVL, jobVR, A->m(), A->data_r(), A->lda(),
           wrx, wix, vlx, ldVL, vrx, ldVR, &work1, -1, info);

    // Call with workspace
    lwork = (int) work1;
    work  = new double[lwork];
    dgeev_(jobVL, jobVR, A->m(), A->data_r(), A->lda(),
           wrx, wix, vlx, ldVL, vrx, ldVR, work, lwork, info);
    delete[] work;

    return info;
}
