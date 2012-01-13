#ifndef PETSC_JDQZ_UTIL_H
#define PETSC_JDQZ_UTIL_H

#include "petscksp.h"
#include "matshellab.h"

using std::complex;

void ProjectedGS(Vec* V, int size_v, Vec* Z, int size_z,
                Vec* Q, PetscScalar* R, int ldR, Vec tempv,
                int type=0, int gamma=1, int maxiter=3);

void InitArnoldi(Vec* V, Vec v0, Mat A, Mat B, PetscScalar* sigma, int minq);

void ProjectedOP(Mat A, Vec* V, int size_v, Vec* W, int size_w, Vec tempv,
                PetscScalar* MA, int ldMA, int imin=0, int jmin=0);

void HarmonicVector(Mat A, Mat B, PetscScalar* sigma, Vec v, Vec w, Vec Av);

/** Utitility functions implementing LAPACK
 *
 */
void ScaleEig(complex<double>* a, complex<double>* b, int n);
void ScaleEig(double* ar, double* ai, double* br, double* bi, int n);

void SortEig(complex<double>* s, complex<double>* t,
             complex<double>* a, complex<double>* b,
             int* index, int n, int sortCD=1);
void SortEig(double* sr, double* si, double* tr, double* ti,
             double ar, double ai, double br, double bi,
             int* index, int n, int sortCD=1);
void SortQZ(int n, complex<double>* Q, int ldQ,
                   complex<double>* Z, int ldZ,
                   complex<double>* S, int ldS,
                   complex<double>* T, int ldT,
                   complex<double>* MA,int ldMA,
                   complex<double>* MB,int ldMB,
                   complex<double> alpha, complex<double> beta, 
                   int gamma=1);
void SwapQZ(int n, complex<double>* Q, int ldQ,
                   complex<double>* Z, int ldZ,
                   complex<double>* S, int ldS,
                   complex<double>* T, int ldT,
                   int* index, int nsort);

#endif
