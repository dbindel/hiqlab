#ifndef PETSC_JDQZ_H
#define PETSC_JDQZ_H

#include <complex>
#include "petscksp.h"

using std::complex;

class JDQZ_Data;
class JDQZ_ParameterList;

JDQZ_Data* JDQZ(Mat A, Mat B, PC pc, Vec v0, 
          complex<double>* shift, int neigs, JDQZ_ParameterList* jpl);
JDQZ_Data* JDQZ(Mat A, Mat B, PC pc, Vec v0, 
          complex<double>* shift, int neigs, int maxiter, double atol, int initq, int minq, int maxq);
JDQZ_Data* JDQZ(Mat A, Mat B, PC pc, Vec v0, 
          double shift, int neigs, JDQZ_ParameterList* jpl);
JDQZ_Data* JDQZ(Mat A, Mat B, PC pc, Vec v0, 
          double shift, int neigs, int maxiter, double atol, int initq, int minq, int maxq);

#endif /* PETSC_JDQZ_H */
