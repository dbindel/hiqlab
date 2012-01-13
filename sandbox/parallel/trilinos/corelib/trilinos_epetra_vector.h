#ifndef _TRILINOS_EPETRA_VECTOR_H
#define _TRILINOS_EPETRA_VECTOR_H

/** @file trilinos_epetra_vector.h
 *
 *  Classes and functions related to Epetra
 */

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

#include "qcomplex.h"

/** Complex data structure for Epetra_MultiVector
 */
class Epetra_MultiVector_Complex : public Epetra_MultiVector {
 public:

    /** Create a new vector */
    Epetra_MultiVector_Complex(const Epetra_BlockMap& Map, int NumVectors,
                               int is_real=0, bool zeroOut = true);

    /** Copy constructors */
    Epetra_MultiVector_Complex(const Epetra_MultiVector_Complex& Source);
    Epetra_MultiVector_Complex(const Epetra_MultiVector& Vx,
                               const Epetra_MultiVector& Vz);

    Epetra_MultiVector_Complex(Epetra_DataAccess CV, const Epetra_BlockMap& Map,
                       double *Vx, double *Vz,
                       int MyLDA, int NumVectors);
    Epetra_MultiVector_Complex(const Epetra_BlockMap& Map,
                       dcomplex *V,
                       int MyLDA, int NumVectors);

    Epetra_MultiVector_Complex(Epetra_DataAccess CV,
                       const Epetra_MultiVector_Complex& Source,
                       int *Indices, int NumVectors);
    Epetra_MultiVector_Complex(Epetra_DataAccess CV,
                       const Epetra_MultiVector& Vx, const Epetra_MultiVector& Vz,
                       int *Indices, int NumVectors);

    Epetra_MultiVector_Complex(Epetra_DataAccess CV,
                       const Epetra_MultiVector_Complex& Source,
                       int StartIndex, int NumVectors);
    Epetra_MultiVector_Complex(Epetra_DataAccess CV,
                       const Epetra_MultiVector& Vx, const Epetra_MultiVector& Vz,
                       int StartIndex,int NumVectors);

    virtual ~Epetra_MultiVector_Complex();

    /** View or get imag part of vector */
    const Epetra_MultiVector& view_Vz() const { return *Vz; };
    Epetra_MultiVector* get_Vz()              { return Vz;  };

    /** Is real vector ? */
    int is_real_vector() const { return is_real; };

    /** Fill with random numbers */
    int Random();

    /** Extract copy */
    int ExtractCopy(double* vr, double* vi, int MyLDA);
    int ExtractCopy(dcomplex* v, int MyLDA);

    /** Scale vector */
    int Scale(double sr, double si);
    int Scale(dcomplex s);
    int Scale(double ar, double ai, const Epetra_MultiVector_Complex& A);
    int Scale(dcomplex a, const Epetra_MultiVector_Complex& A);

    /** Update vector */
    int Update(double ar, double ai, const Epetra_MultiVector_Complex&A,
               double sr, double si);
    int Update(dcomplex a,const Epetra_MultiVector_Complex&A, dcomplex s);
    int Update(double ar, double ai, const Epetra_MultiVector_Complex&A,
               double br, double bi, const Epetra_MultiVector_Complex&B,
               double sr, double si);
    int Update(dcomplex a, const Epetra_MultiVector_Complex&A,
               dcomplex b, const Epetra_MultiVector_Complex&B,
               dcomplex s);

    /** Set scalar value */
    int PutScalar(double ar, double ai);
    int PutScalar(dcomplex a);

    /** Compute dot product */
    int Dot(const Epetra_MultiVector_Complex& A, double* dr, double* di) const;
    int Dot(const Epetra_MultiVector_Complex& A, dcomplex* d) const;

    /** Compute 2-Norm of vectors */
    int Norm2(double* dr, double* di) const;
    int Norm2(dcomplex* d) const;
    int Norm2(double* d) const;

    /** Multiply vectors */
    int Multiply (char TransA, char TransB, double sabr, double sabi,
                  const Epetra_MultiVector_Complex& A,
                  const Epetra_MultiVector_Complex& B,
                  double sr, double si );
    int Multiply (char TransA, char TransB, dcomplex sab,
                  const Epetra_MultiVector_Complex& A,
                  const Epetra_MultiVector_Complex& B,
                  dcomplex s);
    int Multiply (double sabr, double sabi,
                  const Epetra_MultiVector_Complex& A,
                  const Epetra_MultiVector_Complex& B,
                  double sr, double si );
    int Multiply (dcomplex sab,
                  const Epetra_MultiVector_Complex& A,
                  const Epetra_MultiVector_Complex& B,
                  dcomplex s);

 private:
    int is_real;
    Epetra_MultiVector* Vz;   /* Imag part     */

};


/** Complex data structure for Epetra_Vector
 */
class Epetra_Vector_Complex : public Epetra_Vector {
 public:

    /** Create a new vector */
    Epetra_Vector_Complex(const Epetra_BlockMap& Map, int is_real=0,
                                             bool zeroOut = true);

    /** Copy constructors */
    Epetra_Vector_Complex(const Epetra_Vector_Complex& Source);
    Epetra_Vector_Complex(const Epetra_Vector& Vx, const Epetra_Vector& Vz);

    Epetra_Vector_Complex(Epetra_DataAccess CV, const Epetra_BlockMap& Map,
                          double *Vx, double *Vz);
    Epetra_Vector_Complex(const Epetra_BlockMap& Map, dcomplex *V);

    Epetra_Vector_Complex(Epetra_DataAccess CV,
                    const Epetra_MultiVector_Complex& Source, int Index);
    Epetra_Vector_Complex(Epetra_DataAccess CV, const Epetra_MultiVector& Vx,
                    const Epetra_MultiVector& Vz, int Index);

    virtual ~Epetra_Vector_Complex();

    /** View or get imag part of vector */
    const Epetra_Vector& view_Vz() const { return *Vz; };
    Epetra_Vector* get_Vz()              { return Vz;  };

    /** Is real matrix ? */
    int is_real_vector() const { return is_real; };

    /** Fill with random numbers */
    int Random();

    /** Scale vector */
    int Scale(double sr, double si);
    int Scale(dcomplex s);
    int Scale(double ar, double ai, Epetra_Vector_Complex& A);
    int Scale(dcomplex a, Epetra_Vector_Complex& A);

    /** Update vector */
    int Update(double ar, double ai, Epetra_Vector_Complex&A,
               double sr, double si);
    int Update(dcomplex a, Epetra_Vector_Complex&A, dcomplex s);
    int Update(double ar, double ai, Epetra_Vector_Complex&A,
               double br, double bi, Epetra_Vector_Complex&B,
               double sr, double si);
    int Update(dcomplex a, Epetra_Vector_Complex&A,
               dcomplex b, Epetra_Vector_Complex&B,
               dcomplex s);

    /** Set scalar value */
    int PutScalar(double ar, double ai);
    int PutScalar(dcomplex a);

    /** Compute dot product */
    int Dot(const Epetra_Vector_Complex& A, double* dr, double* di) const;
    int Dot(const Epetra_Vector_Complex& A, dcomplex* d) const;

    /** Compute 2-Norm of vectors */
    int Norm2(double* dr, double* di);
    int Norm2(dcomplex* d);
    int Norm2(double* dr);

 private:
    int is_real;
    Epetra_Vector* Vz;   /* Imag part     */

};

#endif /* _TRILINOS_EPETRA_VECTOR_H */
