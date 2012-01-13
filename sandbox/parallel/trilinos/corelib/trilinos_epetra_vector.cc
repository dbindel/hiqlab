/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */
#include <iostream>

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "trilinos_epetra_vector.h"

// -- Epetra_MultiVector_Complex
// -- Create new object
Epetra_MultiVector_Complex::Epetra_MultiVector_Complex(const Epetra_BlockMap& Map,
                                       int NumVectors, int is_real, bool zeroOut)
 : Epetra_MultiVector(Map, NumVectors, zeroOut), is_real(is_real), Vz(NULL)
{
    if (!is_real)
        Vz = new Epetra_MultiVector(Map, NumVectors, zeroOut);
}

// -- Copy constructors
Epetra_MultiVector_Complex::Epetra_MultiVector_Complex(
                                const Epetra_MultiVector_Complex& Source)
 : Epetra_MultiVector(Source), is_real(Source.is_real_vector()), Vz(NULL)
{
    if (!is_real)
        Vz = new Epetra_MultiVector( Source.view_Vz() );
}

Epetra_MultiVector_Complex::Epetra_MultiVector_Complex(
                const Epetra_MultiVector& xVx, const Epetra_MultiVector& xVz)
 : Epetra_MultiVector(xVx), is_real(0), Vz(NULL)
{
        Vz = new Epetra_MultiVector(xVz);
}

Epetra_MultiVector_Complex::Epetra_MultiVector_Complex(Epetra_DataAccess CV,
                const Epetra_BlockMap& Map,
                double *xVx, double *xVz, int MyLDA, int NumVectors)
 : Epetra_MultiVector(CV,Map,xVx,MyLDA,NumVectors), is_real(0), Vz(NULL)
{
    // FIXME: THIS IS NOT CORRECT(IMPLEMENT METHOD BELOW)
    Vz = new Epetra_MultiVector(CV, Map, xVz, MyLDA, NumVectors);
}

Epetra_MultiVector_Complex::Epetra_MultiVector_Complex(const Epetra_BlockMap& Map,
                dcomplex *V, int MyLDA, int NumVectors)
 : Epetra_MultiVector(Map, NumVectors, true), is_real(0), Vz(NULL)
{
    double* Vr;
    double* Vi;
    dcomplex* VV;
//    double* Vr = new double[MyLDA*NumVectors];
//    double* Vi = new double[MyLDA*NumVectors];
    int lda;

    Vz = new Epetra_MultiVector(Map, NumVectors, true);

    this->ExtractView(&Vr,&lda);
    Vz->ExtractView(&Vi,&lda);

    for (int i = 0; i < NumVectors; ++i) {
        for (int j = 0; j < lda; ++j) {
            Vr[j+lda*i] = real(V[j+MyLDA*i]);
            Vi[j+lda*i] = imag(V[j+MyLDA*i]);
        }
    }
    // This will not work when NumRow neq LDA
    //copy_complex(V, Vr, Vi, MyLDA*NumVectors);
}

Epetra_MultiVector_Complex::Epetra_MultiVector_Complex(Epetra_DataAccess CV,
      const Epetra_MultiVector_Complex& Source, int *Indices, int NumVectors)
 : Epetra_MultiVector(CV, Source, Indices, NumVectors), is_real(0), Vz(NULL)
{
    Vz = new Epetra_MultiVector(CV, Source.view_Vz(), Indices, NumVectors);
}

Epetra_MultiVector_Complex::Epetra_MultiVector_Complex(Epetra_DataAccess CV,
      const Epetra_MultiVector& xVx, const Epetra_MultiVector& xVz,
      int *Indices, int NumVectors)
 : Epetra_MultiVector(CV, xVx, Indices, NumVectors), is_real(0), Vz(NULL)
{
    Vz = new Epetra_MultiVector(CV, xVz, Indices, NumVectors);
}

Epetra_MultiVector_Complex::Epetra_MultiVector_Complex(Epetra_DataAccess CV,
      const Epetra_MultiVector_Complex& Source,
      int StartIndex,int NumVectors)
 : Epetra_MultiVector(CV, Source, StartIndex, NumVectors), is_real(0), Vz(NULL)
{
    Vz = new Epetra_MultiVector(CV, Source.view_Vz(), StartIndex, NumVectors);
}

Epetra_MultiVector_Complex::Epetra_MultiVector_Complex(Epetra_DataAccess CV,
      const Epetra_MultiVector& xVx, const Epetra_MultiVector& xVz,
      int StartIndex,int NumVectors)
 : Epetra_MultiVector(CV, xVx, StartIndex, NumVectors), is_real(0), Vz(NULL)
{
    Vz = new Epetra_MultiVector(CV, xVz, StartIndex, NumVectors);
}

Epetra_MultiVector_Complex::~Epetra_MultiVector_Complex()
{
    if (Vz)
         delete Vz;
}

int Epetra_MultiVector_Complex::Random()
{
    int ierr = 0;

    ierr += Epetra_MultiVector::Random();
    if (Vz)
        ierr += Vz->Random();

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::ExtractCopy(double* vr, double* vi, int MyLDA)
{
    int ierr = 0;

    ierr += Epetra_MultiVector::ExtractCopy(vr, MyLDA);
    ierr += Vz->ExtractCopy(vi, MyLDA);

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::ExtractCopy(dcomplex* v, int MyLDA)
{
    double* vr;
    double* vi;

    int n = Epetra_MultiVector::NumVectors();
    int ierr = 0;
    vr = new double[n*MyLDA];
    vi = new double[n*MyLDA];
    ierr += ExtractCopy(vr, vi, MyLDA);

    copy_complex(vr, vi, v, n*MyLDA);

    delete[] vr;
    delete[] vi;

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Scale(double sr, double si)
{
    int ierr = 0;

    Epetra_MultiVector Sr = Epetra_MultiVector(*this);
    Epetra_MultiVector Si = Epetra_MultiVector(*Vz);

    // -- Real part
    ierr += Epetra_MultiVector::Update(-si, Si, sr);

    // -- Imag part
    ierr += Vz->Update( si, Sr,  sr);

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Scale(dcomplex s)
{
    return Scale(real(s), imag(s));
}

int Epetra_MultiVector_Complex::Scale(double ar, double ai,
                                      const Epetra_MultiVector_Complex& A)
{
    int ierr = 0;

    // -- Real part
    ierr += Epetra_MultiVector::Update(-ai, A.view_Vz(), ar, A, 0.0);

    // -- Imag part
    ierr += Vz->Update(ai, A, ar, A.view_Vz(), 0.0);

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Scale(dcomplex a, const Epetra_MultiVector_Complex& A)
{
    return Scale(real(a), imag(a), A);
}

int Epetra_MultiVector_Complex::Update(double ar, double ai,
                                       const Epetra_MultiVector_Complex&A,
                                       double sr, double si)
{
    int ierr = 0;
    // FIXME: Copying is expensive. Referece EMV::Update
    Epetra_MultiVector Sr = Epetra_MultiVector(*this);
    Epetra_MultiVector Si = Epetra_MultiVector(*Vz);
    Epetra_MultiVector Ar = Epetra_MultiVector(A);
    Epetra_MultiVector Ai = Epetra_MultiVector(A.view_Vz());

    // -- Real part
    ierr += Epetra_MultiVector::Update( ar, Ar,  sr);
    ierr += Epetra_MultiVector::Update(-ai, Ai, -si, Si, 1.0);

    // -- Imag part
    ierr += Vz->Update( ar, Ai,  sr);
    ierr += Vz->Update( ai, Ar,  si, Sr, 1.0);

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Update(dcomplex a, const Epetra_MultiVector_Complex&A,
                                       dcomplex s)
{
    return Update(real(a), imag(a), A, real(s), imag(s));
}

int Epetra_MultiVector_Complex::Update(double ar, double ai,
   const Epetra_MultiVector_Complex&A, double br, double bi,
   const Epetra_MultiVector_Complex&B, double sr, double si)
{
    int ierr = 0;
    // FIXME: Copying is expensive. Referece EMV::Update
    Epetra_MultiVector_Complex Az = Epetra_MultiVector_Complex(A);
    Epetra_MultiVector_Complex Bz = Epetra_MultiVector_Complex(B);

    ierr += Update(ar, ai, Az, sr, si);
    ierr += Update(br, bi, Bz, 1.0, 0.0);

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Update(dcomplex a, const Epetra_MultiVector_Complex&A,
                                       dcomplex b, const Epetra_MultiVector_Complex&B,
                                       dcomplex s)
{
    return Update(real(a), imag(a), A, real(b), imag(b), B, real(s), imag(s));
}

int Epetra_MultiVector_Complex::PutScalar(double ar, double ai)
{
    int ierr = 0;

    ierr += Epetra_MultiVector::PutScalar(ar);
    ierr += Vz->PutScalar(ai);

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::PutScalar(dcomplex a)
{
    return PutScalar(real(a), imag(a));
}

int Epetra_MultiVector_Complex::Dot(const Epetra_MultiVector_Complex& A,
                                    double* dr, double* di) const
{
    int ierr = 0;

    int NumVectors = A.NumVectors();
    double* tempr = new double[NumVectors];
    double* tempi = new double[NumVectors];

    // -- Real part
    ierr += Epetra_MultiVector::Dot(A, dr);
    ierr += Vz->Dot(A.view_Vz(), tempr);

    // -- Imag part
    ierr += Epetra_MultiVector::Dot(A.view_Vz(), tempi);
    ierr += Vz->Dot(A, di);
    //ierr += Epetra_MultiVector::Dot(A.view_Vz(), di);
    //ierr += Vz->Dot(A, tempi);

    for (int i = 0; i < NumVectors; ++i) {
        dr[i]+=tempr[i];
        di[i]-=tempi[i];
    }

    delete[] tempr;
    delete[] tempi;

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Dot(const Epetra_MultiVector_Complex& A,
                                    dcomplex* d) const
{
    int ierr = 0;

    int NumVectors = A.NumVectors();
    double* tempr = new double[NumVectors];
    double* tempi = new double[NumVectors];

    ierr += Dot(A, tempr, tempi);
    copy_complex(tempr, tempi, d, NumVectors);

    delete[] tempr;
    delete[] tempi;

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Norm2(double* dr, double* di) const
{
    int ierr = 0;

    int NumVectors = Epetra_MultiVector::NumVectors();

    ierr += Dot(*this, dr, di);
    for (int i = 0; i < NumVectors; ++i) {
        dr[i]=sqrt(dr[i]);
        di[i]=0.0;
    }

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Norm2(double* dr) const
{
    int ierr = 0;

    int NumVectors = Epetra_MultiVector::NumVectors();
    double* tempi = new double[NumVectors];

    ierr += Dot(*this, dr, tempi);
    for (int i = 0; i < NumVectors; ++i) {
        dr[i]=sqrt(dr[i]);
    }

    delete[] tempi;

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Norm2(dcomplex* d) const
{
    int ierr = 0;

    int NumVectors = Epetra_MultiVector::NumVectors();
    double* tempr = new double[NumVectors];
    double* tempi = new double[NumVectors];

    ierr += Norm2(tempr, tempi);
    copy_complex(tempr, tempi, d, NumVectors);

    delete[] tempr;
    delete[] tempi;

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Multiply (char TransA, char TransB,
                  double sabr, double sabi,
                  const Epetra_MultiVector_Complex& A,
                  const Epetra_MultiVector_Complex& B,
                  double sr, double si )
{
    int ierr = 0;
    // FIXME: Copying is expensive. Referece EMV::Update
    Epetra_MultiVector Ar = Epetra_MultiVector(A);
    Epetra_MultiVector Ai = Epetra_MultiVector(A.view_Vz());
    Epetra_MultiVector Br = Epetra_MultiVector(B);
    Epetra_MultiVector Bi = Epetra_MultiVector(B.view_Vz());

    // Must take conjugate transpose when 'T'
    if (TransA=='T')
        Ai.Scale(-1.0);
    if (TransB=='T')
        Bi.Scale(-1.0);
    // -- Scale first
    ierr += Scale(sr, si);
    // -- real part
    ierr += Epetra_MultiVector::Multiply(TransA, TransB,
                  sabr, Ar         , Br         , 1.0);
    ierr += Epetra_MultiVector::Multiply(TransA, TransB,
                 -sabr, Ai         , Bi         , 1.0);
    ierr += Epetra_MultiVector::Multiply(TransA, TransB,
                 -sabi, Ar         , Bi         , 1.0);
    ierr += Epetra_MultiVector::Multiply(TransA, TransB,
                 -sabi, Ai         , Br         , 1.0);

    // -- imag part
    ierr += Vz->Multiply(TransA, TransB,
                  sabi, Ar         , Br         , 1.0);
    ierr += Vz->Multiply(TransA, TransB,
                 -sabi, Ai         , Bi         , 1.0);
    ierr += Vz->Multiply(TransA, TransB,
                  sabr, Ar         , Bi         , 1.0);
    ierr += Vz->Multiply(TransA, TransB,
                  sabr, Ai         , Br         , 1.0);

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Multiply (char TransA, char TransB,
                  dcomplex sab,
                  const Epetra_MultiVector_Complex& A,
                  const Epetra_MultiVector_Complex& B,
                  dcomplex s)
{
    return Multiply(TransA, TransB, real(sab), imag(sab), A, B,
                                    real(s)  , imag(s));
}

int Epetra_MultiVector_Complex::Multiply (double sabr, double sabi,
                  const Epetra_MultiVector_Complex& A,
                  const Epetra_MultiVector_Complex& B,
                  double sr, double si )
{
    int ierr = 0;

    // FIXME: Copying is expensive. Referece EMV::Update
    Epetra_MultiVector Ar = Epetra_MultiVector(A);
    Epetra_MultiVector Ai = Epetra_MultiVector(A.view_Vz());
    Epetra_MultiVector Br = Epetra_MultiVector(B);
    Epetra_MultiVector Bi = Epetra_MultiVector(B.view_Vz());

    // -- Scale first
    ierr += Scale(sr, si);

    // -- real part
    ierr += Epetra_MultiVector::Multiply(
                  sabr, Ar         , Br         , 1.0);
    ierr += Epetra_MultiVector::Multiply(
                 -sabr, Ai         , Bi         , 1.0);
    ierr += Epetra_MultiVector::Multiply(
                 -sabi, Ar         , Bi         , 1.0);
    ierr += Epetra_MultiVector::Multiply(
                 -sabi, Ai         , Br         , 1.0);

    // -- imag part
    ierr += Vz->Multiply(
                  sabi, Ar         , Br         , 1.0);
    ierr += Vz->Multiply(
                 -sabi, Ai         , Bi         , 1.0);
    ierr += Vz->Multiply(
                  sabr, Ar         , Bi         , 1.0);
    ierr += Vz->Multiply(
                  sabr, Ai         , Br         , 1.0);

    return (ierr==0)?0:1;
}

int Epetra_MultiVector_Complex::Multiply (dcomplex sab,
                  const Epetra_MultiVector_Complex& A,
                  const Epetra_MultiVector_Complex& B,
                  dcomplex s)
{
    return Multiply(real(sab), imag(sab), A, B, real(s), imag(s));
}

// -- Epetra_Vector_Complex
// -- Create new object
Epetra_Vector_Complex::Epetra_Vector_Complex(const Epetra_BlockMap& Map,
                                             int is_real, bool zeroOut)
 : Epetra_Vector(Map, zeroOut), is_real(is_real), Vz(NULL)
{
    if (!is_real)
        Vz = new Epetra_Vector(Map, zeroOut);
}

// -- Copy constructors
Epetra_Vector_Complex::Epetra_Vector_Complex(const Epetra_Vector_Complex& Source)
 : Epetra_Vector(Source), is_real(Source.is_real_vector()), Vz(NULL)
{
    if (!is_real)
        Vz = new Epetra_Vector( Source.view_Vz() );
}

Epetra_Vector_Complex::Epetra_Vector_Complex(const Epetra_Vector& xVx,
                                             const Epetra_Vector& xVz)
 : Epetra_Vector(xVx), is_real(0), Vz(NULL)
{
    Vz = new Epetra_Vector(xVz);
}

Epetra_Vector_Complex::Epetra_Vector_Complex(Epetra_DataAccess CV,
                                      const Epetra_BlockMap& Map,
                                      double *xVx, double *xVz)
 : Epetra_Vector(CV, Map, xVx), is_real(0), Vz(NULL)
{
    Vz = new Epetra_Vector(CV, Map, xVz);
}

Epetra_Vector_Complex::Epetra_Vector_Complex(const Epetra_BlockMap& Map,
                                      dcomplex *V)
 : Epetra_Vector(Map, true), is_real(0), Vz(NULL)
{
    double* Vr;
    double* Vi;
    //= new double[Map.NumMyElements()];
    //= new double[Map.NumMyElements()];

    Vz = new Epetra_Vector(Map, true);

    this->ExtractView(&Vr);
    Vz->ExtractView(&Vi);

    copy_complex(V, Vr, Vi, Map.NumMyElements());
}

Epetra_Vector_Complex::Epetra_Vector_Complex(Epetra_DataAccess CV,
                              const Epetra_MultiVector_Complex& Source, int Index)
 : Epetra_Vector(CV, Source, Index), is_real(0), Vz(NULL)
{
    Vz = new Epetra_Vector(CV, Source.view_Vz(), Index);
}

Epetra_Vector_Complex::Epetra_Vector_Complex(Epetra_DataAccess CV,
                              const Epetra_MultiVector& xVx,
                              const Epetra_MultiVector& xVz, int Index)
 : Epetra_Vector(CV, xVx, Index), is_real(0), Vz(NULL)
{
    Vz = new Epetra_Vector(CV, xVz, Index);
}

Epetra_Vector_Complex::~Epetra_Vector_Complex()
{
    if (Vz)
         delete Vz;
}

int Epetra_Vector_Complex::Random()
{
    int ierr = 0;

    ierr += Epetra_Vector::Random();
    if (Vz)
        ierr += Vz->Random();

    return (ierr==0)?0:1;
}

int Epetra_Vector_Complex::Scale(double sr, double si)
{
    int ierr = 0;

    Epetra_Vector Sr = Epetra_Vector(*this);
    Epetra_Vector Si = Epetra_Vector(*Vz);

    // -- Real part
    ierr += Epetra_Vector::Update(-si, Si, sr);

    // -- Imag part
    ierr += Vz->Update( si, Sr,  sr);

    return (ierr==0)?0:1;
}

int Epetra_Vector_Complex::Scale(dcomplex s)
{
    return Scale(real(s), imag(s));
}

int Epetra_Vector_Complex::Scale(double ar, double ai, Epetra_Vector_Complex& A)
{
    int ierr = 0;

    // -- Real part
    ierr += Epetra_Vector::Update(-ai, *(A.get_Vz()), ar, A, 0.0);

    // -- Imag part
    ierr += Vz->Update(ai, A, ar, *(A.get_Vz()), 0.0);

    return (ierr==0)?0:1;
}

int Epetra_Vector_Complex::Scale(dcomplex a, Epetra_Vector_Complex& A)
{
    return Scale(real(a), imag(a), A);
}

int Epetra_Vector_Complex::Update(double ar, double ai, Epetra_Vector_Complex&A,
                                  double sr, double si)
{
    int ierr = 0;

    Epetra_Vector Sr = Epetra_Vector(*this);
    Epetra_Vector Si = Epetra_Vector(*Vz);

    // -- Real part
    ierr += Epetra_Vector::Update( ar,             A,  sr);
    ierr += Epetra_Vector::Update(-ai, *(A.get_Vz()), -si, Si, 1.0);

    // -- Imag part
    ierr += Vz->Update( ar, *(A.get_Vz()),  sr);
    ierr += Vz->Update( ai,             A,  si, Sr, 1.0);

    return (ierr==0)?0:1;
}

int Epetra_Vector_Complex::Update(dcomplex a, Epetra_Vector_Complex&A,
                                       dcomplex s)
{
    return Update(real(a), imag(a), A, real(s), imag(s));
}

int Epetra_Vector_Complex::Update(double ar, double ai, Epetra_Vector_Complex&A,
                                  double br, double bi, Epetra_Vector_Complex&B,
                                  double sr, double si)
{
    int ierr = 0;

    ierr += Update(ar, ai, A, sr, si);
    ierr += Update(br, bi, B, 1.0, 0.0);

    return (ierr==0)?0:1;
}

int Epetra_Vector_Complex::Update(dcomplex a, Epetra_Vector_Complex&A,
                                  dcomplex b, Epetra_Vector_Complex&B,
                                  dcomplex s)
{
    return Update(real(a), imag(a), A, real(b), imag(b), B, real(s), imag(s));
}

int Epetra_Vector_Complex::PutScalar(double ar, double ai)
{
    int ierr = 0;

    ierr += Epetra_Vector::PutScalar(ar);
    ierr += Vz->PutScalar(ai);

    return (ierr==0)?0:1;
}

int Epetra_Vector_Complex::PutScalar(dcomplex a)
{
    return PutScalar(real(a), imag(a));
}

int Epetra_Vector_Complex::Dot(const Epetra_Vector_Complex& A,
                               double* dr, double* di) const
{
    int ierr = 0;

    int NumVectors = 1;
    double* tempr = new double[NumVectors];
    double* tempi = new double[NumVectors];

    // -- Real part
    ierr += Epetra_Vector::Dot(A, dr);
    ierr += Vz->Dot(A.view_Vz(), tempr);

    // -- Imag part
    ierr += Epetra_Vector::Dot(A.view_Vz(), di);
    ierr += Vz->Dot(A, tempi);

    for (int i = 0; i < NumVectors; ++i) {
        dr[i]+=tempr[i];
        di[i]-=tempi[i];
    }

    delete[] tempr;
    delete[] tempi;

    return (ierr==0)?0:1;
}

int Epetra_Vector_Complex::Dot(const Epetra_Vector_Complex& A, dcomplex* d) const
{
    int ierr = 0;

    int NumVectors = 1;
    double* tempr = new double[NumVectors];
    double* tempi = new double[NumVectors];

    ierr += Dot(A, tempr, tempi);
    copy_complex(tempr, tempi, d, NumVectors);

    delete[] tempr;
    delete[] tempi;

    return (ierr==0)?0:1;
}

int Epetra_Vector_Complex::Norm2(double* dr, double* di)
{
    int ierr = 0;

    int NumVectors = 1;
    ierr += Dot(*this, dr, di);
    for (int i = 0; i < NumVectors; ++i) {
        dr[i]=sqrt(dr[i]);
        di[i]=0.0;
    }

    return (ierr==0)?0:1;
}

int Epetra_Vector_Complex::Norm2(double* dr)
{
    int ierr = 0;

    int NumVectors = 1;
    double* tempi = new double[NumVectors];

    ierr += Dot(*this, dr, tempi);
    for (int i = 0; i < NumVectors; ++i) {
        dr[i]=sqrt(dr[i]);
    }

    delete[] tempi;

    return (ierr==0)?0:1;
}

int Epetra_Vector_Complex::Norm2(dcomplex* d)
{
    int ierr = 0;

    int NumVectors = 1;
    double* tempr = new double[NumVectors];
    double* tempi = new double[NumVectors];

    ierr += Norm2(tempr, tempi);
    copy_complex(tempr, tempi, d, NumVectors);

    delete[] tempr;
    delete[] tempi;

    return (ierr==0)?0:1;
}
