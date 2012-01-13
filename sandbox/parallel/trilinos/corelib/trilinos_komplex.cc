/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include <iostream>

#include "trilinos_komplex.h"

#include "trilinos_super_matrix.h"
#include "trilinos_super_vector.h"
#include "trilinos_indexmap.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

class IndexMap0 : public IndexMap {
 public:
   IndexMap0(int sop):sop(sop) {};
   ~IndexMap0() {};
   int inline row(int id)
   {
       if (0<=id && id<=sop) return 2*id;
       return -sop;
   };
   int inline col(int id)
   {
       if (0<=id && id<=sop) return 2*id;
       return -sop;
   }
 private:
   int sop;
};

class IndexMap1 : public IndexMap {
 public:
   IndexMap1(int sop):sop(sop) {};
   ~IndexMap1() {};
   int inline row(int id)
   {
       if (0<=id && id<=sop) return 2*id;
       return -sop;
   };
   int inline col(int id)
   {
       if (0<=id && id<=sop) return 2*id+1;
       return -sop;
   }
 private:
   int sop;
};

class IndexMap2 : public IndexMap {
 public:
   IndexMap2(int sop):sop(sop) {};
   ~IndexMap2() {};
   int inline row(int id)
   {
       if (0<=id && id<=sop) return 2*id+1;
       return -sop;
   };
   int inline col(int id)
   {
       if (0<=id && id<=sop) return 2*id;
       return -sop;
   }
 private:
   int sop;
};

class IndexMap3 : public IndexMap {
 public:
   IndexMap3(int sop):sop(sop) {};
   ~IndexMap3() {};
   int inline row(int id)
   {
       if (0<=id && id<=sop) return 2*id+1;
       return -sop;
   };
   int inline col(int id)
   {
       if (0<=id && id<=sop) return 2*id+1;
       return -sop;
   }
 private:
   int sop;
};

class IndexMapE0 : public IndexMap{
 public:
   IndexMapE0(int sop):sop(sop) {};
   ~IndexMapE0() {};
   int inline row(int id)
   {
       if (id%2 == 0) return id/2;
       return -sop;
   };
   int inline col(int id)
   {
       if (id%2 == 0) return id/2;
       return -sop;
   }
 private:
   int sop;
};

class IndexMapE1 : public IndexMap{
 public:
   IndexMapE1(int sop):sop(sop) {};
   ~IndexMapE1() {};
   int inline row(int id)
   {
       if (id%2 == 1) return (id-1)/2;
       return -sop;
   };
   int inline col(int id)
   {
       if (id%2 == 1) return (id-1)/2;
       return -sop;
   }
 private:
   int sop;
};

Epetra_Super_CrsMatrix* Complex2SuperCrsMatrix(Epetra_CrsMatrix_Complex* A,
                                               int form)
{
    Epetra_Super_CrsMatrix* SM;
    Epetra_CrsMatrix* Atemp;

    // -- Construct FullMap
    Epetra_Map FullMap(A->NumGlobalRows()*2,0, A->Comm());

    // Construct super matrix
    SM = new Epetra_Super_CrsMatrix(FullMap);

    if (form==0) {

        IndexMap0 im0(A->NumGlobalRows());
        IndexMap1 im1(A->NumGlobalRows());
        IndexMap2 im2(A->NumGlobalRows());
        IndexMap3 im3(A->NumGlobalRows());

        SM->InsertMatrix(A          ,&im0);
        Atemp = new Epetra_CrsMatrix(*(A->get_Az()));
        Atemp->Scale(-1.0);
        SM->InsertMatrix(Atemp      ,&im1);
        SM->InsertMatrix(A->get_Az(),&im2);
        SM->InsertMatrix(A          ,&im3);

    }
    SM->ConstructMatrix();

    // -- Clean up
    delete Atemp;

    return SM;
}

Epetra_Super_MultiVector* Complex2SuperMultiVector(Epetra_MultiVector_Complex* A,                                               int form)
{
    Epetra_Super_MultiVector* SM;

    // -- Construct FullMap
    Epetra_Map FullMap(A->GlobalLength()*2,0, A->Comm());

    // Construct super multivector
    SM = new Epetra_Super_MultiVector(FullMap,A->NumVectors());

    if (form==0) {
        IndexMap0 im0(A->GlobalLength());
        IndexMap3 im3(A->GlobalLength());

        SM->InsertVector(A,          &im0);
        SM->InsertVector(A->get_Vz(),&im3);
    }
    SM->ConstructVector();

    return SM;
}

void Multi2ComplexMultiVector(Epetra_MultiVector* Source,
                              Epetra_MultiVector_Complex* A,
                                               int form)
{
    if (form==0) {

        IndexMapE0 im0(A->GlobalLength());
        IndexMapE1 im1(A->GlobalLength());

        ExtractMultiVector(Source, A,&im0);
        ExtractMultiVector(Source, A->get_Vz(),&im1);

    }
}

Epetra_Super_MultiVector* Complex2SuperMultiVector(Epetra_Vector_Complex* A,                                               int form)
{
    Epetra_Super_MultiVector* SM;

    // -- Construct FullMap
    Epetra_Map FullMap(A->GlobalLength()*2,0, A->Comm());

    // Construct super multivector
    SM = new Epetra_Super_MultiVector(FullMap,1);

    if (form==0) {
        IndexMap0 im0(A->GlobalLength());
        IndexMap3 im3(A->GlobalLength());

        SM->InsertVector(A,          &im0);
        SM->InsertVector(A->get_Vz(),&im3);
    }
    SM->ConstructVector();

    return SM;
}
