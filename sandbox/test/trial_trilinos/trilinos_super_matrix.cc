/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 * $Id: trilinos_super_matrix.cc,v 1.2 2006/08/08 13:07:20 tkoyama Exp $
 */

#include <algorithm>
#include <iostream>

#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Epetra_MpiComm.h"

#include "trilinos_super_matrix.h"

struct map_pair {
    int l;
    int g;
};

static int compare_map_pair(const map_pair& pair1, const map_pair& pair2)
{
    if (pair1.g  < pair2.g ) return  1;
    if (pair1.g  > pair2.g ) return  0;
    if (pair1.l  < pair2.l ) return  1;
    if (pair1.l  > pair2.l ) return  0;
    return 0;
}

Epetra_Super_CrsMatrix::Epetra_Super_CrsMatrix(Epetra_Map& RowMap) : 
    nsm(0), FullMap(&RowMap)
{
}

Epetra_Super_CrsMatrix::~Epetra_Super_CrsMatrix()
{
  delete C;
}

void Epetra_Super_CrsMatrix::InsertMatrix(Epetra_CrsMatrix* A, IndexMap* im)
{
    nsm++;
    sub_mat.push_back(A);
    indexmap.push_back(im); 
}

void Epetra_Super_CrsMatrix::ConstructMatrix()
{
     CreateMatrix();
     FillMatrix();
}

void Epetra_Super_CrsMatrix::CreateMatrix()
{
    int MyNNZ[FullMap->NumMyElements()];
    memset(MyNNZ, 0, FullMap->NumMyElements() * sizeof(int)); 
    setMyNNZ(MyNNZ);
    C = new Epetra_CrsMatrix(Copy,*FullMap,MyNNZ); 
}

void Epetra_Super_CrsMatrix::FillMatrix()
{
    for (int ism = 0; ism < nsm; ++ism) {
        
        Epetra_CrsMatrix* A;

        // -- Construct target map of submatrix
        Epetra_Map TargetMap = s2f_rowmap(ism);

        // -- Obtain nnz of the submatrix in the new map
        int MyNNZ_i[TargetMap.NumMyElements()];
        NumMyEntriesPerRow( *(sub_mat[ism]), TargetMap, MyNNZ_i);

        // -- Construct matrix with new map
        A = new Epetra_CrsMatrix(Copy, TargetMap, MyNNZ_i);
        Epetra_Export Exporter( sub_mat[ism]->RowMap(), TargetMap);
        A->Export( *(sub_mat[ism]), Exporter, Add);
        A->FillComplete(); 

        // -- Insert nnz structure into supermatrix
        int GID[A->RowMap().NumMyElements()];
        for (int i = 0; i < A->RowMap().NumMyElements(); ++i)
            GID[i] = indexmap[ism]->row(A->RowMap().GID(i));

        // -- Fill matrix
        for (int i = 0; i < A->RowMap().NumMyElements(); ++i) {
            double val_t[A->NumMyEntries(i)];
            int   id_t [A->NumMyEntries(i)];
            int   nnz_t;
            A->ExtractGlobalRowCopy(A->RowMap().GID(i),A->NumMyEntries(i),nnz_t,val_t,id_t);
            for (int j = 0; j < nnz_t; ++j)
               id_t[j] = indexmap[ism]->col(id_t[j]);

            // -- Filter values
            double val_f[A->NumMyEntries(i)];
            int   id_f [A->NumMyEntries(i)];
            int   nnz_f=0;
            for (int j = 0; j < nnz_t; ++j) {
                if (id_t[j] >= 0) {
                    id_f[nnz_f] = id_t[j];
                    val_f[nnz_f]= val_t[j];
                    nnz_f++;
                }
            }
            C->InsertGlobalValues(GID[i],nnz_f,val_f,id_f);
        }
        delete A;
    }
    C->FillComplete();
}


void Epetra_Super_CrsMatrix::setMyNNZ(int* MyNNZ)
{
    for (int ism = 0; ism < nsm; ++ism) {

        Epetra_CrsMatrix* A;

        // -- Construct target map of submatrix
        Epetra_Map TargetMap = s2f_rowmap(ism);

        // -- Obtain nnz of the submatrix in the new map
        int MyNNZ_i[TargetMap.NumMyElements()];
        NumMyEntriesPerRow( *(sub_mat[ism]), TargetMap, MyNNZ_i);

        // -- Construct matrix with new map
        A = new Epetra_CrsMatrix(Copy, TargetMap, MyNNZ_i);
        Epetra_Export Exporter( sub_mat[ism]->RowMap(), TargetMap);
        A->Export( *(sub_mat[ism]), Exporter, Add);
        A->FillComplete(); 

        // -- Insert nnz structure into supermatrix
        int GID[A->RowMap().NumMyElements()];
        for (int i = 0; i < A->RowMap().NumMyElements(); ++i) 
            GID[i] = indexmap[ism]->row(A->RowMap().GID(i));

        // -- Fill matrix
        for (int i = 0; i < A->RowMap().NumMyElements(); ++i) {
            double val_t[A->NumMyEntries(i)];
            int   id_t [A->NumMyEntries(i)];
            int   nnz_t;
            A->ExtractGlobalRowCopy(A->RowMap().GID(i),A->NumMyEntries(i),nnz_t,val_t,id_t);
            for (int j = 0; j < nnz_t; ++j)
               id_t[j] = indexmap[ism]->col(id_t[j]);

            // -- Filter values
            double val_f[A->NumMyEntries(i)];
            int   id_f [A->NumMyEntries(i)];
            int   nnz_f=0;
            for (int j = 0; j < nnz_t; ++j) {
                if (id_t[j] >= 0) {
                    id_f[nnz_f] = id_t[j];
                    val_f[nnz_f]= val_t[j];
                    nnz_f++;
                }
            }
//            C->InsertGlobalValues(GID[i],nnz_f,val_f,id_f);
            MyNNZ[FullMap->LID(GID[i])] += nnz_f;
        }
/*
        // -- Set NNZ structure
        for (int i = 0; i < A->RowMap().NumMyElements(); ++i)
            MyNNZ[FullMap->LID(GID[i])] += A->NumMyEntries(i);
*/
        delete A;
    }
}


void Epetra_Super_CrsMatrix::NumMyEntriesPerRow(const Epetra_CrsMatrix& A, const Epetra_Map& TargetMap, 
                      int* MyNNZ_B)
{
    Epetra_Export Exporter(A.RowMap(),TargetMap);

    // Compute the nonzero structure of current distributed matrix
    Epetra_Vector id_A(A.RowMap());
    Epetra_Vector id_B(TargetMap);
    int NumMyElements_A = A.RowMap().NumMyElements();
    double MyNNZ_A[NumMyElements_A];
    int MyID_A[NumMyElements_A];
    for (int i = 0; i <= NumMyElements_A; ++i) {
        MyNNZ_A[i] = (double)(A.NumMyEntries(i));
        MyID_A[i] = i;
    }
    id_A.ReplaceMyValues(NumMyElements_A,MyNNZ_A,MyID_A);
    id_B.Export(id_A,Exporter,Add);

    // Compute the nonzero structure of new distributed matrix
    int NumMyElements_B = TargetMap.NumMyElements();
    double id_B_d[NumMyElements_B];
    id_B.ExtractCopy(id_B_d); 
    for (int i = 0; i <= NumMyElements_B; ++i) {
        MyNNZ_B[i] = (int)(id_B_d[i]);
    }

} 

Epetra_Map Epetra_Super_CrsMatrix::s2f_rowmap(int ism)
{
    // Select row mapping function

    // Find ids within this processor
    std::vector<map_pair> id_pair;
    map_pair id_pair_t;
    int numid = -indexmap[ism]->row(-1);
    int NumMyElements = 0;
    for (int i = 0; i < numid; ++i) {
        int gid = indexmap[ism]->row(i);
        if (FullMap->MyGID(gid)) {
            NumMyElements++;
            id_pair_t.l = i;
            id_pair_t.g = gid;
            id_pair.push_back(id_pair_t);
        }
    }

    // Sort (localmatrixid, globalmatrixid) pair
    sort(id_pair.begin(),id_pair.end(),compare_map_pair);

    int MappedID[NumMyElements];
    for (int i = 0; i < NumMyElements; ++i)
        MappedID[i] = id_pair[i].l;

    return Epetra_Map(-1,NumMyElements, MappedID, 0, FullMap->Comm());
}
