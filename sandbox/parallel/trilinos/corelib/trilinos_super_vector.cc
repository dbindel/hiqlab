/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include <algorithm>
#include <iostream>

#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Epetra_MpiComm.h"

#include "trilinos_super_vector.h"

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

Epetra_Super_MultiVector::Epetra_Super_MultiVector(Epetra_Map& RowMap, int NumVectors) :
    nsv(0), FullMap(&RowMap)
{
    C = new Epetra_MultiVector(RowMap,NumVectors);
}

Epetra_Super_MultiVector::~Epetra_Super_MultiVector()
{
  delete C;
}

void Epetra_Super_MultiVector::InsertVector(Epetra_MultiVector* A, IndexMap* im)
{
    nsv++;
    sub_vec.push_back(A);
    indexmap.push_back(im);
}

void Epetra_Super_MultiVector::ConstructVector()
{
     FillVector();
}

void Epetra_Super_MultiVector::FillVector()
{
    for (int ism = 0; ism < nsv; ++ism) {

        Epetra_MultiVector* A;

        // -- Construct target map of subvector
        Epetra_Map TargetMap = s2f_rowmap(ism);

        // -- Construct vector with new map
        A = new Epetra_MultiVector(TargetMap, sub_vec[ism]->NumVectors() );
        Epetra_Export Exporter( sub_vec[ism]->Map(), TargetMap);
        A->Export( *(sub_vec[ism]), Exporter, Add);

        // -- Insert nnz structure into supermatrix
        int GID[A->Map().NumMyElements()];
        for (int i = 0; i < A->Map().NumMyElements(); ++i) {
            GID[i] = indexmap[ism]->row(A->Map().GID(i));
        }

        // -- Fill matrix
        double* vals;
        int MyLDA;
        A->ExtractView(&vals,&MyLDA);
        for (int j = 0; j < A->NumVectors(); ++j) {
            for (int i = 0; i < A->Map().NumMyElements(); ++i) {
                C->SumIntoGlobalValue(GID[i],j,vals[i]);
            }
            vals+=MyLDA;
        }

        delete A;
    }
}


Epetra_Map Epetra_Super_MultiVector::s2f_rowmap(int ism)
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

void ExtractMultiVector(Epetra_MultiVector* Source, Epetra_MultiVector* B, IndexMap* im)
{
    std::vector<int> ExtractLID;
    double* svals;
    int* MappedID;

    int NumMyElements = 0;

    // -- Count number of elements in this process to be extracted
    for (int i = 0; i < Source->Map().NumMyElements(); ++i)
        if (im->row(Source->Map().GID(i)) >= 0) {
           ExtractLID.push_back(i);
           NumMyElements++;
        }
    svals    = new double[NumMyElements*Source->NumVectors()];
    MappedID = new int   [NumMyElements];

    // -- Extract values and put in designated array
    double* fvals;
    int MyLDA;
    Source->ExtractView(&fvals,&MyLDA);

    for (int i = 0; i < NumMyElements; ++i)
        MappedID[i] = im->row(Source->Map().GID(ExtractLID[i]));
    for (int j = 0; j < Source->NumVectors(); ++j) {
        for (int i = 0; i < NumMyElements; ++i)
            svals[i+j*NumMyElements] = fvals[ExtractLID[i]];
        fvals+=MyLDA;
    }

    // -- Construct map for the extracted vector
    Epetra_Map ExtMap(-1,NumMyElements, MappedID, 0, Source->Map().Comm());
    Epetra_MultiVector ExtVec(Copy,ExtMap,svals,NumMyElements,Source->NumVectors());

    // -- Export into B
    Epetra_Export Exporter( ExtMap, B->Map());
    B->Export( ExtVec, Exporter, Insert);

    // -- Cleanup
    delete[] svals;
    delete[] MappedID;

}
