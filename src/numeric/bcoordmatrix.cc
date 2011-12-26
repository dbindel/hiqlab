/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include "bcoordmatrix.h"


BlockCoordMatrix::BlockCoordMatrix(CoordMatrix* data, int* idx, int nblocks) :
        idx(idx), nblocks(nblocks)
{
    matrix = data;
    owns_matrix = 0;
}


BlockCoordMatrix::BlockCoordMatrix(int* idx, int nblocks) :
    idx(idx), nblocks(nblocks)
{
    matrix = new CoordMatrix(idx[nblocks]);
    owns_matrix = 1;
}


BlockCoordMatrix::~BlockCoordMatrix()
{
    if (owns_matrix)
        delete matrix;
}


void BlockCoordMatrix::reset()
{
    if (owns_matrix)
        delete matrix;
    owns_matrix = 1;
    matrix = new CoordMatrix(idx[nblocks]);
}


void BlockCoordMatrix::add(BlockCoordMatrix& src,
                           int idest, int jdest,
                           int isrc,  int jsrc, double c)
{
    matrix->add_CoordMatrix(src.matrix,
                            src.idx[isrc],    src.idx[jsrc],
                            src.idx[isrc+1]-1,src.idx[jsrc+1]-1,
                            idx[idest],idx[jdest],
                            c);
}


void BlockCoordMatrix::add(BlockCoordMatrix& src,
                           int idest, int jdest,
                           int isrc,  int jsrc, dcomplex c)
{
    matrix->add_CoordMatrix(src.matrix,
                            src.idx[isrc],    src.idx[jsrc],
                            src.idx[isrc+1]-1,src.idx[jsrc+1]-1,
                            idx[idest],idx[jdest],
                            c);
}


void BlockCoordMatrix::get_block(CoordMatrix* result, int i, int j, dcomplex c)
{
    result->add_CoordMatrix(matrix,
                            idx[i],    idx[j],
                            idx[i+1]-1,idx[j+1]-1,
                            0, 0, c);
}
