/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef BCOORDMATRIX_H
#define BCOORDMATRIX_H

#include "qcomplex.h"
#include "coordmatrix.h"


/** Block matrix accessor
 */
class BlockCoordMatrix {
public:
    BlockCoordMatrix(CoordMatrix* data, int* idx, int nblocks);
    BlockCoordMatrix(int* idx, int nblocks);
    ~BlockCoordMatrix();

    void reset();
    void add(BlockCoordMatrix& src,
             int idest, int jdest,
             int isrc,  int jsrc, double c);
    void add(BlockCoordMatrix& src,
             int idest, int jdest,
             int isrc,  int jsrc, dcomplex c);
    void get_block(CoordMatrix* m, int i, int j, dcomplex c = 1.0);

    CoordMatrix* matrix;
    int* idx;
    int nblocks;
    int owns_matrix;
};


#endif /* BCOORDMATRIX_H */
