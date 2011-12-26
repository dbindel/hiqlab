/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef COORDMATRIX_H
#define COORDMATRIX_H

#include "qcomplex.h"
#include "cscmatrix.h"
#include "qassembly.h"

#include <vector>
#include <algorithm>


/** Matrix assembler for complex mass and stiffnesses.
 */
class CoordMatrix : public QAssembler {
 public:

    struct coord_t {
        int i;
        int j;
        dcomplex Kij;
    };

    /** Create a new assembler for mass and stiffness matrices
     *
     * @param N  dimension of the system
     */
    CoordMatrix(int N);
    CoordMatrix(int M, int N);
    ~CoordMatrix();

    /** Restrict the range of indices to be summed into matrices
     *
     * @param min_r  Minimum row index
     * @param max_r  Maximum row index
     * @param min_c  Minimum column index
     * @param max_c  Maximum column index
     *
     * By default this is set to [ 0, M-1, 0, N-1] in order
     */
    void restrict_range(int min_r, int max_r, int min_c, int max_c);
    void restrict_range(int min_r, int max_r);

    /** Add local mass and stiffness contributions
     *
     * @param eltid  Global identifiers for element dofs (n)
     * @param n      Number of element dofs
     * @param Ke     Element matrix (n-by-n)
     */
    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double*   Ke);
    void add(CoordMatrix* A, dcomplex coeff = 1.0);

    void add_dense(dcomplex* data, int i1, int j1, int i2, int j2,
                   dcomplex coeff = 1.0);

    /** Add a submatrix of another CoordMatrix to THIS CoordMatrix
     *
     * @param A_sub   CoordMatrix being added to THIS
     * @param i_init  Initial Row index of the submatrix
     * @param j_init  Initial Col index of the submatrix
     * @param i_end   Ending  Row index of the submatrix
     * @param j_end   Ending  Col index of the submatrix
     * @param i_dest  Initial Row index of destination
     * @param j_dest  Initial Row index of destination
     * @param mult    Scalar multiplier of the adding matrix
     */
     void add_CoordMatrix(CoordMatrix* A_sub, int i_init, int j_init,
                        int i_end, int j_end, int i_dest, int j_dest,
                        dcomplex coeff);

    /** Pack the mass and stiffness matrices so there is only
     *  one entry per coordinate.  After packing, the coordinates
     *  are ordered in increasing column-major order.
     */
    void pack();

    /** Generate CSC or standard coordinate representations.
     */
    template<class idx_t, class real_t>
    void to_sparse(idx_t* jc, idx_t* ir, real_t* pr, real_t* pi);

    CSCMatrix* to_sparse();

    /** Generate row-major standard coordinate representations.
     */
    void to_sparse_row(int* ir, int* jc, double* pr, double* pi);

    int      get_ncoord() { return coord.size(); }
    int      get_M()      { return M; }
    int      get_N()      { return N; }

    std::vector<coord_t>::iterator get_coord_begin() { return coord.begin(); }
    std::vector<coord_t>::iterator get_coord_end()   { return coord.end();   }

 private:
    std::vector<coord_t> coord;
    int M, N;
    int min_r, max_r, min_c, max_c;
    int is_packed;
    int is_real;
};


template<class idx_t, class real_t>
void CoordMatrix::to_sparse(idx_t* jc, idx_t* ir, real_t* pr, real_t* pi)
{
    pack();
    std::vector<coord_t>::iterator acoord = coord.begin();

    std::fill(jc, jc+N+1, 0);
    int ncoord = coord.size();
    for (int i = 0; i < ncoord; ++i) {
        if (acoord->i >= min_r && acoord->i <= max_r &&
            acoord->j >= min_c && acoord->j <= max_c) {
            jc[acoord->j+1]++;
            ir[i] = (idx_t) acoord->i;
            pr[i] = (real_t) real(acoord->Kij);
            if (pi)
                pi[i] = (real_t) imag(acoord->Kij);
            ++acoord;
        }
    }
    for (int i = 1; i < N+1; ++i) {
        jc[i] += jc[i-1];
    }
}


#endif /* COORDMATRIX_H */
