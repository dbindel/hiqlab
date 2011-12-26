/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <vector>
#include <algorithm>

#include "coordmatrix.h"

using std::vector;
using std::fill;


static int compare_coord_col(const CoordMatrix::coord_t& xcoord, 
                             const CoordMatrix::coord_t& ycoord)
{
    if (xcoord.j < ycoord.j) return  1;
    if (xcoord.j > ycoord.j) return  0;
    if (xcoord.i < ycoord.i) return  1;
    if (xcoord.i > ycoord.i) return  0;
    return 0;
}

static int compare_coord_row(const CoordMatrix::coord_t& xcoord, 
                             const CoordMatrix::coord_t& ycoord)
{
    if (xcoord.i < ycoord.i) return  1;
    if (xcoord.i > ycoord.i) return  0;
    if (xcoord.j < ycoord.j) return  1;
    if (xcoord.j > ycoord.j) return  0;
    return 0;
}

CoordMatrix::CoordMatrix(int M, int N) :
    M(M), N(N), min_r(0), max_r(M-1), min_c(0), max_c(N-1)
{
    is_packed = 1;
    is_real = 1;
}


CoordMatrix::CoordMatrix(int N) :
    M(N), N(N), min_r(0), max_r(M-1), min_c(0), max_c(N-1)
{
    is_packed = 1;
    is_real = 1;
}


CoordMatrix::~CoordMatrix()
{
}

void CoordMatrix::restrict_range(int min_r, int max_r, int min_c, int max_c)
{
    this->min_r = min_r; this->max_r = max_r;
    this->min_c = min_c; this->max_c = max_c;
}

void CoordMatrix::restrict_range(int min_r, int max_r)
{
    this->min_r = min_r; this->max_r = max_r;
}

void CoordMatrix::add(int* eltid, int n, dcomplex* Ke)
{
    coord_t acoord;
    is_packed = 0;
    is_real = 0;

    for (int j = 0; j < n; ++j) {
        if (eltid[j] >= min_c && eltid[j] <= max_c) {
            for (int i = 0; i < n; ++i) {
                if (eltid[i] >= min_r && eltid[i] <= max_r) {
                    acoord.i = eltid[i];
                    acoord.j = eltid[j];
                    acoord.Kij = 0.0;
                    if (Ke != NULL)
                        acoord.Kij = Ke[j*n + i];
                    coord.push_back(acoord);
                }
            }
        }
    }
}


void CoordMatrix::add(int* eltid, int n, double* Ke)
{
    coord_t acoord;
    is_packed = 0;

    for (int j = 0; j < n; ++j) {
        if (eltid[j] >= min_c && eltid[j] <= max_c) {
            for (int i = 0; i < n; ++i) {
                if (eltid[i] >= min_r && eltid[i] <= max_r) {
                    acoord.i = eltid[i];
                    acoord.j = eltid[j];
                    acoord.Kij = 0.0;
                    if (Ke != NULL)
                        acoord.Kij = Ke[j*n + i];
                    coord.push_back(acoord);
                }
            }
        }
    }
}


void CoordMatrix::add(CoordMatrix* A, dcomplex coeff)
{
    if (coeff == 0.0)
        return;

    is_packed = 0;
    is_real   = 0;

    int ncoord                       = A->coord.size();
    vector<coord_t>::iterator acoord = A->coord.begin();

    for (int k = 0; k < ncoord; ++k) {
        coord_t dcoord = *acoord++;
        dcoord.Kij *= coeff;
        coord.push_back(dcoord);
    }
}


void CoordMatrix::add_dense(dcomplex* data, int i1, int j1, int i2, int j2,
                            dcomplex coeff)
{
    for (int j = j1; j <= j2; ++j) {
        for (int i = i1; i <= i2; ++i) {
            coord_t dcoord;
            dcoord.i = i;
            dcoord.j = j;
            dcoord.Kij = coeff*(*data++);
            coord.push_back(dcoord);
        }
    }
}


void CoordMatrix::add_CoordMatrix(CoordMatrix* A_sub,
                   int i_init, int j_init, int i_end, int j_end,
                   int i_dest, int j_dest, dcomplex coeff)
{
    if (coeff == 0.0)
        return;

    coord_t dcoord;
    is_packed = 0;
    is_real   = 0;

    int k = 0;
    int n = get_ncoord();

    A_sub->pack();
    int ncoord = A_sub->get_ncoord();
    vector<coord_t>::iterator acoord = A_sub->get_coord_begin();

    while ( k < ncoord && acoord->j <  j_init ) { ++acoord; ++k; }
    while ( k < ncoord && acoord->j <= j_end  ) {
        if ( acoord->i >= i_init && acoord->i <=i_end ) {
            dcoord.i   = acoord->i - i_init + i_dest;
            dcoord.j   = acoord->j - j_init + j_dest;
            dcoord.Kij = coeff * acoord->Kij;
            coord.push_back(dcoord);
            ++n;
        }
        ++acoord;
        ++k;
    }
}

void CoordMatrix::pack()
{
    if (is_packed)
        return;
    is_packed = 1;

    vector<coord_t>::iterator acoord = coord.begin();
    vector<coord_t>::iterator dest   = coord.begin();

    int k = 0;
    int n = 0;
    int i, j;

    // Sort the coordinates in column-major order
    sort(coord.begin(), coord.end(), compare_coord_col);

    int ncoord = coord.size();
    while (k < ncoord) {

        // Loop through elements with the same coordinate
        dcomplex Kij = 0.0;
        do {
            i   =  acoord->i;
            j   =  acoord->j;
            Kij += acoord->Kij;
            ++acoord;
            ++k;
        } while (k < ncoord && acoord->i == i && acoord->j == j);

        // Record element
        //FIXME: This check for Kij==0 is removed since it results in
        //       a nonzero structure that can be slightly different from
        //       what is expected without.
//        if (Kij != 0.0) {
            dest->i   = i;
            dest->j   = j;
            dest->Kij = Kij;
            ++dest;
            ++n;
//        }

    }
    coord.resize(n);
}


void CoordMatrix::to_sparse_row(int* ir, int* jc, double* pr, double* pi)
{
    pack();

    // Sort the coordinates in row-major order
    sort(coord.begin(), coord.end(), compare_coord_row);

    vector<coord_t>::iterator acoord = coord.begin();

    fill(ir, ir+N+1, 0);
    int ncoord = coord.size();
    for (int i = 0; i < ncoord; ++i) {
        if (acoord->i >= min_r && acoord->i <= max_r &&
            acoord->j >= min_c && acoord->j <= max_c) {
            ir[acoord->i+1]++;
            jc[i] = acoord->j;
            pr[i] = real(acoord->Kij);
            if (pi)
                pi[i] = imag(acoord->Kij);
            ++acoord;
        }
    }
    for (int i = 1; i < N+1; ++i) {
        ir[i] += ir[i-1];
    }
}

CSCMatrix* CoordMatrix::to_sparse()
{
    pack();
    CSCMatrix* mat = new CSCMatrix(M, N, get_ncoord(), is_real);
    to_sparse(mat->get_jc(), mat->get_ir(),
              mat->get_Ax(), (is_real ? NULL : mat->get_Az()));
    return mat;
}
