/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include "qiarray.h"
#include <cassert>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <algorithm>

using namespace std;

QIArray::QIArray(int m, int n, int base, int* data) :
    m_(m), n_(n), base_(base), data_(data)
{
    lda_ = m_;
    is_owner = 1;
    if (data_ != 0) {
        is_owner = 0;
    } else {
        data_ = new int[m*n];
        std::fill(data_, data_+m*n, 0);
    }
}


QIArray::QIArray(const QIArray& array) :
    m_(array.m_), n_(array.n_), base_(array.base_), data_(array.data_)
{
    is_owner = 0;
}


QIArray::~QIArray()
{
    if (is_owner) {
        if (data_) delete[] data_;
    }
}


QIArray& QIArray::operator=(const QIArray& array)
{
    if (is_owner) {
        if (data_) delete[] data_;
    }
    m_ = array.m_;
    n_ = array.n_;
    base_ = array.base_;
    data_ = array.data_;
    is_owner = 0;
    return *this;
}


void QIArray::set(int i, int j, int x)
{
    int k = (j-base_)*lda_ + (i-base_);
    assert(k >= 0 && k < lda_*n_);
    data_[k] = x;
}


int QIArray::get(int i, int j)
{
    int k = (j-base_)*lda_ + (i-base_);
    assert(k >= 0 && k < lda_*n_);
    return data_[k];
}


void QIArray::wipe()
{
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i)
            set(i+base_,j+base_, 0);
}


void QIArray::copy(QIArray* A)
{
    assert(m_ == A->m_ && n_ == A->n_);
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i) {
            int Ax = A->get (i+A->base_,j+A->base_);
            set(i+base_,j+base_, Ax);
        }
}


QIArray* QIArray::clone()
{
    QIArray* array = new QIArray(m_, n_, base_, data_);
    array->make_owner();
    return array;
}


QIArray* QIArray::view(int i1, int i2, int j1, int j2)
{
    assert(i1 >= base_ && i2 >= i1 && i2 < m_+base_);
    assert(j1 >= base_ && j2 >= j1 && j2 < n_+base_);
    QIArray* array = new QIArray(m_, n_, base_, data_);

    array->lda_ = lda_;
    array->m_ = i2-i1+1;
    array->n_ = j2-j1+1;
    array->data_ += ((j1-base_)*lda_ + (i1-base_));

    return array;
}


void QIArray::make_owner()
{
    if (is_owner)
        return;

    int* old_data = data_;
    data_ = new int[m_*n_];
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i)
            data_[j*m_+i] = old_data[j*lda_+i];

    lda_ = m_;
}


void QIArray::print(char* fmt)
{
    char buf[256];
    int len = 0;
    int no_data = 1;

    if (fmt == NULL)
        fmt = "%d";

    // Get the lengths of all the real and imag
    for (int j = base_; j < n_+base_; ++j)
        for (int i = base_; i < m_+base_; ++i) {
            int x = get(i,j);
            if (x != 0) no_data = 0;
            sprintf(buf, fmt, x);
            int xl = strlen(buf);
            if (len < xl)
                len = xl;
    }

    // Print the array type information
    printf("Integer array %d:%d x %d:%d\n",
           base_, m_+base_-1,
           base_, n_+base_-1);

    if (no_data) {

        // All zero array
        printf("  All zeros\n");

    } else {

        // Figure out how many columns fit on an 80-char terminal
        // Allow room for at least two space between columns + 2 at EOL
        int ncol = 78/(len+2);

        // Given the number of columns, how much room each left for spaces?
        int nspace = (78-ncol*len) / ncol;

        // Figure out where the padding should go
        len += nspace;

        // Now rip through and print the columns
        for (int j = base_; j < n_+base_; j += ncol) {

            // Print one column group
            if (j+ncol < n_+base_)
                printf(" Columns %d to %d\n", j, j+ncol-1);
            else
                printf(" Columns %d to %d\n", j, n_+base_-1);

            for (int i = base_; i < m_+base_; ++i) {

                // Print one row
                for (int jj = j; jj < j+ncol && jj < n_+base_; ++jj) {
                    sprintf(buf, fmt, get(i,jj));
                    printf("%*s", len, buf);
                }
                printf("\n");

            }
            printf("\n");

        }
    }
}
