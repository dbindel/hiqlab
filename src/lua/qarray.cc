/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include "qarray.h"
#include <cassert>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <algorithm>

using namespace std;

QArray::QArray(int m, int n, int type, int base,
               double* data_r, double* data_i, int lda) :
    m_(m), n_(n), type_(type), base_(base),
    data_r_(data_r), data_i_(data_i), lda_(lda)
{
    if (lda==0)
        lda_ = m_;
    is_owner = 1;
    if (data_r_ != 0 || data_i_ != 0) {
        is_owner = 0;
        assert(data_r_ != 0);
        assert(type != 2 || data_i_ != 0);
    } else if (type_ == 0) {
        data_r_ = new double[m*n];
        data_i_ = NULL;
        std::fill(data_r_, data_r_+m*n, 0);
    } else if (type_ == 1) {
        data_r_ = new double[2*m*n];
        data_i_ = NULL;
        std::fill(data_r_, data_r_+2*m*n, 0);
    } else if (type_ == 2) {
        data_r_ = new double[m*n];
        data_i_ = new double[m*n];
        std::fill(data_r_, data_r_+m*n, 0);
        std::fill(data_i_, data_i_+m*n, 0);
    }
}


QArray::QArray(const QArray& array) :
    m_(array.m_), n_(array.n_), type_(array.type_), base_(array.base_),
    data_r_(array.data_r_), data_i_(array.data_i_)
{
    is_owner = 0;
}


QArray::~QArray()
{
    if (is_owner) {
        if (data_r_) delete[] data_r_;
        if (data_i_) delete[] data_i_;
    }
}


QArray& QArray::operator=(const QArray& array)
{
    if (is_owner) {
        if (data_r_) delete[] data_r_;
        if (data_i_) delete[] data_i_;
    }
    m_ = array.m_;
    n_ = array.n_;
    type_ = array.type_;
    base_ = array.base_;
    data_r_ = array.data_r_;
    data_i_ = array.data_i_;
    is_owner = 0;
    return *this;
}


void QArray::set(int i, int j, double xr, double xi)
{
    int k = (j-base_)*lda_ + (i-base_);
    assert(k >= 0 && k < lda_*n_);

    if (type_ == 0) {
        data_r_[k] = xr;
    } else if (type_ == 1) {
        data_r_[2*k+0] = xr;
        data_r_[2*k+1] = xi;
    } else if (type_ == 2) {
        data_r_[k] = xr;
        data_i_[k] = xi;
    }
}


double QArray::get(int i, int j)
{
    int k = (j-base_)*lda_ + (i-base_);
    assert(k >= 0 && k < lda_*n_);

    if (type_ == 0) {
        return data_r_[k];
    } else if (type_ == 1) {
        return data_r_[2*k+0];
    } else if (type_ == 2) {
        return data_r_[k];
    } else {
        assert(0);
        return 0;
    }
}


double QArray::geti(int i, int j)
{
    int k = (j-base_)*lda_ + (i-base_);
    assert(k >= 0 && k < lda_*n_);

    if (type_ == 0) {
        return 0;
    } else if (type_ == 1) {
        return data_r_[2*k+1];
    } else if (type_ == 2) {
        return data_i_[k];
    } else {
        assert(0);
        return 0;
    }
}


void QArray::wipe()
{
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i)
            set(i+base_,j+base_, 0,0);
}


void QArray::copy(QArray* A)
{
    assert(m_ == A->m_ && n_ == A->n_);
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i) {
            double Ar = A->get (i+A->base_,j+A->base_);
            double Ai = A->geti(i+A->base_,j+A->base_);
            set(i+base_,j+base_, Ar,Ai);
        }
}


QArray* QArray::clone()
{
    QArray* array = new QArray(m_, n_, type_, base_, data_r_, data_i_);
    array->make_owner();
    return array;
}


QArray* QArray::view(int i1, int i2, int j1, int j2)
{
    assert(i1 >= base_ && i2 >= i1 && i2 < m_+base_);
    assert(j1 >= base_ && j2 >= j1 && j2 < n_+base_);
    QArray* array = new QArray(m_, n_, type_, base_, data_r_, data_i_);

    array->lda_ = lda_;
    array->m_ = i2-i1+1;
    array->n_ = j2-j1+1;

    if (type_ == 0) {
        array->data_r_ += ((j1-base_)*lda_ + (i1-base_));
    } else if (type_ == 1) {
        array->data_r_ += 2*((j1-base_)*lda_ + (i1-base_));
    } else if (type_ == 2) {
        array->data_r_ += ((j1-base_)*lda_ + (i1-base_));
        array->data_i_ += ((j1-base_)*lda_ + (i1-base_));
    }

    return array;
}


QArray* QArray::transpose()
{
    QArray* array = new QArray(n_, m_, type_, base_);
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i) {
            double Ar = get (i+base_, j+base_);
            double Ai = geti(i+base_, j+base_);
            array->set(j+base_,i+base_, Ar,Ai);
        }
    return array;
}


void QArray::make_owner()
{
    if (is_owner)
        return;

    is_owner = 1;
    double* old_data_r = data_r_;
    double* old_data_i = data_i_;

    if (type_ == 0) {
        data_r_ = new double[m_*n_];
        data_i_ = NULL;
        for (int j = 0; j < n_; ++j)
            for (int i = 0; i < m_; ++i)
                data_r_[j*m_+i] = old_data_r[j*lda_+i];
    } else if (type_ == 1) {
        data_r_ = new double[2*m_*n_];
        data_i_ = NULL;
        for (int j = 0; j < n_; ++j)
            for (int i = 0; i < m_; ++i) {
                data_r_[2*(j*m_+i)+0] = old_data_r[2*(j*lda_+i)+0];
                data_r_[2*(j*m_+i)+1] = old_data_r[2*(j*lda_+i)+1];
            }
    } else if (type_ == 2) {
        data_r_ = new double[m_*n_];
        data_i_ = new double[m_*n_];
        for (int j = 0; j < n_; ++j)
            for (int i = 0; i < m_; ++i) {
                data_r_[j*m_+i] = old_data_r[j*lda_+i];
                data_i_[j*m_+i] = old_data_i[j*lda_+i];
            }
    }

    lda_ = m_;
}


double QArray::normf()
{
    double sum = 0;
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i) {
            double Ar = get (i+base_,j+base_);
            double Ai = geti(i+base_,j+base_);
            sum += Ar*Ar + Ai*Ai;
        }
    return sqrt(sum);
}


void QArray::print(char* fmt, FILE* stream)
{
    char buf[256];
    int rlen = 0;
    int ilen = 0;
    int no_real = 1;
    int no_imag = 1;

    if (fmt == NULL)
        fmt = "%g";

    // Get the lengths of all the real and imag
    for (int j = base_; j < n_+base_; ++j)
        for (int i = base_; i < m_+base_; ++i) {

            double rx = get(i,j);
            if (rx != 0) no_real = 0;
            sprintf(buf, fmt, rx);
            int rxl = strlen(buf);
            if (rlen < rxl)
                rlen = rxl;

            double ix = geti(i,j);
            if (ix != 0) no_imag = 0;
            sprintf(buf, fmt, ix);
            int ixl = strlen(buf);
            if (rlen < ixl)
                ilen = ixl;
    }

    // Print the array type information
    if (type_ == 0) {
        fprintf(stream,"Real array %d:%d x %d:%d\n",
               base_, m_+base_-1,
               base_, n_+base_-1);
    } else if (type_ == 1) {
        fprintf(stream,"Complex array %d:%d x %d:%d (interlaced)\n",
               base_, m_+base_-1,
               base_, n_+base_-1);
    } else if (type_ == 2) {
        fprintf(stream,"Complex array %d:%d x %d:%d (separate)\n",
               base_, m_+base_-1,
               base_, n_+base_-1);
    }

    if (no_real && no_imag) {

        // All zero array
        fprintf(stream,"  All zeros\n");

    } else {

        // Figure out the length of one field
        int tlen = 0;
        if (no_real) {
            tlen = ilen+1;
        } else if (no_imag) {
            tlen = rlen;
        } else {
            tlen = rlen+ilen+4;
        }

        // Figure out how many columns fit on an 80-char terminal
        // Allow room for at least two space between columns + 2 at EOL
        int ncol = 78/(tlen+2);

        // Given the number of columns, how much room each left for spaces?
        int nspace = (78-ncol*tlen) / ncol;

        // Figure out where the padding should go
        if   (no_real) ilen += nspace;
        else           rlen += nspace;

        // Now rip through and print the columns
        for (int j = base_; j < n_+base_; j += ncol) {

            // Print one column group
            if (j+ncol < n_+base_)
                fprintf(stream," Columns %d to %d\n", j, j+ncol-1);
            else
                fprintf(stream," Columns %d to %d\n", j, n_+base_-1);

            for (int i = base_; i < m_+base_; ++i) {

                // Print one row
                for (int jj = j; jj < j+ncol && jj < n_+base_; ++jj) {
                    if (!no_real) {
                        sprintf(buf, fmt, get(i,jj));
                        fprintf(stream,"%*s", rlen, buf);
                    }
                    if (!no_real && !no_imag) {
                        fprintf(stream," + ");
                    }
                    if (!no_imag) {
                        sprintf(buf, fmt, geti(i,jj));
                        fprintf(stream,"%*si", ilen, buf);
                    }
                }
                fprintf(stream,"\n");

            }
            fprintf(stream,"\n");

        }
    }
}


void QArray::dump(const char* fname, char* fmt)
{
    FILE* fp = fopen(fname, "w+");
    print(fmt,fp);
    fclose(fp);
}


void QArray::add(QArray* A)
{
    assert(m_ == A->m_ && n_ == A->n_);
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i) {
            double mr = get (i+base_, j+base_);
            double mi = geti(i+base_, j+base_);
            double Ar = A->get (i+A->base_,j+A->base_);
            double Ai = A->geti(i+A->base_,j+A->base_);
            set(i+base_,j+base_, mr+Ar,mi+Ai);
        }
}


void QArray::sub(QArray* A)
{
    assert(m_ == A->m_ && n_ == A->n_);
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i) {
            double mr = get (i+base_, j+base_);
            double mi = geti(i+base_, j+base_);
            double Ar = A->get (i+A->base_,j+A->base_);
            double Ai = A->geti(i+A->base_,j+A->base_);
            set(i+base_,j+base_, mr-Ar,mi-Ai);
        }
}


void QArray::add(double x)
{
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i) {
            double mr = get (i+base_, j+base_);
            double mi = geti(i+base_, j+base_);
            set(i+base_,j+base_, mr+x,mi);
        }
}


void QArray::sub(double x)
{
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i) {
            double mr = get (i+base_, j+base_);
            double mi = geti(i+base_, j+base_);
            set(i+base_,j+base_, mr-x,mi);
        }
}


void QArray::mul(double x)
{
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i) {
            double mr = get (i+base_, j+base_);
            double mi = geti(i+base_, j+base_);
            set(i+base_,j+base_, mr*x,mi*x);
        }
}


void QArray::div(double x)
{
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < m_; ++i) {
            double mr = get (i+base_, j+base_);
            double mi = geti(i+base_, j+base_);
            set(i+base_,j+base_, mr/x,mi/x);
        }
}
