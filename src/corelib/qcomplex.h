/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef QCOMPLEX_H
#define QCOMPLEX_H

#include <complex>
typedef std::complex<double> dcomplex;


inline void copy_complex(dcomplex* src, double* destr, double* desti, int n)
{
    for (int i = 0; i < n; ++i) {
        destr[i] = real(src[i]);
        desti[i] = imag(src[i]);
    }
}


inline void copy_complex(double* srcr, double* srci, dcomplex* dest, int n)
{
    for (int i = 0; i < n; ++i)
        dest[i] = dcomplex(srcr[i], srci[i]);
}


#endif /* QCOMPLEX_H */
