/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef QARRAY_H
#define QARRAY_H

#include <cstdio>

class QArray {
public:
    QArray(int m, int n, int type = 0, int base = 1,
           double* data_r = 0, double* data_i = 0, int lda = 0);
    QArray(const QArray& array);
    ~QArray();

    QArray& operator=(const QArray& array);

    void    set(int i, int j, double xr, double xi = 0);
    double  get(int i, int j);
    double  geti(int i, int j);

    void    wipe();
    void    copy(QArray* A);
    QArray* view(int i1, int i2, int j1, int j2);
    QArray* clone();
    QArray* transpose();
    void    make_owner();

    double  normf();
    void    print(char* fmt = 0, FILE* stream=stdout);
    void    dump(const char* fname, char* fmt=0);

    void    sub(QArray* array);
    void    add(QArray* array);
    void    sub(double x);
    void    add(double x);
    void    mul(double x);
    void    div(double x);

    int     lda()    { return lda_;    }
    int     m()      { return m_;      }
    int     n()      { return n_;      }
    int     type()   { return type_;   }
    int     base()   { return base_;   }
    double* data_r() { return data_r_; }
    double* data_i() { return data_i_; }

private:
    int lda_, m_, n_;
    int type_;
    int base_;
    double* data_r_;
    double* data_i_;
    int is_owner;
};

#endif /* QARRAY_H */
