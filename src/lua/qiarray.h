/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef QIARRAY_H
#define QIARRAY_H

class QIArray {
public:
    QIArray(int m, int n = 1, int base = 1, int* data = 0);
    QIArray(const QIArray& array);
    ~QIArray();

    QIArray& operator=(const QIArray& array);

    void    set(int i, int j, int x);
    int     get(int i, int j);

    void     wipe();
    void     copy(QIArray* A);
    QIArray* view(int i1, int i2, int j1, int j2);
    QIArray* clone();
    void     make_owner();
    void     print(char* fmt = 0);

    int     lda()    { return lda_;    }
    int     m()      { return m_;      }
    int     n()      { return n_;      }
    int     type()   { return type_;   }
    int     base()   { return base_;   }
    int*    data()   { return data_; }

private:
    int lda_, m_, n_;
    int type_;
    int base_;
    int* data_;
    int is_owner;
};

#endif /* QIARRAY_H */
