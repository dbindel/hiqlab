/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef QMATRIX_H
#define QMATRIX_H


template<class T, int m, int n>
class QMatrix1 {
public:

    QMatrix1()  {}
    ~QMatrix1() {}

    T& operator()(int i, int j) { return data[m*j+i]; }
    T& operator()(int i)        { return data[i];     }
    T& operator[](int i)        { return data[i];     }

    QMatrix1<T,m,n>& operator=(const T* Bdata) {
        for (int i = 0; i < m*n; ++i)
            data[i] = Bdata[i];
        return *this;
    }

    QMatrix1<T,m,n>& operator=(const QMatrix1<T,m,n>& B) {
        for (int i = 0; i < m*n; ++i)
            data[i] = B.data[i];
        return *this;
    }

    void clear() {
        for (int i = 0; i < m*n; ++i)
            data[i] = 0;
    }

    T data[m*n];
};


template<class T>
class QMatrix {
public:

    QMatrix(int m, int n = 1) :
        m(m), n(n) {
        owns_data = 1;
        this->data = new T[m*n];
        clear();
    }

    QMatrix(T* data, int m, int n = 1) :
        data(data), m(m), n(n) {
        if (data)
            owns_data = 0;
        else {
            owns_data = 1;
            this->data = new T[m*n];
            clear();
        }
    }

    ~QMatrix() {
        if (owns_data)
            delete[] data;
    }

    T& operator()(int i, int j) { return data[m*j+i]; }
    T& operator()(int i)        { return data[i];     }
    T& operator[](int i)        { return data[i];     }

    QMatrix<T>& operator=(const T* Bdata) {
        for (int i = 0; i < m*n; ++i)
            data[i] = Bdata[i];
        return *this;
    }

    QMatrix<T>& operator=(const QMatrix<T>& B) {
        for (int i = 0; i < m*n; ++i)
            data[i] = B.data[i];
        return *this;
    }

    void clear() {
        for (int i = 0; i < m*n; ++i)
            data[i] = 0;
    }

    T* data;
    int m, n;
    int owns_data;
};


#endif /* QMATRIX_H */
