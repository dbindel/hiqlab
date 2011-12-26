/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef QASSEMBLY_H
#define QASSEMBLY_H

#include "qcomplex.h"
#include "qmatrix.h"
#include <vector>

class Mesh;


class QAssembler {
public:
    virtual ~QAssembler();
    virtual void add(int* eltid, int n, dcomplex* Ke) = 0;
    virtual void add(int* eltid, int n, double*   Ke) = 0;
    virtual void add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke) {};
    virtual void add(int* eltidm, int m, int* eltidn, int n, double*   Ke) {};
};


class QGlobalAssembler : public QAssembler {
public:
    QGlobalAssembler(Mesh* mesh, QAssembler& assembler) :
        mesh(mesh), assembler(assembler) {}
    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double* Ke);
private:
    Mesh* mesh;
    QAssembler& assembler;
    int get_Vlocal(QMatrix<double>& Vl, int* eltid, int n);
    template<class T>
        void add_helper(int* eltid, int n, T* Ke);
};


class QReduceAssembler : public QAssembler {
public:
    QReduceAssembler(Mesh* mesh, QAssembler& assembler) :
        mesh(mesh), assembler(assembler) {}
    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double* Ke);
private:
    Mesh* mesh;
    QAssembler& assembler;
};


class QStructAssembler {
public:
    virtual ~QStructAssembler();
    virtual void add(int i, int j) = 0;
    virtual void add(int* eltid, int n);
    virtual void add(int* eltidm, int m, int* eltidn, int n);
};


class QGlobalStructAssembler : public QStructAssembler {
public:
    QGlobalStructAssembler(Mesh* mesh, QStructAssembler& assembler) :
        mesh(mesh), assembler(assembler) {}
    void add(int i, int j);
    void add(int* eltid, int n);
private:
    Mesh* mesh;
    QStructAssembler& assembler;
    int get_Vlocal(QMatrix<double>& Vl, int* eltid, int n);
};


class QReduceStructAssembler : public QStructAssembler {
public:
    QReduceStructAssembler(Mesh* mesh, QStructAssembler& assembler) :
        mesh(mesh), assembler(assembler) {}
    void add(int i, int j);
    void add(int* eltid, int n);
private:
    Mesh* mesh;
    QStructAssembler& assembler;
};


class QVecAssembler : public QAssembler {
public:
    QVecAssembler(double* vr, double* vi, Mesh* mesh, int reduced = 1,
                  int stride = 1);
    ~QVecAssembler();
    void add(int* eltid, int n, dcomplex* Ve);
    void add(int* eltid, int n, double*   Ve);

    void add(int i, double eltr, double elti = 0);
    void set(int i, double eltr, double elti = 0);

    void add(int i, dcomplex elt) { add(i, real(elt), imag(elt)); }
    void set(int i, dcomplex elt) { set(i, real(elt), imag(elt)); }

private:
    double* vr;
    double* vi;
    Mesh* mesh;
    int reduced;
    int stride;

    int map(int i);
};


class QBCAssembler : public QAssembler {
public:
    QBCAssembler(Mesh* mesh, char type) : mesh(mesh), type(type) {}
    QBCAssembler(Mesh* mesh, const char* type) : mesh(mesh), type(*type) {}
    ~QBCAssembler();
    void add(int* eltid, int n, dcomplex* Ve);
    void add(int* eltid, int n, double*   Ve);

    void add(int i, double eltr, double elti = 0);
    void set(int i, double eltr, double elti = 0);

    void add(int i, dcomplex elt) { add(i, real(elt), imag(elt)); }
    void set(int i, dcomplex elt) { set(i, real(elt), imag(elt)); }

private:
    Mesh* mesh;
    char type;
};


#endif /* QASSEMBLY_H */
