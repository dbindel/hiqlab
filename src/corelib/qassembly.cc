/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstring>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>

#include "qassembly.h"
#include "mesh.h"


QAssembler::~QAssembler()
{
}


void QGlobalAssembler::add(int* eltid, int n, dcomplex* Ke)
{
    add_helper(eltid, n, Ke);
}


void QGlobalAssembler::add(int* eltid, int n, double* Ke)
{
    add_helper(eltid, n, Ke);
}


int QGlobalAssembler::get_Vlocal(QMatrix<double>& Vl, int* eltid, int n)
{
    int ng = mesh->numglobals();
    int nzflag = 0;
    for (int j = 0; j < ng; ++j)
        for (int i = 0; i < n; ++i) {
            Vl(i,j) = mesh->shapeg(eltid[i],j);
            nzflag = nzflag || (Vl(i,j) != 0);
        }
    return nzflag;
}


template<class T>
void QGlobalAssembler::add_helper(int* eltid, int n, T* Ke)
{
    int ng = mesh->numglobals();

    if (ng == 0) {
        assembler.add(eltid, n, Ke);
        return;
    }

    QMatrix<double>   Vl(n,ng);
    if (!get_Vlocal(Vl,eltid,n)) {
        assembler.add(eltid, n, Ke);
        return;
    }

    QMatrix<T> Ke2(n+ng,n+ng);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            Ke2(i,j) = Ke[j*n+i];

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < ng; ++j)
            for (int l = 0; l < n; ++l)
                Ke2(i,j+n) += Ke2(i,l) * Vl(l,j);

    for (int i = 0; i < ng; ++i)
        for (int j = 0; j < n; ++j)
            for (int l = 0; l < n; ++l)
                Ke2(i+n,j) += Vl(l,i) * Ke2(l,j);

    for (int i = 0; i < ng; ++i)
        for (int j = 0; j < ng; ++j)
            for (int l = 0; l < n; ++l)
                Ke2(i+n,j+n) += Ke2(i+n,l) * Vl(l,j);

    std::vector<int> myid(n+ng);
    for (int i = 0; i < n;  ++i)   myid[i]   = eltid[i];
    for (int i = 0; i < ng; ++i)   myid[n+i] = mesh->globalid(i);
    assembler.add(&(myid[0]), n+ng, Ke2.data);
}


void QReduceAssembler::add(int* eltid, int n, dcomplex* Ke)
{
    std::vector<int> mapped_eltid(n);
    for (int i = 0; i < n; ++i)
        mapped_eltid[i] = mesh->id(eltid[i]);
    assembler.add(&(mapped_eltid[0]), n, Ke);
}


void QReduceAssembler::add(int* eltid, int n, double* Ke)
{
    std::vector<int> mapped_eltid(n);
    for (int i = 0; i < n; ++i)
        mapped_eltid[i] = mesh->id(eltid[i]);
    assembler.add(&mapped_eltid[0], n, Ke);
}


QStructAssembler::~QStructAssembler()
{
}


void QStructAssembler::add(int* eltid, int n)
{
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            add(eltid[i], eltid[j]);
}


void QStructAssembler::add(int* eltidm, int m, int* eltidn, int n)
{
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; ++i)
            add(eltidm[i], eltidn[j]);
}


int QGlobalStructAssembler::get_Vlocal(QMatrix<double>& Vl, int* eltid, int n)
{
    int ng = mesh->numglobals();
    int nzflag = 0;
    for (int j = 0; j < ng; ++j)
        for (int i = 0; i < n; ++i) {
            Vl(i,j) = mesh->shapeg(eltid[i],j);
            nzflag = nzflag || (Vl(i,j) != 0);
        }
    return nzflag;
}


void QGlobalStructAssembler::add(int i, int j)
{
    assert(0);
}


void QGlobalStructAssembler::add(int* eltid, int n)
{
    int ng = mesh->numglobals();

    if (ng == 0) {
        assembler.add(eltid, n);
        return;
    }

    //FIXME!?: probably should be unmapped eltid not reduced
    //         which is the case when ReduceStructAssembler
    //         calls GlobalStructAssembler -> get_Vlocal
    QMatrix<double>   Vl(n,ng);
    if (!get_Vlocal(Vl,eltid,n)) { 
        assembler.add(eltid, n);
        return;
    }

    std::vector<int> myid(n+ng);
    for (int i = 0; i < n;  ++i)   myid[i]   = eltid[i];
    for (int i = 0; i < ng; ++i)   myid[n+i] = mesh->globalid(i);
    assembler.add(&myid[0], n+ng);
}


void QReduceStructAssembler::add(int i, int j)
{
    assembler.add(mesh->id(i), mesh->id(j));
}


void QReduceStructAssembler::add(int* eltid, int n)
{
    std::vector<int> mapped_eltid(n);
    for (int i = 0; i < n; ++i)
        mapped_eltid[i] = mesh->id(eltid[i]);
    assembler.add(&mapped_eltid[0], n);
}


QVecAssembler::QVecAssembler(double* vr, double* vi, Mesh* mesh,
                             int reduced, int stride) :
    vr(vr), vi(vi), mesh(mesh), reduced(reduced), stride(stride)
{
}


QVecAssembler::~QVecAssembler()
{
}


int QVecAssembler::map(int i)
{
    return (reduced && mesh) ? mesh->id(i) : i;
}


void QVecAssembler::add(int* eltid, int n, dcomplex* Ve)
{
    for (int i = 0; i < n; ++i) {
        int k = map(eltid[i]);
        if (k >= 0) {
            if (vr) vr[k*stride] += real(Ve[i]);
            if (vi) vi[k*stride] += imag(Ve[i]);
        }
    }
}


void QVecAssembler::add(int* eltid, int n, double* Ve)
{
    for (int i = 0; i < n; ++i) {
        int k = map(eltid[i]);
        if (k >= 0) {
            if (vr) vr[k*stride] += Ve[i];
        }
    }
}


void QVecAssembler::add(int i, double eltr, double elti)
{
    int k = map(i);
    if (k >= 0) {
        if (vr) vr[k*stride] += eltr;
        if (vi) vi[k*stride] += elti;
    }
}


void QVecAssembler::set(int i, double eltr, double elti)
{
    int k = map(i);
    if (k >= 0) {
        if (vr) vr[k*stride] = eltr;
        if (vi) vi[k*stride] = elti;
    }
}


QBCAssembler::~QBCAssembler()
{
}


void QBCAssembler::add(int* eltid, int n, dcomplex* Ve)
{
    for (int i = 0; i < n; ++i) {
        mesh->bc (eltid[i]) =  type;
        mesh->bv (eltid[i]) += real(Ve[i]);
        mesh->bvi(eltid[i]) += imag(Ve[i]);
    }
}


void QBCAssembler::add(int* eltid, int n, double* V)
{
    for (int i = 0; i < n; ++i) {
        mesh->bc (eltid[i]) =  type;
        mesh->bv (eltid[i]) += V[i];
    }
}


void QBCAssembler::add(int i, double eltr, double elti)
{
    mesh->bc(i)  =  type;
    mesh->bv(i)  += eltr;
    mesh->bvi(i) += elti;
}


void QBCAssembler::set(int i, double eltr, double elti)
{
    mesh->bc(i)  = type;
    mesh->bv(i)  = eltr;
    mesh->bvi(i) = elti;
}
