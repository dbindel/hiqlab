/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstring>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>

#include "qassemblyp.h"
#include "mesh_partition.h"
#include "mesh_add_block.h"

using std::vector;


/*@t ----------
 * \section{Methods for QGlobalPAssembler, QGlobalPStructAssembler, 
 *          QReducePAssembler, QReducePStructAssembler}
 *
 *@c*/


void QGlobalPAssembler::add(int* eltid, int n, dcomplex* Ke)
{
    QGlobalPAssembler::add(eltid, n, eltid, n, Ke);
}


void QGlobalPAssembler::add(int* eltid, int n, double* Ke)
{
    QGlobalPAssembler::add(eltid, n, eltid, n, Ke);
}


void QGlobalPAssembler::add(int* eltidm, int m, int* eltidn, int n, 
                            dcomplex* Ke)
{
    assembler.add(eltidm, m, eltidn, n, Ke);
}


void QGlobalPAssembler::add(int* eltidm, int m, int* eltidn, int n, double* Ke)
{
    assembler.add(eltidm, m, eltidn, n, Ke);
}


void QReducePAssembler::add(int* eltid, int n, dcomplex* Ke)
{
    QReducePAssembler::add(eltid, n, eltid, n, Ke);
}


void QReducePAssembler::add(int* eltid, int n, double* Ke)
{
    QReducePAssembler::add(eltid, n, eltid, n, Ke);
}


void QReducePAssembler::add(int* eltidm, int m, int* eltidn, int n, 
                            dcomplex* Ke)
{
    vector<int> mapped_eltidm(m);
    vector<int> mapped_eltidn(n);
    for (int i = 0; i < m; ++i)
        mapped_eltidm[i] = prolongator->get_to()->id(eltidm[i]);
    for (int i = 0; i < n; ++i)
        mapped_eltidn[i] = prolongator->get_from()->id(eltidn[i]);
    assembler.add(&(mapped_eltidm[0]), m, &(mapped_eltidn[0]), n, Ke);
}


void QReducePAssembler::add(int* eltidm, int m, int* eltidn, int n, double* Ke)
{
    vector<int> mapped_eltidm(m);
    vector<int> mapped_eltidn(n);
    for (int i = 0; i < m; ++i)
        mapped_eltidm[i] = prolongator->get_to()->id(eltidm[i]);
    for (int i = 0; i < n; ++i)
        mapped_eltidn[i] = prolongator->get_from()->id(eltidn[i]);
    assembler.add(&(mapped_eltidm[0]), m, &(mapped_eltidn[0]), n, Ke);
}


void QGlobalPStructAssembler::add(int i, int j)
{
    assert(0);
}


void QGlobalPStructAssembler::add(int* eltid, int n)
{
    QGlobalPStructAssembler::add(eltid, n, eltid, n);
}


void QGlobalPStructAssembler::add(int* eltidm, int m, int* eltidn, int n)
{
    assembler.add(eltidm, m, eltidn, n);
}


void QReducePStructAssembler::add(int i, int j)
{
    assembler.add(prolongator->get_to()->id(i), prolongator->get_from()->id(j));
}


void QReducePStructAssembler::add(int* eltid, int n)
{
    QReducePStructAssembler::add(eltid, n, eltid, n);
}


void QReducePStructAssembler::add(int* eltidm, int m, int* eltidn, int n)
{
    vector<int> mapped_eltidm(m);
    vector<int> mapped_eltidn(n);
    for (int i = 0; i < m; ++i)
        mapped_eltidm[i] = prolongator->get_to()->id(eltidm[i]);
    for (int i = 0; i < n; ++i)
        mapped_eltidn[i] = prolongator->get_from()->id(eltidn[i]);
    assembler.add(&(mapped_eltidm[0]), m, &(mapped_eltidn[0]), n);
}


/*@t ----------
 * \section{Methods for QReduceAssemblerP, QReduceStructAssemblerP}
 *
 *@c*/


void QReduceAssemblerP::add(int* eltid, int n, dcomplex* Ke)
{
    vector<int> mapped_eltidm(n);
    vector<int> mapped_eltidn(n);

    // -- Row map(exclude off partition rows)
    for (int i = 0; i < n; ++i) {

        mapped_eltidm[i] = mesh->idg(eltid[i]);

        if (mapped_eltidm[i] < mesh->get_pstartid() || 
            mapped_eltidm[i] >= mesh->get_pend1id())
            mapped_eltidm[i] = -1;
    }

    // -- Column map(add all)
    for (int i = 0; i < n; ++i)
        mapped_eltidn[i] = mesh->idg(eltid[i]);

    assembler.add(&(mapped_eltidm[0]), n, &(mapped_eltidn[0]), n, Ke);
}


void QReduceAssemblerP::add(int* eltid, int n, double* Ke)
{
    vector<int> mapped_eltidm(n);
    vector<int> mapped_eltidn(n);

    // -- Row map(exclude off partition rows)
    for (int i = 0; i < n; ++i) {

        mapped_eltidm[i] = mesh->idg(eltid[i]);

        if (mapped_eltidm[i] < mesh->get_pstartid() || 
            mapped_eltidm[i] >= mesh->get_pend1id())
            mapped_eltidm[i] = -1;
    }

    // -- Column map(add all)
    for (int i = 0; i < n; ++i)
        mapped_eltidn[i] = mesh->idg(eltid[i]);

    assembler.add(&(mapped_eltidm[0]), n, &(mapped_eltidn[0]), n, Ke);
}


void QReduceStructAssemblerP::add(int i, int j)
{
    int idi = mesh->idg(i);
    int idj = mesh->idg(j);

    if (idi < mesh->get_pstartid() || idi >= mesh->get_pend1id())
        idi = -1;
    assembler.add(idi, idj);
}


void QReduceStructAssemblerP::add(int* eltid, int n)
{
    vector<int> mapped_eltidm(n);
    vector<int> mapped_eltidn(n);

    // -- Row map(exclude off partition rows)
    for (int i = 0; i < n; ++i) {

        mapped_eltidm[i] = mesh->idg(eltid[i]);

        if (mapped_eltidm[i] < mesh->get_pstartid() || 
            mapped_eltidm[i] >= mesh->get_pend1id())
            mapped_eltidm[i] = -1;
    }

    // -- Column map(add all)
    for (int i = 0; i < n; ++i)
        mapped_eltidn[i] = mesh->idg(eltid[i]);

    assembler.add(&(mapped_eltidm[0]), n, &(mapped_eltidn[0]), n);
}


/*@t ----------
 * \section{Methods for QReducePAssemblerP, QReducePStructAssemblerP}
 *
 *@c*/


void QReducePAssemblerP::add(int* eltid, int n, dcomplex* Ke)
{
    QReducePAssemblerP::add(eltid, n, eltid, n, Ke);
}


void QReducePAssemblerP::add(int* eltid, int n, double* Ke)
{
    QReducePAssemblerP::add(eltid, n, eltid, n, Ke);
}


void QReducePAssemblerP::add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke)
{
    vector<int> mapped_eltidm(m);
    vector<int> mapped_eltidn(n);

    // -- Row map(exclude off partition rows)
    for (int i = 0; i < m; ++i) {

        mapped_eltidm[i] = prolongator->get_to()->idg(eltidm[i]);

        if (mapped_eltidm[i] <  prolongator->get_to()->get_pstartid() || 
            mapped_eltidm[i] >= prolongator->get_to()->get_pend1id())
            mapped_eltidm[i] = -1;
    }

    // -- Column map(add all)
    for (int i = 0; i < n; ++i)
        mapped_eltidn[i] = prolongator->get_from()->idg(eltidn[i]);

    assembler.add(&(mapped_eltidm[0]), m, &(mapped_eltidn[0]), n, Ke);
}


void QReducePAssemblerP::add(int* eltidm, int m, int* eltidn, int n, double* Ke)
{
    vector<int> mapped_eltidm(m);
    vector<int> mapped_eltidn(n);

    // -- Row map(exclude off partition rows)
    for (int i = 0; i < m; ++i) {

        mapped_eltidm[i] = prolongator->get_to()->idg(eltidm[i]);

        if (mapped_eltidm[i] <  prolongator->get_to()->get_pstartid() || 
            mapped_eltidm[i] >= prolongator->get_to()->get_pend1id())
            mapped_eltidm[i] = -1;
    }

    // -- Column map(add all)
    for (int i = 0; i < n; ++i)
        mapped_eltidn[i] = prolongator->get_from()->idg(eltidn[i]);

    assembler.add(&(mapped_eltidm[0]), m, &(mapped_eltidn[0]), n, Ke);
}


void QReducePStructAssemblerP::add(int i, int j)
{
    int idi = prolongator->get_to()->idg(i);
    int idj = prolongator->get_from()->idg(j);

    if (idi < prolongator->get_to()->get_pstartid() || 
        idi >= prolongator->get_to()->get_pend1id())
        idi = -1;

    assembler.add(idi, idj);
}


void QReducePStructAssemblerP::add(int* eltid, int n)
{
    QReducePStructAssemblerP::add(eltid, n, eltid, n);
}


void QReducePStructAssemblerP::add(int* eltidm, int m, int* eltidn, int n)
{
    vector<int> mapped_eltidm(m);
    vector<int> mapped_eltidn(n);

    // -- Row map(exclude off partition rows)
    for (int i = 0; i < m; ++i) {

        mapped_eltidm[i] = prolongator->get_to()->idg(eltidm[i]);

        if (mapped_eltidm[i] <  prolongator->get_to()->get_pstartid() || 
            mapped_eltidm[i] >= prolongator->get_to()->get_pend1id())
            mapped_eltidm[i] = -1;
    }

    // -- Column map(add all)
    for (int i = 0; i < n; ++i)
        mapped_eltidn[i] = prolongator->get_from()->idg(eltidn[i]);

    assembler.add(&(mapped_eltidm[0]), m, &(mapped_eltidn[0]), n);
}
