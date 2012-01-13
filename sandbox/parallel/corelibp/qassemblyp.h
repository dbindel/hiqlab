/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef QASSEMBLYP_H
#define QASSEMBLYP_H

#include "qcomplex.h"
#include "qmatrix.h"
#include <vector>
#include "qassembly.h"

class Mesh_Partition;
class QAddBlockProlongator;
class QAddBlockProlongator_Partition;


/*@t ----------
 * \section{QGlobalPAssembler, QReducePAssembler, QGlobalPStructAssembler, & QReducePStructAssmbler}
 *
 * These classes are equivalent in its purpose to [QGlobalAssembler],[QGlobalStructAssembler],
 * [QReduceAssember],[QReduceStructAssembler] but instead take a [QAddBlockProlongator] object and allows one to
 * assemble a prolongator. The row indices are mapped using the TO mesh indices and the column indices
 * are mapped using the FROM mesh indices.
 *
 *@q*/


class QGlobalPAssembler : public QAssembler {
  public:
    QGlobalPAssembler(QAddBlockProlongator* prologator, QAssembler& assembler) :
        prolongator(prolongator), assembler(assembler) {}
    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, double*   Ke);
  private:
    QAddBlockProlongator* prolongator;
    QAssembler& assembler;
    template<class T>
        void add_helper(int* eltid, int n, T* Ke);
};


class QReducePAssembler : public QAssembler {
  public:
    QReducePAssembler(QAddBlockProlongator* prolongator, QAssembler& assembler) :
        prolongator(prolongator), assembler(assembler) {}
    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, double*   Ke);
  private:
    QAddBlockProlongator* prolongator;
    QAssembler& assembler;
};


class QGlobalPStructAssembler : public QStructAssembler {
  public:
    QGlobalPStructAssembler(QAddBlockProlongator* prolongator, QStructAssembler& assembler) :
        prolongator(prolongator), assembler(assembler) {}
    void add(int i, int j);
    void add(int* eltid, int n);
    void add(int* eltidm, int m, int* eltidn, int n);
  private:
    QAddBlockProlongator* prolongator;
    QStructAssembler& assembler;
};


class QReducePStructAssembler : public QStructAssembler { 
  public:
    QReducePStructAssembler(QAddBlockProlongator* prolongator, QStructAssembler& assembler) :
        prolongator(prolongator), assembler(assembler) {}
    void add(int i, int j);
    void add(int* eltid, int n);
    void add(int* eltidm, int m, int* eltidn, int n);
  private:
    QAddBlockProlongator* prolongator;
    QStructAssembler& assembler;
};


/*@t ----------
 * \section{QReduceAssemblerP & QReduceStructAssemblerP object}
 * 
 * These classes are equivalent in its purpose to [QReduceAssember],[QReduceStructAssembler]
 * but instead take a [Mesh_Partition] object and allows one to
 * assemble a local partition into a global object by mapping local
 * indices into global indices.
 *@q*/


class QReduceAssemblerP : public QAssembler {
  public:
    QReduceAssemblerP(Mesh_Partition* mesh, QAssembler& assembler) :
        mesh(mesh), assembler(assembler) {}
    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double* Ke);
  private:
    Mesh_Partition* mesh;
    QAssembler& assembler;
};


class QReduceStructAssemblerP : public QStructAssembler {
  public:
    QReduceStructAssemblerP(Mesh_Partition* mesh, QStructAssembler& assembler) :
        mesh(mesh), assembler(assembler) {}
    void add(int i, int j);
    void add(int* eltid, int n);
  private:
    Mesh_Partition* mesh;
    QStructAssembler& assembler;
};


/*@t ----------
 * \section{QReducePAssemblerP & QReducePStructAssmblerP}
 *
 * These classes are equivalent in its purpose to [QReducePAssember],[QReducePStructAssembler] 
 * but instead take a [QAddBlockProlongator_Partition] object and allows one to
 * assemble a partition of the prolongator. The row indices are mapped using the TO mesh indices 
 * and the column indices are mapped using the FROM mesh indices.
 *
 *@q*/


class QReducePAssemblerP : public QAssembler {
  public:
    QReducePAssemblerP(QAddBlockProlongator_Partition* prolongator, QAssembler& assembler) :
        prolongator(prolongator), assembler(assembler) {}
    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, double*   Ke);
  private:
    QAddBlockProlongator_Partition* prolongator;
    QAssembler& assembler;
};


class QReducePStructAssemblerP : public QStructAssembler {
  public:
    QReducePStructAssemblerP(QAddBlockProlongator_Partition* prolongator, QStructAssembler& assembler) :
        prolongator(prolongator), assembler(assembler) {}
    void add(int i, int j);
    void add(int* eltid, int n);
    void add(int* eltidm, int m, int* eltidn, int n);
  private:
    QAddBlockProlongator_Partition* prolongator;
    QStructAssembler& assembler;
};


#endif /* QASSEMBLYP_H */
