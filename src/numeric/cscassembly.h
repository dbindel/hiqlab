#ifndef CSCASSEMBLY_H
#define CSCASSEMBLY_H

#include "mesh.h"
#include "cscmatrix.h"

/*@T
 * \section{Compressed sparse column assembly}
 *
 * The [[CSCAssembler]] class assembles the tangent stiffness into an
 * existing compressed sparse column matrix.  It is a checked error
 * to assemble into a slot that's not allocated in the index structure
 * going in.~\footnote{
 *   The compressed sparse column matrix class should eventually have
 *   a cache for entries that aren't supposed to go into preallocated slots.
 * }
 *
 *@c*/

class CSCAssembler : public QAssembler {
public:
    CSCAssembler(CSCMatrix* matrix);
    CSCAssembler(int* jc, int* ir, double* pr, double* pi);
    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, dcomplex* Ke);
    void add(int* eltidm, int m, int* eltidn, int n, double* Ke);
private:
    int mm;
    int nn;
    int* jc;
    int* ir;
    double* pr;
    double* pi;
    void add(int i, int j, dcomplex Ke);
    void add(int i, int j, double Ke);
};

void assemble_dR(Mesh* mesh, CSCMatrix* K, 
                 double cx = 1, double cv = 0, double ca = 0,
                 int reduced = 1);

#endif /* CSCASSEMBLY_H */
