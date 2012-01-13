/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cmath>
#include <memory>

#include "qcomplex.h"
#include "cscmatrix.h"
#include "umfmatrix.h"
#include "mesh_csc_dR.h"
#include "mesh.h"
#include "qmatrix.h"
#include "pmlmode.h"
#include "areigs.h"

using std::auto_ptr;


class ArpackPML : public ArpackZN {
public:
    ArpackPML(CSCMatrix& Kshift, CSCMatrix& M) :
        Kshift(&Kshift), M(M), scratch(Kshift.get_n()) {
        set_n(Kshift.get_n());
    }
protected:
    void times_OP1(dcomplex* x, dcomplex* opx);
private:
    UMFMatrix Kshift;
    CSCMatrix& M;
    QMatrix<dcomplex> scratch;
};


void ArpackPML::times_OP1(dcomplex* x, dcomplex* opx)
{
    M.apply(x, scratch.data);
    Kshift.solve(opx, scratch.data);
}


int compute_eigs(CSCMatrix& Kshift, CSCMatrix& M,  int nev, int ncv,
                 dcomplex* d, dcomplex* v, int ldv)
{
    // Check to see if ncv follows criterion
    int n = M.get_n();
    if (ncv-nev < 2) ncv = nev + 2;
    if (ncv     > n) ncv = n;

    ArpackPML arpack(Kshift, M);
    arpack.set_nev(nev);
    arpack.set_ncv(ncv);
    arpack.set_maxitr(30);
    arpack.set_which("LM");
    arpack.set_mode(1);
    return arpack.compute_eigs(d, v);
}


int compute_eigs(Mesh* mesh, double w0, int nev, int ncv,
                 double* dr, double* di, double* vr, double* vi)
{
    int n = mesh->get_numid();
    QMatrix<dcomplex> d(ncv);
    QMatrix<dcomplex> v(n,ncv);
    auto_ptr<CSCMatrix> Kshift( assemble_dR(mesh, 1.0, 0.0, -w0*w0) );
    auto_ptr<CSCMatrix> M     ( assemble_dR(mesh, 0.0, 0.0, 1.0) );
    int status = compute_eigs(*Kshift, *M, nev, ncv, d.data, v.data, n);
    for (int i = 0; i < nev; ++i)
        d[i]  = sqrt(w0*w0 + 1.0/d[i]);
    copy_complex(d.data, dr, di, nev);
    if (vr && vi)
        copy_complex(v.data, vr, vi, n*nev);
    return status;
}
