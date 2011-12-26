/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include "pmlmode.h"
#include "pzlinear.h"
#include "mesh_csc_dR.h"
#include "meshlabeler.h"
#include "coordmatrix.h"
#include "bcoordmatrix.h"
#include "qmatrix.h"

#include <cmath>

using std::auto_ptr;


int pz_block_mesh(Mesh* mesh)
{
    MeshLabeler labeler(mesh);
    int ndf = mesh->get_ndf();
    int numid_m = labeler.label_nodal(0,ndf-1);  // Presumed mechanical
    labeler.label_nodal(ndf-1,ndf);              // Presumed potential
    labeler.label_global();
    labeler.label_branch();
    return numid_m;
}


int em_block_mesh(Mesh* mesh)
{
    MeshLabeler labeler(mesh);
    int ndf     = mesh->get_ndf();
    int numid_m = labeler.label_nodal(0,ndf-1);  // Presumed mechanical
    labeler.label_nodal(ndf-1,ndf);              // Presumed potential
    labeler.label_global();
    labeler.label_branch();
    return numid_m;
}


int em_circuit_block_mesh(Mesh* mesh)
{
    MeshLabeler labeler(mesh);
    int ndf     = mesh->get_ndf();
    int numid_m = labeler.label_nodal(0,ndf-2);  // Presumed mechanical
    labeler.label_nodal(ndf-2,ndf-1);            // Presumed potential
    labeler.label_nodal(ndf-1,ndf);              // Presumed voltage
    labeler.label_global();
    labeler.label_branch();
    return numid_m;
}


int pz_compute_eigs(Mesh* mesh, double w0, int nev, int ncv,
                    double* dr, double* di, double* vr, double* vi)
{
    int n = mesh->get_numid();
    auto_ptr<CSCMatrix> Kshift( assemble_dR(mesh, 1.0, 0.0, -w0*w0) );
    auto_ptr<CSCMatrix> M     ( assemble_dR(mesh, 0.0, 0.0, 1.0) );
    QMatrix<dcomplex> d(ncv);
    QMatrix<dcomplex> v(n,ncv);

    int status = compute_eigs(*Kshift, *M, nev, ncv, d.data, v.data, n);
    for (int i = 0; i < nev; ++i)
        d[i]  = sqrt(w0*w0 + 1.0/d[i]);
    copy_complex(d.data, dr, di, nev);
    if (vr && vi)
        copy_complex(v.data, vr, vi, n*nev);
    return status;
}


int pz_compute_eigs(Mesh* mesh, double w0, int nev, int ncv,
                    dcomplex* d, dcomplex* v)
{
    int n = mesh->get_numid();
    QMatrix<double> dr(nev);
    QMatrix<double> di(nev);
    QMatrix<double> vr(n,nev);
    QMatrix<double> vi(n,nev);

    int status = pz_compute_eigs(mesh, w0, nev, ncv, dr.data, di.data,
                                 vr.data, vi.data);
    copy_complex(dr.data, di.data, d, nev);
    if (v)
        copy_complex(vr.data, vi.data, v, n*nev);
    return status;
}


int pz_compute_eigs_mech(Mesh* mesh, double w0, int nev, int ncv,
                         double* dr, double* di, double* vr, double* vi)
{
    int n = mesh->get_numid();
    int nm = pz_block_mesh(mesh);
    int isrc[3]  = {0,nm,n};

    // Assemble and pack global K, M matrices
    BlockCoordMatrix bcoordK(isrc,2);
    BlockCoordMatrix bcoordM(isrc,2);
    mesh->assemble_dR(bcoordK.matrix, 1,0,0); bcoordK.matrix->pack();
    mesh->assemble_dR(bcoordM.matrix, 0,0,1); bcoordM.matrix->pack();

    // Form and solve the mechanical problems
    CoordMatrix coordKuu(nm);
    bcoordK.get_block(&coordKuu, 0,0, 1);
    bcoordM.get_block(&coordKuu, 0,0, -w0*w0);
    auto_ptr<CSCMatrix> Kuushift( coordKuu.to_sparse() );

    CoordMatrix coordMuu(nm);
    bcoordM.get_block(&coordMuu, 0,0, 1);
    auto_ptr<CSCMatrix> Muu( coordMuu.to_sparse() );

    QMatrix<dcomplex> d(ncv);
    QMatrix<dcomplex> v(nm,ncv);

    int status = compute_eigs(*Kuushift, *Muu, nev, ncv, d.data, v.data, nm);
    for (int i = 0; i < nev; ++i)
        d[i]  = sqrt(w0*w0 + 1.0/d[i]);
    copy_complex(d.data, dr, di, nev);
    if (vr && vi)
        copy_complex(v.data, vr, vi, nm*nev);
    return status;
}


int pz_compute_eigs_mech(Mesh* mesh, double w0, int nev, int ncv,
                         dcomplex* d, dcomplex* v)
{
    int nm = pz_block_mesh(mesh);
    QMatrix<double> dr(nev);
    QMatrix<double> di(nev);
    QMatrix<double> vr(nm,nev);
    QMatrix<double> vi(nm,nev);
    int status = pz_compute_eigs_mech(mesh, w0, nev, ncv,
                                      dr.data, di.data, vr.data, vi.data);
    copy_complex(dr.data, di.data, d, nev);
    if (v)
        copy_complex(vr.data, vi.data, v, nm*nev);
    return status;
}
