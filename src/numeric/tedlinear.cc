/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include "pmlmode.h"
#include "tedlinear.h"
#include "meshlabeler.h"
#include "mesh_csc_dR.h"
#include "qmatrix.h"
#include "cscmatrix.h"
#include "umfmatrix.h"
#include "bcoordmatrix.h"

#include <cmath>
#include <memory>

using std::vector;
using std::auto_ptr;


int ted_block_mesh(Mesh* mesh)
{
    MeshLabeler labeler(mesh);
    int ndf = mesh->get_ndf();
    int numid_m = labeler.label_nodal(0,ndf-1);  // Presumed mechanical
    labeler.label_nodal(ndf-1,ndf);              // Presumed thermal
    labeler.label_global();
    labeler.label_branch();
    return numid_m;
}


void ted_assemble(CoordMatrix* AB, Mesh* mesh, double T0, double cT,
                  dcomplex ca, dcomplex cb)
{
    int n  = mesh->get_numid();
    int nm = ted_block_mesh(mesh);
    int idest[4] = {0,nm,nm+nm,nm+n};
    int isrc[3]  = {0,nm,n};

    BlockCoordMatrix bcoordAB(AB, idest,3);
    BlockCoordMatrix bcoordK(isrc,2);

    if (T0 == 0) {

        // -------- Linearization 1 -----------------------------

        mesh->assemble_dR(bcoordK.matrix, 0,0,1);
        bcoordK.matrix->pack();
        bcoordAB.add(bcoordK, 0,0, 0,0,  cb/cT/cT);      // Buu +=  Muu
        bcoordAB.add(bcoordK, 1,1, 0,0,  cb/cT/cT);      // Bvv +=  Muu
        bcoordAB.add(bcoordK, 0,1, 0,0, -ca/cT/cT);      // Auv += -Muu
        bcoordK.reset();

        mesh->assemble_dR(bcoordK.matrix, 0,1,0);
        bcoordK.matrix->pack();
        bcoordAB.add(bcoordK, 2,1, 1,0,  ca/cT);      // Atv +=  Ctu
        bcoordAB.add(bcoordK, 2,2, 1,1,  cb/cT);      // Btt +=  Ctt
        bcoordK.reset();

        mesh->assemble_dR(bcoordK.matrix, 1,0,0);
        bcoordK.matrix->pack();
        bcoordAB.add(bcoordK, 1,0, 0,0,  ca);      // Avu +=  Kuu
        bcoordAB.add(bcoordK, 1,2, 0,1,  ca);      // Avt +=  Kut
        bcoordAB.add(bcoordK, 2,2, 1,1,  ca);      // Att +=  Ktt

    } else {

        // -------- Linearization 2 -----------------------------

        mesh->assemble_dR(bcoordK.matrix, 0,0,1);
        bcoordAB.add(bcoordK, 1,1, 0,0,  cb/cT/cT);      // Bvv +=  Muu
        bcoordK.reset();

        mesh->assemble_dR(bcoordK.matrix, 0,1,0);
        bcoordAB.add(bcoordK, 2,1, 1,0,  -ca/T0/cT);  // Atv += -Ctu/T0
        bcoordAB.add(bcoordK, 2,2, 1,1,  -cb/T0/cT);  // Btt += -Ctt/T0
        bcoordK.reset();

        mesh->assemble_dR(bcoordK.matrix, 1,0,0);
        bcoordAB.add(bcoordK, 0,0, 0,0, -cb);      // Buu += -Kuu
        bcoordAB.add(bcoordK, 0,1, 0,0,  ca);      // Auv +=  Kuu
        bcoordAB.add(bcoordK, 1,0, 0,0,  ca);      // Avu +=  Kuu
        bcoordAB.add(bcoordK, 1,2, 0,1,  ca);      // Avt +=  Kut
        bcoordAB.add(bcoordK, 2,2, 1,1, -ca/T0);   // Att += -Ktt/T0

    }
    AB->pack();
}


CSCMatrix* ted_assemble(Mesh* mesh, double T0, double cT,
                        dcomplex ca, dcomplex cb)
{
    int n  = mesh->get_numid();
    int nm = ted_block_mesh(mesh);
    int nt = n - nm;
    CoordMatrix coordAB(nm+n);
    vector<double> uscale(nm+n);
    vector<double> fscale(nm+n);
    ted_assemble(&coordAB, mesh, T0, cT, ca, cb);
    mesh->get_scaling_vector(&(uscale[0]),  &(fscale[0]),  true);
    mesh->get_scaling_vector(&(uscale[nm]), &(fscale[nm]), true);
    CSCMatrix* result = coordAB.to_sparse();
    result->scale(&(uscale[0]), &(fscale[0]));
    return result;
}


int ted_compute_eigs(Mesh* mesh, double w0, double T0, double cT,
                     int nev, int ncv,
                     double* dr, double* di, double* vr, double* vi)
{
    int n = mesh->get_numid();
    int nm = ted_block_mesh(mesh);

    // Assemble shifted A matrix, B matrix
    auto_ptr<CSCMatrix> Ashift( ted_assemble(mesh, T0, cT, 1,
                                             dcomplex(0,w0*cT)) );
    auto_ptr<CSCMatrix> B     ( ted_assemble(mesh, T0, cT, 0, 1) );
    QMatrix<dcomplex> d(ncv);
    QMatrix<dcomplex> v(n+nm,ncv);
    int status = compute_eigs(*Ashift, *B, nev, ncv, d.data, v.data, n+nm);

    // FIXME: This should be done once, not three times!
    vector<double> uscale(nm+n);
    mesh->get_scaling_vector(&(uscale[0]), NULL, true);
    mesh->get_scaling_vector(&(uscale[nm]), NULL, true);

    for (int i = 0; i < nev; ++i) {
        d[i]  = w0 + dcomplex(0,1)/d[i]/cT;
        for (int j = 0; j < nm+n; ++j)
            v(j,i) *= uscale[j];
    }

    copy_complex(d.data, dr, di, nev);
    if (vr && vi) {
        for (int j = 0; j < nev; ++j) {
            copy_complex(&(v(0,j)),    vr+   j*n, vi+   j*n, nm);
            copy_complex(&(v(2*nm,j)), vr+nm+j*n, vi+nm+j*n, n-nm);
        }
    }
    return status;
}


int ted_compute_eigs(Mesh* mesh, double w0, double T0, double cT,
                     int nev, int ncv,
                     dcomplex* d, dcomplex* v)
{
    int n = mesh->get_numid();
    QMatrix<double> dr(nev);
    QMatrix<double> di(nev);
    QMatrix<double> vr(n,nev);
    QMatrix<double> vi(n,nev);
    int status = ted_compute_eigs(mesh, w0, T0, cT, nev, ncv,
                                  dr.data, di.data, vr.data, vi.data);
    copy_complex(dr.data, di.data, d, nev);
    if (v)
        copy_complex(vr.data, vi.data, v, n*nev);
    return status;
}


int ted_compute_eigsp(Mesh* mesh, double w0, int nev, int ncv,
                      double* dr, double* di, double* vr, double* vi)
{
    int n = mesh->get_numid();
    int nm = ted_block_mesh(mesh);
    int nt = n-nm;
    int isrc[3]  = {0,nm,n};

    // Assemble and pack global K, C, M matrices
    BlockCoordMatrix bcoordK(isrc,2);
    BlockCoordMatrix bcoordC(isrc,2);
    BlockCoordMatrix bcoordM(isrc,2);
    mesh->assemble_dR(bcoordK.matrix, 1,0,0); bcoordK.matrix->pack();
    mesh->assemble_dR(bcoordC.matrix, 0,1,0); bcoordC.matrix->pack();
    mesh->assemble_dR(bcoordM.matrix, 0,0,1); bcoordM.matrix->pack();

    // Form and solve the mechanical problems
    CoordMatrix* coordKuu = new CoordMatrix(nm);
    bcoordK.get_block(coordKuu, 0,0, 1);
    bcoordM.get_block(coordKuu, 0,0, -w0*w0);
    CSCMatrix* Kuushift = coordKuu->to_sparse();
    delete coordKuu;

    CoordMatrix* coordMuu = new CoordMatrix(nm);
    bcoordM.get_block(coordMuu, 0,0, 1);
    CSCMatrix* Muu = coordMuu->to_sparse();
    delete coordMuu;

    QMatrix<dcomplex> d(ncv);
    QMatrix<dcomplex> v(n,ncv);
    QMatrix<dcomplex> scratch(2*n);
    compute_eigs(*Kuushift, *Muu, nev, ncv, d.data, v.data, n);
    for (int i = 0; i < nev; ++i)
        d[i] = sqrt(w0*w0 + 1.0/d[i]);

    delete Kuushift;

    // Compute the resulting temperature fields and perturb
    CoordMatrix* coordCtu = new CoordMatrix(nt, nm);
    bcoordC.get_block(coordCtu, 1,0, -1);
    CSCMatrix* Ctu = coordCtu->to_sparse();
    delete coordCtu;

    CoordMatrix* coordKut = new CoordMatrix(nm, nt);
    bcoordK.get_block(coordKut, 0,1, -1);
    CSCMatrix* Kut = coordKut->to_sparse();
    delete coordKut;

    for (int i = 0; i < nev; ++i) {

        // Compute the resulting temperature fields
        CoordMatrix* coordCtt = new CoordMatrix(nt);
        bcoordC.get_block(coordCtt, 1,1, 1);
        bcoordK.get_block(coordCtt, 1,1, -dcomplex(0,1/real(d[i])));
        CSCMatrix* Ctt = coordCtt->to_sparse();
        UMFMatrix* Cttf = new UMFMatrix(Ctt);
        delete coordCtt;

        Ctu->apply(&(v(0,i)),  scratch.data);
        Cttf->solve(&(v(nm,i)), scratch.data);
        delete Cttf;
        delete Ctt;

        // Solve bordered system for the perturbation
        CoordMatrix* coordKuub = new CoordMatrix(nm+1);
        bcoordK.get_block(coordKuub, 0,0, 1);
        bcoordM.get_block(coordKuub, 0,0, -w0*w0);
        Muu->apply(&(v[i*n]), scratch.data);
        coordKuub->add_dense(&(v(0,i)),    nm,0, nm,nm-1,  1e-6  );
        coordKuub->add_dense(scratch.data, 0,nm, nm-1,nm, -2.0*d[i]);
        CSCMatrix* Kuub = coordKuub->to_sparse();
        UMFMatrix* Kuubf = new UMFMatrix(Kuub);
        delete coordKuub;

        scratch[nm] = 0;
        Kut->apply(&(v(nm,i)),     scratch.data);
        Kuubf->solve(&(scratch[n]), scratch.data);
        delete Kuubf;
        delete Kuub;

        // Update u and omega
        for (int j = 0; j < nm; ++j)
            v[i*n+j] += scratch[n+j];
        d[i] += scratch[n+nm];

        dr[i] = real(d[i]);
        di[i] = imag(d[i]);
        if (vr && vi)
            copy_complex(&(v(0,i)), vr+i*n, vi+i*n, n);
    }

    delete Kut;
    delete Ctu;
    delete Muu;

    return 0;
}


int ted_compute_eigsp(Mesh* mesh, double w0, int nev, int ncv,
                      dcomplex* d, dcomplex* v)
{
    int n = mesh->get_numid();
    QMatrix<double> dr(nev);
    QMatrix<double> di(nev);
    QMatrix<double> vr(n,nev);
    QMatrix<double> vi(n,nev);
    int status = ted_compute_eigsp(mesh, w0, nev, ncv,
                                   dr.data, di.data, vr.data, vi.data);
    copy_complex(dr.data, di.data, d, nev);
    if (v)
        copy_complex(vr.data, vi.data, v, n*nev);
    return status;
}


int ted_compute_eigs_mech(Mesh* mesh, double w0, int nev, int ncv,
                          double* dr, double* di, double* vr, double* vi)
{
    int n = mesh->get_numid();
    int nm = ted_block_mesh(mesh);
    int isrc[3]  = {0,nm,n};

    // Assemble and pack global K, M matrices
    BlockCoordMatrix bcoordK(isrc,2);
    BlockCoordMatrix bcoordM(isrc,2);
    mesh->assemble_dR(bcoordK.matrix, 1,0,0); bcoordK.matrix->pack();
    mesh->assemble_dR(bcoordM.matrix, 0,0,1); bcoordM.matrix->pack();

    // Form and solve the mechanical problems
    CoordMatrix* coordKuu = new CoordMatrix(nm);
    bcoordK.get_block(coordKuu, 0,0, 1);
    bcoordM.get_block(coordKuu, 0,0, -w0*w0);
    CSCMatrix* Kuushift = coordKuu->to_sparse();
    delete coordKuu;

    CoordMatrix* coordMuu = new CoordMatrix(nm);
    bcoordM.get_block(coordMuu, 0,0, 1);
    CSCMatrix* Muu = coordMuu->to_sparse();
    delete coordMuu;

    QMatrix<dcomplex> d(ncv);
    QMatrix<dcomplex> v(nm,ncv);
    int status = compute_eigs(*Kuushift, *Muu, nev, ncv, d.data, v.data, nm);

    for (int i = 0; i < nev; ++i)
        d[i] = sqrt(w0*w0 + 1.0/d[i]);
    copy_complex(d.data, dr, di, nev);
    if (vr && vi)
        copy_complex(v.data, vr, vi, nm*nev);

    delete Kuushift;
    delete Muu;
    return status;
}


int ted_compute_eigs_mech(Mesh* mesh, double w0, int nev, int ncv,
                          dcomplex* d, dcomplex* v)
{
    int nm = ted_block_mesh(mesh);
    QMatrix<double> dr(nev);
    QMatrix<double> di(nev);
    QMatrix<double> vr(nm,nev);
    QMatrix<double> vi(nm,nev);
    int status = ted_compute_eigs_mech(mesh, w0, nev, ncv,
                                       dr.data, di.data, vr.data, vi.data);
    for (int i = 0; i < nev; ++i)
        d[i] = dcomplex(dr[i], di[i]);
    copy_complex(dr.data, di.data, d, nev);
    if (v)
        copy_complex(vr.data, vi.data, v, nm*nev);
    return status;
}
