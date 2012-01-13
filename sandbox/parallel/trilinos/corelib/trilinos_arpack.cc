/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#include <cstdio>
#include <cstring>
#include <cmath>
#include <memory>
#include <iostream>

#include "qcomplex.h"
#include "cscmatrix.h"
#include "mesh.h"
#include "qmatrix.h"
#include "pmlmode.h"
#include "pareigs.h"
#include "qpassembly_trilinos.h"

#include "Epetra_MpiComm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Export.h"

#include "trilinos_arpack.h"
#include "Amesos_Operator.h"

#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_vector.h"
#include "Amesos_Operator_Complex.h"

extern Epetra_MpiComm* HIQLAB_Comm;

using namespace std;

class ArpackDS_trilinos : public PArpackDS {
public:
    ArpackDS_trilinos(Epetra_RowMatrix* Kshift, Epetra_RowMatrix* M) :
        Kshift(Kshift), M(M), Op(new Amesos_Operator(Kshift,M)) {
        set_n(Kshift->NumGlobalRows());
        make_root_map();
        x_g   = new Epetra_MultiVector(M->Map(), 1);
        opx_g = new Epetra_MultiVector(M->Map(), 1);
        export_l2g = new Epetra_Export(*root_map, M->Map());
    }
    ArpackDS_trilinos::~ArpackDS_trilinos();
protected:
    void times_OP1(double* x, double* opx);
    void broadcast(int* x, int n, int root);
    void broadcast(double* x, int n, int root);
    void barrier();
    int  mypid() { return HIQLAB_Comm->MyPID(); };
private:
    Epetra_RowMatrix* Kshift;
    Epetra_RowMatrix* M;
    Epetra_Operator*  Op;

    Epetra_Map*       root_map;
    Epetra_Export*    export_l2g;
    Epetra_MultiVector* x_g;
    Epetra_MultiVector* opx_g;
    void make_root_map();
};

ArpackDS_trilinos::~ArpackDS_trilinos()
{
    delete export_l2g;
    delete x_g;
    delete opx_g;
    delete root_map;
}

void ArpackDS_trilinos::make_root_map()
{
    if (HIQLAB_Comm->MyPID()==0)
        root_map = new Epetra_Map(-1, get_n(), 0, *HIQLAB_Comm);
    else
        root_map = new Epetra_Map(-1,       0, 0, *HIQLAB_Comm);
}

void ArpackDS_trilinos::times_OP1(double* x, double* opx)
{
    int n = root_map->NumMyElements();
    Epetra_MultiVector     x_l(View, *root_map,   x, n, 1);
    Epetra_MultiVector   opx_l(View, *root_map, opx, n, 1);
    x_g->Export(x_l, *export_l2g, Insert);
    Op->Apply(*x_g, *opx_g);
    opx_l.Import(*opx_g, *export_l2g, Insert);
}

void ArpackDS_trilinos::broadcast(int* x, int n, int root)
{
    HIQLAB_Comm->Broadcast(x, n, root);
}

void ArpackDS_trilinos::broadcast(double* x, int n, int root)
{
    HIQLAB_Comm->Broadcast(x, n, root);
}

void ArpackDS_trilinos::barrier()
{
    HIQLAB_Comm->Barrier();
}

int compute_eigs_arpack(Epetra_RowMatrix* Kshift, Epetra_RowMatrix* M,
                        int nev, int ncv,  double* d,
                        Epetra_MultiVector* vr)
{

    double* v;

    // Check to see if ncv follows criterion
    int n = M->NumGlobalRows();
    if (ncv-nev < 2) ncv = nev + 2;
    if (ncv     > n) ncv = n;

    ArpackDS_trilinos arpack(Kshift, M);
    arpack.set_nev(nev);
    arpack.set_ncv(ncv);
    arpack.set_maxitr(30);
    arpack.set_which("LM");
    arpack.set_mode(1);

    int ierr = 0;
    int NumMyElements;
    if (HIQLAB_Comm->MyPID()==0) {
        NumMyElements = n;
        v = new double[n*ncv];
    } else {
        NumMyElements = 0;
    }

    // -- Call ARPACK
    HIQLAB_Comm->Barrier();
    ierr = arpack.compute_eigs(d, v);

    // -- Distribute eigenvalues
    HIQLAB_Comm->Barrier();
    HIQLAB_Comm->Broadcast(d, nev, 0);

    // -- Distribute eigenvectors
    HIQLAB_Comm->Barrier();
    Epetra_Map Map(n, NumMyElements, 0, *HIQLAB_Comm);
    Epetra_MultiVector v_loc(View, Map, v, NumMyElements, vr->NumVectors() );

    Epetra_Export Exporter(Map, vr->Map());
    vr->Export(v_loc, Exporter, Insert);

    // -- Delete temporary
    HIQLAB_Comm->Barrier();
    if (HIQLAB_Comm->MyPID()==0)
        delete[] v;

    HIQLAB_Comm->Barrier();
    return ierr;
}


// -- ArpackDN_trilinos class

class ArpackDN_trilinos : public PArpackDN {
public:
    ArpackDN_trilinos(Epetra_RowMatrix* Kshift, Epetra_RowMatrix* M) :
        Kshift(Kshift), M(M), Op(new Amesos_Operator(Kshift,M)) {
        set_n(Kshift->NumGlobalRows());
        make_root_map();
        x_g   = new Epetra_MultiVector(M->Map(), 1);
        opx_g = new Epetra_MultiVector(M->Map(), 1);
        export_l2g = new Epetra_Export(*root_map, M->Map());
    }
    ArpackDN_trilinos::~ArpackDN_trilinos();
protected:
    void times_OP1(double* x, double* opx);
    void broadcast(int* x, int n, int root);
    void broadcast(double* x, int n, int root);
    void barrier();
    int  mypid() { return HIQLAB_Comm->MyPID(); };
private:
    Epetra_RowMatrix* Kshift;
    Epetra_RowMatrix* M;
    Epetra_Operator*  Op;

    Epetra_Map*       root_map;
    Epetra_Export*    export_l2g;
    Epetra_MultiVector* x_g;
    Epetra_MultiVector* opx_g;
    void make_root_map();
};

ArpackDN_trilinos::~ArpackDN_trilinos()
{
    delete export_l2g;
    delete x_g;
    delete opx_g;
    delete root_map;
}

void ArpackDN_trilinos::make_root_map()
{
    if (HIQLAB_Comm->MyPID()==0)
        root_map = new Epetra_Map(-1, get_n(), 0, *HIQLAB_Comm);
    else
        root_map = new Epetra_Map(-1,       0, 0, *HIQLAB_Comm);
}

void ArpackDN_trilinos::times_OP1(double* x, double* opx)
{
    int n = root_map->NumMyElements();
    Epetra_MultiVector     x_l(View, *root_map,   x, n, 1);
    Epetra_MultiVector   opx_l(View, *root_map, opx, n, 1);
    x_g->Export(x_l, *export_l2g, Insert);
    Op->Apply(*x_g, *opx_g);
    opx_l.Import(*opx_g, *export_l2g, Insert);
}

void ArpackDN_trilinos::broadcast(int* x, int n, int root)
{
    HIQLAB_Comm->Broadcast(x, n, root);
}

void ArpackDN_trilinos::broadcast(double* x, int n, int root)
{
    HIQLAB_Comm->Broadcast(x, n, root);
}

void ArpackDN_trilinos::barrier()
{
    HIQLAB_Comm->Barrier();
}

int compute_eigs_arpack(Epetra_RowMatrix* Kshift, Epetra_RowMatrix* M,
                        int nev, int ncv,  double* dr, double* di,
                        Epetra_MultiVector* vri)
{

    double* v;

    // Check to see if ncv follows criterion
    int n = M->NumGlobalRows();
    if (ncv-nev < 2) ncv = nev + 2;
    if (ncv     > n) ncv = n;

    ArpackDN_trilinos arpack(Kshift, M);
    arpack.set_nev(nev);
    arpack.set_ncv(ncv);
    arpack.set_maxitr(30);
    arpack.set_which("LM");
    arpack.set_mode(1);

    int ierr = 0;
    int NumMyElements;
    if (HIQLAB_Comm->MyPID()==0) {
        NumMyElements = n;
        v = new double[n*ncv];
    } else {
        NumMyElements = 0;
    }

    // -- Call ARPACK
    HIQLAB_Comm->Barrier();
    ierr = arpack.compute_eigs(dr, di, v);

    // -- Distribute eigenvalues
    HIQLAB_Comm->Barrier();
    HIQLAB_Comm->Broadcast(dr, nev, 0);
    HIQLAB_Comm->Broadcast(di, nev, 0);

    // -- Distribute eigenvectors
    HIQLAB_Comm->Barrier();
    Epetra_Map Map(n, NumMyElements, 0, *HIQLAB_Comm);
    Epetra_MultiVector v_loc(View, Map, v, NumMyElements, vri->NumVectors() );

    Epetra_Export Exporter(Map, vri->Map());
    vri->Export(v_loc, Exporter, Insert);

    // -- Delete temporary
    HIQLAB_Comm->Barrier();
    if (HIQLAB_Comm->MyPID()==0)
        delete[] v;

    HIQLAB_Comm->Barrier();
    return ierr;
}


// -- ArpackZN_trilinos class

class ArpackZN_trilinos : public PArpackZN {
public:
    ArpackZN_trilinos(Epetra_CrsMatrix_Complex* Kshift, Epetra_CrsMatrix_Complex* M) :
        Kshift(Kshift), M(M), Op(new Amesos_Operator_Complex(Kshift,M)) {
        set_n(Kshift->NumGlobalRows());
        make_root_map();
        x_g        = new Epetra_MultiVector_Complex(M->Map(), 1);
        opx_g      = new Epetra_MultiVector_Complex(M->Map(), 1);
        export_l2g = new Epetra_Export(*root_map, M->Map());
    }
    ArpackZN_trilinos::~ArpackZN_trilinos();
protected:
    void times_OP1(dcomplex* x, dcomplex* opx);
    void broadcast(int* x, int n, int root);
    void broadcast(double* x, int n, int root);
    void broadcast(dcomplex* x, int n, int root);
    void barrier();
    int  mypid() { return HIQLAB_Comm->MyPID(); };
private:
    Epetra_CrsMatrix_Complex* Kshift;
    Epetra_CrsMatrix_Complex* M;
    Epetra_Operator_Complex*  Op;

    Epetra_Map*       root_map;
    Epetra_Export*    export_l2g;
    Epetra_MultiVector_Complex* x_g;
    Epetra_MultiVector_Complex* opx_g;
    void make_root_map();
};

ArpackZN_trilinos::~ArpackZN_trilinos()
{
    delete Op;
    delete export_l2g;
    delete x_g;
    delete opx_g;
    delete root_map;
}

void ArpackZN_trilinos::make_root_map()
{
    if (HIQLAB_Comm->MyPID()==0)
        root_map = new Epetra_Map(-1, get_n(), 0, *HIQLAB_Comm);
    else
        root_map = new Epetra_Map(-1,       0, 0, *HIQLAB_Comm);
}

void ArpackZN_trilinos::times_OP1(dcomplex* x, dcomplex* opx)
{
    int n = root_map->NumMyElements();
    Epetra_MultiVector_Complex     x_l(*root_map,   x, n, 1);
    Epetra_MultiVector_Complex   opx_l(*root_map, opx, n, 1);
    x_g->Export          (x_l,           *export_l2g, Insert);
    x_g->get_Vz()->Export(x_l.view_Vz(), *export_l2g, Insert);
    Op->Apply(*x_g, *opx_g);
    opx_l.Import          (*opx_g,           *export_l2g, Insert);
    opx_l.get_Vz()->Import(opx_g->view_Vz(), *export_l2g, Insert);
    opx_l.ExtractCopy(opx, n);
}

void ArpackZN_trilinos::broadcast(int* x, int n, int root)
{
    HIQLAB_Comm->Broadcast(x, n, root);
}

void ArpackZN_trilinos::broadcast(double* x, int n, int root)
{
    HIQLAB_Comm->Broadcast(x, n, root);
}

void ArpackZN_trilinos::broadcast(dcomplex* x, int n, int root)
{
    double* dr = new double[n];
    double* di = new double[n];

    copy_complex(x, dr, di, n);
    broadcast(dr, n, root);
    broadcast(di, n, root);
    copy_complex(dr, di, x, n);

    barrier();
    delete[] dr;
    delete[] di;
}

void ArpackZN_trilinos::barrier()
{
    HIQLAB_Comm->Barrier();
}

int compute_eigs_arpack(Epetra_CrsMatrix_Complex* Kshift, Epetra_CrsMatrix_Complex* M,
                        int nev, int ncv,  dcomplex* dz,
                        Epetra_MultiVector_Complex* vz)
{

    dcomplex* v;
    double* dr = new double[nev];
    double* di = new double[nev];

    // Check to see if ncv follows criterion
    int n = M->NumGlobalRows();
    if (ncv-nev < 2) ncv = nev + 2;
    if (ncv     > n) ncv = n;

    ArpackZN_trilinos arpack(Kshift, M);
    arpack.set_nev(nev);
    arpack.set_ncv(ncv);
    arpack.set_maxitr(30);
    arpack.set_which("LM");
    arpack.set_mode(1);

    int ierr = 0;
    int NumMyElements;
    if (HIQLAB_Comm->MyPID()==0) {
        NumMyElements = n;
        v = new dcomplex[n*ncv];
    } else {
        NumMyElements = 0;
    }

    // -- Call ARPACK
    HIQLAB_Comm->Barrier();
    ierr = arpack.compute_eigs(dz, v);

    // -- Distribute eigenvalues
    HIQLAB_Comm->Barrier();
    if (HIQLAB_Comm->MyPID()==0)
        copy_complex(dz, dr, di, nev);
    HIQLAB_Comm->Broadcast(dr, nev, 0);
    HIQLAB_Comm->Broadcast(di, nev, 0);
    copy_complex(dr, di, dz, nev);

    // -- Distribute eigenvectors
    HIQLAB_Comm->Barrier();
    Epetra_Map Map(n, NumMyElements, 0, *HIQLAB_Comm);
    Epetra_MultiVector_Complex v_loc(Map, v, NumMyElements, vz->NumVectors() );

    Epetra_Export Exporter(Map, vz->Map());
    vz->Export          (v_loc,           Exporter, Insert);
    vz->get_Vz()->Export(v_loc.view_Vz(), Exporter, Insert);

    // -- Delete temporary
    HIQLAB_Comm->Barrier();
    delete[] dr;
    delete[] di;
    if (HIQLAB_Comm->MyPID()==0)
        delete[] v;

    HIQLAB_Comm->Barrier();
    return ierr;
}

int compute_eigs_arpack(Epetra_CrsMatrix_Complex* Kshift,
                        Epetra_CrsMatrix_Complex* M,
                        int nev, int ncv,  double* dr, double* di,
                        Epetra_MultiVector_Complex* vz)
{
    dcomplex* dz = new dcomplex[nev];

    int ierr = compute_eigs_arpack(Kshift, M, nev, ncv, dz, vz);
    copy_complex(dz, dr, di, nev);

    delete[] dz;

    return ierr;
}

// -- Solvers which treat mesh directly

int compute_eigs_arpack(Mesh* mesh, double w0, int nev, int ncv,
                        double* dr, double* di,
                        Epetra_MultiVector* vri)
{
    Epetra_CrsMatrix* Kshift;
    Epetra_CrsMatrix* M;
    dcomplex* d = new dcomplex[nev];

    Kshift = Mesh_assemble_dR_trilinos(mesh, 1.0, 0.0, -w0*w0, 1);
    M      = Mesh_assemble_dR_trilinos(mesh, 0.0, 0.0,    1.0, 1);
    int status = compute_eigs_arpack(Kshift, M, nev, ncv, dr, di, vri);
    copy_complex(dr, di, d, nev);
    for (int i = 0; i < nev; ++i) {
        d[i] = sqrt(w0*w0 + 1.0/d[i]);
    }
    copy_complex(d, dr, di, nev);

    delete[] d;
    delete M;
    delete Kshift;

    return status;

}

int compute_eigs_arpack(Mesh* mesh, double w0, int nev, int ncv,
                        double* dr, double* di,
                        Epetra_MultiVector_Complex* vz)
{
    Epetra_CrsMatrix_Complex* Kshift;
    Epetra_CrsMatrix_Complex* M;
    dcomplex* d;

    d      = new dcomplex[nev];
    Kshift = Mesh_assemble_dRz_trilinos(mesh, 1.0, 0.0, -w0*w0, 1);
    M      = Mesh_assemble_dRz_trilinos(mesh, 0.0, 0.0,    1.0, 1);
    int status = compute_eigs_arpack(Kshift, M, nev, ncv, d, vz);
    for (int i = 0; i < nev; ++i)
        d[i] = sqrt(w0*w0 + 1.0/d[i]);
    copy_complex(d, dr, di, nev);

    delete[] d;
    delete M;
    delete Kshift;

    return status;

}
