/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "qmatrix.h"
#include "qcomplex.h"
#include "model_circuit.h"
#include "mesh.h"

#define NEN        2
#define NENE       1


/* ==== Resistor Model ==== */


Resistor::Resistor(double resist, double resist_i) :
    Element(1), resist(resist), resist_i(resist_i)
{
}


Resistor::~Resistor()
{
}


void Resistor::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                           double cx, double cv, double ca)
{
    int id[NEN];
    set_local_id(mesh, eltid, id);
    QMatrix<dcomplex> Ke(NULL, NEN, NEN);
    dcomplex resist_c(resist,resist_i);
    if (cx) {
        dcomplex rval = 1.0/resist_c * cx;
        Ke(0,0) += rval;  Ke(0,1) -= rval;
        Ke(1,0) -= rval;  Ke(1,1) += rval;
    }
    K->add(id, NEN, Ke.data);
}


void Resistor::assemble_R(Mesh* mesh, int eltid)
{
    QMatrix<dcomplex> Fe(NULL, 1, NEN);
    dcomplex resist_c(resist,resist_i);
    dcomplex nodeu[NEN];
    set_local_u(mesh, eltid, nodeu);
    Fe(0) += ( nodeu[  0]-nodeu[1])/resist_c;
    Fe(1) += (-nodeu[  0]+nodeu[1])/resist_c;
    add_local_f(mesh,eltid, Fe.data);
}


/* ==== Capacitor Model ==== */


Capacitor::Capacitor(double cap, double cap_i) :
    Element(1), cap(cap), cap_i(cap_i)
{
}


Capacitor::~Capacitor()
{
}


void Capacitor::initialize(Mesh* mesh, int eltid)
{
}


void Capacitor::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                            double cx, double cv, double ca)
{
    int id[NEN];
    set_local_id(mesh, eltid, id);
    QMatrix<dcomplex> Ke(NULL, NEN, NEN);
    dcomplex cap_c(cap,cap_i);
    if (cv) {
        dcomplex rval = cap_c * cv;
        Ke(0,0) += rval;  Ke(0,1) -= rval;
        Ke(1,0) -= rval;  Ke(1,1) += rval;
    }
    K->add(id, NEN, Ke.data);
}


void Capacitor::assemble_R(Mesh* mesh, int eltid)
{
    QMatrix<dcomplex> Fe(NULL, 1, NEN);
    dcomplex cap_c(cap,cap_i);
    dcomplex nodeu[NEN];
    set_local_u(mesh, eltid, nodeu);
    Fe(0) +=  cap_c*(  nodeu[0] - nodeu[1]);
    Fe(1) +=  cap_c*( -nodeu[0] + nodeu[1]);
    add_local_f(mesh,eltid, Fe.data);
}


/* ==== Capacitor Model ==== */


Capacitor1::Capacitor1(double cap, double cap_i) :
    Element(1), cap(cap), cap_i(cap_i)
{
}


Capacitor1::~Capacitor1()
{
}


void Capacitor1::initialize(Mesh* mesh, int eltid)
{
    mesh->nbranch(eltid) = 1;
}


void Capacitor1::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                             double cx, double cv, double ca)
{
    int id[NEN+1];
    set_local_id(mesh, eltid, id);
    QMatrix<dcomplex> Ke(NULL, NEN+1, NEN+1);
    dcomplex cap_c(cap, cap_i);
    if (cx) {
        Ke(2,0) +=  cap_c * cx;
        Ke(2,1) += -cap_c * cx;
        Ke(2,2) += -1.0   * cx;
    }
    if (cv) {
        Ke(0,2) +=  1.0 * cv;
        Ke(1,2) += -1.0 * cv;
    }
    K->add(id, NEN+1, Ke.data);
}


void Capacitor1::assemble_R(Mesh* mesh, int eltid)
{
    QMatrix<dcomplex> Fe(NULL, 1, NEN+1);
    dcomplex cap_c(cap, cap_i);
    dcomplex nodeu[NEN+1],nodev[NEN+1];
    set_local_u(mesh, eltid, nodeu);
    set_local_v(mesh, eltid, nodev);
    Fe(0) += nodev[2];
    Fe(1) -= nodev[2];
    Fe(2) += (cap_c*(nodeu[0] - nodeu[1]) - nodeu[2]);
    add_local_f(mesh,eltid, Fe.data);
}


/* ==== Capacitor2 Model ==== */


Capacitor2::Capacitor2(double cap, double cap_i) :
    Element(1), cap(cap), cap_i(cap_i)
{
}


Capacitor2::~Capacitor2()
{
}


void Capacitor2::initialize(Mesh* mesh, int eltid)
{
    mesh->nbranch(eltid) = 2;
}


void Capacitor2::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                             double cx, double cv, double ca)
{
    int id[NEN+2];
    set_local_id(mesh, eltid, id);
    QMatrix<dcomplex> Ke(NULL, NEN+2, NEN+2);
    dcomplex cap_c(cap, cap_i);
    if (cx) {
        Ke(0,3) +=  1.0 * cx;
        Ke(1,3) += -1.0 * cx;
        Ke(2,3) += -1.0 * cx;

        Ke(3,0) +=  cap_c * cx;
        Ke(3,1) += -cap_c * cx;
        Ke(3,2) += -1.0 * cx;
    }
    if (cv) {
        Ke(2,2) +=  1.0 * cv;
    }
    K->add(id, NEN+2, Ke.data);

}


void Capacitor2::assemble_R(Mesh* mesh, int eltid)
{
    QMatrix<dcomplex> Fe(NULL, 1, NEN+2);
    dcomplex cap_c(cap, cap_i);
    dcomplex nodeu[NEN+2],nodev[NEN+2];
    set_local_u(mesh, eltid, nodeu);
    set_local_v(mesh, eltid, nodev);

    Fe(0) +=   nodeu[3];
    Fe(1) +=  -nodeu[3];
    Fe(2) += (-nodeu[3] + nodev[2]);
    Fe(3) += (cap_c*(nodeu[0] - nodeu[1]) - nodeu[2]);

    add_local_f(mesh,eltid, Fe.data);

}


/* ==== Inductor Model ==== */


Inductor::Inductor(double induct, double induct_i) :
    Element(1), induct(induct), induct_i(induct_i)
{
}


Inductor::~Inductor()
{
}


void Inductor::initialize(Mesh* mesh, int eltid)
{
    mesh->nbranch(eltid) = 1;
}


void Inductor::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                           double cx, double cv, double ca)
{
    int id[NEN+1];
    set_local_id(mesh, eltid, id);
    QMatrix<dcomplex> Ke(NULL, NEN+1, NEN+1);
    dcomplex induct_c(induct, induct_i);
    if (cx) {
        Ke(0,2) +=  1.0 * cx;
        Ke(1,2) += -1.0 * cx;
        Ke(2,0) +=  1.0 * cx;
        Ke(2,1) += -1.0 * cx;
    }
    if (cv) {
        Ke(2,2) += -induct_c * cv;
    }
    K->add(id, NEN+1, Ke.data);
}


void Inductor::assemble_R(Mesh* mesh, int eltid)
{
    QMatrix<dcomplex> Fe(NULL, 1, NEN+1);
    dcomplex induct_c(induct, induct_i);
    dcomplex nodeu[NEN+1],nodev[NEN+1];
    set_local_u(mesh, eltid, nodeu);
    set_local_v(mesh, eltid, nodev);
    Fe(0) +=  nodeu[2];
    Fe(1) += -nodeu[2];
    Fe(2) += (nodeu[0] - nodeu[1] - induct_c * nodev[2]);
    add_local_f(mesh,eltid, Fe.data);

}


/* ==== Voltage or current source Model ==== */


VIsrc::VIsrc() :
    Element(1)
{
}


VIsrc::~VIsrc()
{
}


void VIsrc::initialize(Mesh* mesh, int eltid)
{
    mesh->nbranch(eltid) = 1;
}


void VIsrc::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                        double cx, double cv, double ca)
{
    int id[NEN+1];
    set_local_id(mesh, eltid, id);
    QMatrix<dcomplex> Ke(NULL, NEN+1, NEN+1);
    if (cx) {
        Ke(0,2) +=  1.0 * cx;
        Ke(1,2) += -1.0 * cx;
        Ke(2,0) +=  1.0 * cx;
        Ke(2,1) += -1.0 * cx;
    }
    K->add(id, NEN+1, Ke.data);
}


void VIsrc::assemble_R(Mesh* mesh, int eltid)
{
    QMatrix<dcomplex> Fe(NULL, 1, NEN+1);
    dcomplex nodeu[NEN+1];
    set_local_u(mesh, eltid, nodeu);
    Fe(0) +=  nodeu[2];
    Fe(1) += -nodeu[2];
    Fe(2) += (nodeu[0] - nodeu[1]);
    add_local_f(mesh,eltid, Fe.data);
}


/* ==== Plate Electrode Model ==== */


Electrode::Electrode(int vglobalid, double lt):
    Element(1), vglobalid(vglobalid), lt(lt)
{
}


Electrode::~Electrode()
{
}


void Electrode::initialize(Mesh* mesh, int eltid)
{
    mesh->nbranch(eltid) = 1;
}


void Electrode::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                            double cx, double cv, double ca)
{
    int id[3];
    id[0] = mesh->inode(id_slots[0],mesh->ix(0,eltid));
    id[1] = mesh->iglobal(vglobalid);
    id[2] = mesh->ibranch(0,eltid);

    QMatrix<dcomplex> Ke(NULL, NENE+2, NENE+2);
    if (cx) {
        Ke(1,2) +=  -1.0 * cx;
        Ke(2,0) +=   1.0 * cx;
        Ke(2,1) +=  -1.0 * cx;
    }
    if (cv) {
        Ke(0,2) +=   lt * 1.0 * cv;
    }
    K->add(id, NENE+2, Ke.data);

}

// FIXME
void Electrode::assemble_R(Mesh* mesh, int eltid)
{
    int id[3];
    id[0] = mesh->inode(id_slots[0],mesh->ix(0,eltid));
    id[1] = mesh->iglobal(vglobalid);
    id[2] = mesh->ibranch(0,eltid);

    QMatrix<dcomplex> Fe(NULL, 1, 3);
    dcomplex nodeu[3], nodev[3];
    for (int j = 0; j < 3; ++j) {
        nodeu[j] = dcomplex(mesh->u(id[j]), mesh->ui(id[j]));
        if (mesh->is_harmonic()) {
            dcomplex iw = dcomplex(0,1)*mesh->harmonic_freq();
            nodev[j] = iw*nodeu[j];
        } else
            nodeu[j] = mesh->v(id[j]);
    }

    Fe(0) =  lt * nodev[2];
    Fe(1) = -nodeu[2];
    Fe(2) = (nodeu[0] - nodeu[1]);

    for (int j = 0; j < 3; ++j) {
        mesh->f (id[j]) += real(Fe(j));
        mesh->fi(id[j]) += imag(Fe(j));
    }
}


/* ==== Plate Electrode2 Model ==== */


Electrode2::Electrode2(int vglobalid, double lt):
    Element(1), vglobalid(vglobalid), lt(lt)
{
}


Electrode2::~Electrode2()
{
}


void Electrode2::initialize(Mesh* mesh, int eltid)
{
    mesh->nbranch(eltid) = 2;
}


void Electrode2::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                             double cx, double cv, double ca)
{
    int id[4];
    id[0] = mesh->inode(id_slots[0],mesh->ix(0,eltid));
    id[1] = mesh->iglobal(vglobalid);
    id[2] = mesh->ibranch(0,eltid);
    id[3] = mesh->ibranch(1,eltid);

    QMatrix<dcomplex> Ke(NULL, 4,4);
    if (cx) {
        Ke(0,3) +=   1.0 * cx;
        Ke(1,2) +=  -1.0 * cx;
        Ke(2,0) +=   1.0 * cx;
        Ke(2,1) +=  -1.0 * cx;
        Ke(3,3) +=  -1.0 * cx;
    }
    if (cv) {
        Ke(3,2) +=   lt * 1.0 * cv;
    }
    K->add(id, 4, Ke.data);

}

// FIXME
void Electrode2::assemble_R(Mesh* mesh, int eltid)
{
    int id[4];
    id[0] = mesh->inode(id_slots[0],mesh->ix(0,eltid));
    id[1] = mesh->iglobal(vglobalid);
    id[2] = mesh->ibranch(0,eltid);
    id[3] = mesh->ibranch(1,eltid);

    QMatrix<dcomplex> Fe(NULL, 1, 4);
    dcomplex nodeu[4],nodev[4];
    for (int j = 0; j < 4; ++j) {
        nodeu[j] = dcomplex(mesh->u(id[j]), mesh->ui(id[j]));
        if (mesh->is_harmonic()) {
            dcomplex iw = dcomplex(0,1)*mesh->harmonic_freq();
            nodev[j] = iw*nodeu[j];
        } else {
            nodev[j] = mesh->v(id[j]);
        }
    }

    Fe(0) =  nodeu[3];
    Fe(1) = -nodeu[2];
    Fe(2) = (nodeu[0] - nodeu[1]);
    Fe(3) = (lt * nodev[2] - nodeu[3]);

    for (int j = 0; j < 4; ++j) {
        mesh->f (id[j]) += real(Fe(j));
        mesh->fi(id[j]) += imag(Fe(j));
    }
}
