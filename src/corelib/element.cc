/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstdio>

#include "qmatrix.h"
#include "element.h"
#include "mesh.h"

#define MAXNEN       125


/*@T --------------------------------------------------------------
 * \subsection{Element construction and destruction}
 *
 * The only internal data structure common to all element types is the
 * variable slot map.  By default, we initialize this map to an
 * appropriately-sized identity.  We defer to a higher-level routine
 * (typically in the Lua code) the coordination of how slots should be
 * assigned in a multiphysics problem
 *
 * Because this is a base class, we need a virtual destructor.  Said
 * destructor doesn't need to do anything, though.
 *
 *@c*/

#define ME Element


ME::ME(int nslots) :
    id_slots(nslots)
{
    for (int i = 0; i < nslots; ++i)
        id_slots[i] = i;
}


ME::~ME()
{
}


/*@T --------------------------------------------------------------
 * \subsection{Element initialization}
 *
 * The element initialization phase is where the element tells
 * HiQLab how many branch variables and history variables it will
 * need.  It does this by setting [[mesh->nbranch(eltid)]] and
 * [[mesh->nhist(eltid)]].  An element without branch or history
 * variables can use the default (empty) initialization routine.
 *
 *@c*/

void ME::initialize(Mesh* mesh, int eltid)
{
}


/*@T --------------------------------------------------------------
 * \subsection{Allocating variables}
 *
 * The [[assign_ids]] method marks the nodal and branch variables
 * that it will use.  By default, we assume that this includes all the
 * nodal variables listed in the slots array at each node attached to
 * the element, plus any branch variables that were allocated.  Some
 * element types might not use all variables at all nodes, though; for
 * example, a mixed pressure-displacement formulation could have
 * pressure only associated with an internal node.
 *
 *@c*/

void ME::assign_ids(Mesh* mesh, int eltid)
{
    assign_ids_default(mesh, eltid);
}


void ME::assign_ids_default(Mesh* mesh, int eltid)
{
    int nen = mesh->get_nen(eltid);
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        for (int i = 0; i < num_id_slots(); ++i)
            mesh->id(id_slots[i], nodeid) = 1;
    }
    int nbranch = mesh->nbranch_id(eltid);
    for (int j = 0; j < nbranch; ++j)
        mesh->branchid(j, eltid) = 1;
}


/*@T --------------------------------------------------------------
 * \subsection{Assembling the residual and tangent}
 *
 * The tangent stiffness looks like
 * \[
 *   K = c_u \frac{\partial R}{\partial u} +
 *       c_v \frac{\partial R}{\partial v} +
 *       c_a \frac{\partial R}{\partial a}
 * \]
 * where $u$, $v$, and $a$ are the displacement, velocity, and
 * acceleration vectors, respectively.  The [[assemble_dR]] method
 * adds this element's contribution to $K$.  The
 * [[assemble_struct]] method just assembles the nonzero structure
 * of the tangent stiffness.
 *
 * The residual vector $R(u,v,a)$ is assembled into the $F$ vector in
 * the mesh object.  [[assemble_R]] adds this element's
 * contribution.
 * 
 *@c*/


void ME::assemble_struct(Mesh* mesh, int eltid, QStructAssembler* K)
{
    int nen = mesh->get_nen(eltid);
    int nslots = num_id_slots();
    int nbranch = mesh->nbranch_id(eltid);
    std::vector<int> id(nen*nslots+nbranch);
    set_local_id(mesh, eltid, &(id[0]));
    K->add(&(id[0]), nen*nslots+nbranch);
}


void ME::assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx, double cv, double ca)
{
}


void ME::assemble_R(Mesh* mesh, int eltid)
{
}


void ME::apply_dR(Mesh* mesh, int eltid, double cx, double cv, double ca)
{
}


void ME::apply_harmonic(Mesh* mesh, int eltid, dcomplex w)
{
}


/*@T --------------------------------------------------------------
 * \subsection{Lumped $L^2$ projections}
 *
 * Computing a lumped $L^2$ projection of a quantity defined at the
 * Gauss points is a two-step procedure.  First, call the element
 * [[project_L2]] method on each element in the mesh.  The output
 * of this calculation should be
 * \begin{eqnarray*}
 *   [[Mdiag(i)]] & = & \int_{\Omega} N_i \, d\Omega \\
 *   [[xfields(i,j)]] & = & \int_{\Omega} N_i(x) v_j(x) \, d\Omega \\
 * \end{eqnarray*}
 * where $v_j(x)$ is the $j$th scalar component of the fields to be
 * evaluated, and $N_i(x)$ is the shape associated with a unit displcement
 * of variable $i$.
 *
 * At the end of the initial loop, we scale [[xfields(i,j)]] by
 * [[Mdiag(i)]] in order to get the approximate $L^2$ projection of
 * the fields onto the nodal shape functions.
 *@c*/


void ME::project_L2(Mesh* mesh, int eltid, int nfields,
                    double* Mdiag, double* xfields,
                    FieldEval& func)
{
}


/*@T --------------------------------------------------------------
 * \subsection{Gauss point quantities}
 *
 * Right now, there are two types of quantities that we can compute in
 * an element (typically at the Gauss points).  The first is the
 * stress components; the second is the mean energy flux in a
 * time-harmonic simulation.
 *
 * The latter does seem fairly specific, and it's not coded for all
 * elements.  This piece of the infrastructure should probably be
 * handled more cleanly.
 *
 *@c*/


double* ME::stress(Mesh* mesh, int eltid, double* X, double* stress)
{
    return stress;
}


double* ME::mean_power(Mesh* mesh, int eltid, double* X, double* EX)
{
    return EX;
}


/*@T --------------------------------------------------------------
 * \subsection{Counting the number of ID slots in use}
 *
 * The [[max_id_slot]] routine returns the last (global) variable
 * index listed in the slot array.  Taking the maximum of the
 * [[max_id_slot]] return values over all elements tells us the
 * maximum number of variables we need to allocate per node.
 * I cannot think of any good reason that a user should call this
 * routine -- it only seems useful inside the mesh [[initialize]]
 * method.
 *
 *@c*/


int ME::max_id_slot()
{
    int max_slot = 0;
    for (int i = 0; i < num_id_slots(); ++i)
        if (max_slot < id_slots[i])
            max_slot = id_slots[i];
    return max_slot;
}


/*@T --------------------------------------------------------------
 * \subsection{Transferring data from global arrays}
 *
 * For element calculations, we typically will want slices of the
 * global [[ID]], [[X]], [[U]], [[V]], and [[A]]
 * arrays.  That's what the [[set_local_*]] routines do.  Of these,
 * the only one that is less than obvious is the [[set_local_id]]
 * function, which returns in its output array the list of all nodal
 * variables for the element, followed by the list of all branch
 * variables.
 *
 *@c*/


void ME::set_local_arrays(Mesh* mesh, int eltid,
                          int* id, double* nodex, int ndm)
{
    set_local_x (mesh, eltid, nodex, ndm);
    set_local_id(mesh, eltid, id);
}


void ME::set_local_id(Mesh* mesh, int eltid, int* id1)
{
    int nen = mesh->get_nen(eltid);
    QMatrix<int> id(id1, num_id_slots(), nen);
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        for (int i = 0; i < num_id_slots(); ++i)
            id(i,j) = mesh->inode(id_slots[i], nodeid);
    }
    int nbranch = mesh->nbranch_id(eltid);
    for (int i = 0; i < nbranch; ++i)
        id1[nen*num_id_slots()+i] = mesh->ibranch(i,eltid);
}


void ME::set_local_x(Mesh* mesh, int eltid, double* nodex1, int ndm)
{
    int nen = mesh->get_nen(eltid);
    QMatrix<double> nodex(nodex1, ndm,nen);
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        for (int i = 0; i < ndm; ++i)
            nodex(i,j) = mesh->v(i,nodeid);
    }
}


void ME::set_local_vec(Mesh* mesh, int v, int eltid, double* nodex1)
{
    int nen = mesh->get_nen(eltid);
    QMatrix<double> nodex(nodex1, num_id_slots(), nen);
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        for (int i = 0; i < num_id_slots(); ++i)
            nodex(i,j) = mesh->vec(v,id_slots[i],nodeid);
    }
    int nbranch = mesh->nbranch_id(eltid);
    for (int i = 0; i < nbranch; ++i)
        nodex1[nen*num_id_slots()+i] = mesh->branchvec(v,i,eltid);
}


void ME::set_local_vec(Mesh* mesh, int v, int vi, int eltid, dcomplex* nodex1)
{
    int nen = mesh->get_nen(eltid);
    QMatrix<dcomplex> nodex(nodex1, num_id_slots(), nen);
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        for (int i = 0; i < num_id_slots(); ++i) {
            dcomplex z(mesh->vec(v, id_slots[i],nodeid),
                       mesh->vec(vi,id_slots[i],nodeid));
            nodex(i,j) = z;
        }
    }
    int nbranch = mesh->nbranch_id(eltid);
    for (int i = 0; i < nbranch; ++i) {
        dcomplex z(mesh->branchvec(v, i,eltid),
                   mesh->branchvec(vi,i,eltid));
        nodex1[nen*num_id_slots()+i] = z;
    }
}


void ME::set_local_uvec(Mesh* mesh, int eltid, dcomplex c, dcomplex* nodex1)
{
    int nen = mesh->get_nen(eltid);
    QMatrix<dcomplex> nodex(nodex1, num_id_slots(), nen);
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        for (int i = 0; i < num_id_slots(); ++i) {
            dcomplex z(mesh->vec(Mesh::VEC_U, id_slots[i],nodeid),
                       mesh->vec(Mesh::VEC_UI,id_slots[i],nodeid));
            nodex(i,j) = c*z;
        }
    }
    int nbranch = mesh->nbranch_id(eltid);
    for (int i = 0; i < nbranch; ++i) {
        dcomplex z(mesh->branchvec(Mesh::VEC_U, i,eltid),
                   mesh->branchvec(Mesh::VEC_UI,i,eltid));
        nodex1[nen*num_id_slots()+i] = c*z;
    }
}


void ME::set_local_u(Mesh* mesh, int eltid, double* nodeu)
{
    set_local_vec(mesh, Mesh::VEC_U, eltid, nodeu);
}

void ME::set_local_v(Mesh* mesh, int eltid, double* nodeu)
{
    set_local_vec(mesh, Mesh::VEC_V, eltid, nodeu);
}

void ME::set_local_a(Mesh* mesh, int eltid, double* nodeu)
{
    set_local_vec(mesh, Mesh::VEC_A, eltid, nodeu);
}

void ME::set_local_u(Mesh* mesh, int eltid, dcomplex* nodeu)
{
    // set_local_vec(mesh, Mesh::VEC_U, Mesh::VEC_UI, eltid, nodeu);
    set_local_uvec(mesh, eltid, 1.0, nodeu);
}

void ME::set_local_v(Mesh* mesh, int eltid, dcomplex* nodeu)
{
    // set_local_vec(mesh, Mesh::VEC_V, Mesh::VEC_VI, eltid, nodeu);
    dcomplex w0 = 0.0;
    if (mesh->is_harmonic())
        w0 = mesh->harmonic_freq();
    w0 *= dcomplex(0,1);
    set_local_uvec(mesh, eltid, w0, nodeu);
}

void ME::set_local_a(Mesh* mesh, int eltid, dcomplex* nodeu)
{
    // set_local_vec(mesh, Mesh::VEC_A, Mesh::VEC_AI, eltid, nodeu);
    dcomplex w0 = 0.0;
    if (mesh->is_harmonic())
        w0 = mesh->harmonic_freq();
    w0 *= dcomplex(0,1);
    set_local_uvec(mesh, eltid, w0*w0, nodeu);
}


/*@T --------------------------------------------------------------
 * \subsection{Assembling data from global arrays}
 *
 * The output of some of the assembly loops goes to external storage,
 * but the [[assemble_R]] functions assemble their results into the
 * mesh [[F]] array.  The [[add_local_f]] functions take the
 * local element residual contribution and add it into this array in
 * the default way.
 *
 *@c*/


void ME::add_local_f(Mesh* mesh, int eltid, double* nodef1)
{
    int nen = mesh->get_nen(eltid);
    QMatrix<double> nodef(nodef1, num_id_slots(),nen);
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        for (int i = 0; i < num_id_slots(); ++i)
            mesh->f(id_slots[i],nodeid) += nodef(i,j);
    }
    int nbranch = mesh->nbranch_id(eltid);
    for (int i = 0; i < nbranch; ++i)
        mesh->f(mesh->ibranch(i,eltid)) += nodef1[nen*num_id_slots()+i];
}


void ME::add_local_f(Mesh* mesh, int eltid, dcomplex* nodef1)
{
    int nen = mesh->get_nen(eltid);
    QMatrix<dcomplex> nodef(nodef1, num_id_slots(),nen);
    for (int j = 0; j < nen; ++j) {
        int nodeid = mesh->ix(j,eltid);
        for (int i = 0; i < num_id_slots(); ++i) {
            mesh->f (id_slots[i],nodeid) += real(nodef(i,j));
            mesh->fi(id_slots[i],nodeid) += imag(nodef(i,j));
        }
    }
    int nbranch = mesh->nbranch_id(eltid);
    for (int i = 0; i < nbranch; ++i) {
        mesh->f (mesh->ibranch(i,eltid)) += real(nodef1[nen*num_id_slots()+i]);
        mesh->fi(mesh->ibranch(i,eltid)) += imag(nodef1[nen*num_id_slots()+i]);
    }
}
