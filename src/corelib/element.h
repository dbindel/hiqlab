/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include "qassembly.h"
#include "fieldeval.h"

class Mesh;

/*@T --------------------------------------------------------------
 * \section{Element interface}
 *
 * The [[Element]] class in HiQLab really refers to an element
 * type.  Individual elements in the mesh are not represented by
 * separate objects; rather, we use a ``flyweight'' pattern in which
 * characteristics of generic element types are stored inside an
 * element object, and characteristics of particular instances are
 * stored externally (in the mesh object).  This external information
 * is passed to the element methods as a pointer to the mesh and an
 * element identifier.  In order to do element-by-element
 * initialization and assembly tasks, then, we loop over each element
 * in the mesh and make calls of the form
 * [[mesh->etypes[eltid]->method(mesh, eltid, ...)]].
 *
 *@T --------------------------------------------------------------
 * \subsection{Variable slots}
 *
 * Assembling element contributions into a global system requires that
 * all elements share a common idea of what nodal variable number is
 * assigned to what type of variable.  The variable slot mechanism
 * provides a way to reconcile the global assignment of nodal variable
 * numbers with a local assignment specific to each element type.
 * 
 * For example, suppose we wanted to perform a simulation that
 * combined mechanical, electrical, and thermal calculations.  We
 * might decide at the global level that variables 0-2 at each node
 * correspond to displacements, variable 3 corresponds to
 * electrostatic potential, and variable 4 corresponds to temperature.
 * But a thermomechanical element would not know about electrostatic
 * potential; from the perspective of that element, perhaps variables
 * 0-2 should correspond to displacement and variable 3 should
 * correspond to temperature.  In this case, the slot array for the
 * coupling element would look like [[{0, 1, 2, 4}]].
 * 
 *@q*/

/*@T --------------------------------------------------------------
 * \subsection{Element declarations}
 *
 *@c*/

class Element {
public:
    Element(int nslots);
    virtual ~Element();

    /** Initialize element eltid of the mesh.  Allocate element vars. */
    virtual void initialize(Mesh* mesh, int eltid);

    /** Mark element IDs as active */
    virtual void assign_ids(Mesh* mesh, int eltid);

    /** Assemble structure */
    virtual void assemble_struct(Mesh* mesh, int eltid, QStructAssembler* K);

    /** Assemble element matrices */
    virtual void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                             double cx=1, double cv=0, double ca=0);

    /** Assemble element residual */
    virtual void assemble_R(Mesh* mesh, int eltid);

    // FIXME
    /** Apply tangent operator */
    virtual void apply_dR(Mesh* mesh, int eltid, 
                          double cx=1, double cv=0, double ca=0);

    // FIXME
    /** Apply tangent operator in time-harmonic case */
    virtual void apply_harmonic(Mesh* mesh, int eltid, dcomplex w);

    /** Compute the (lumped) L2 projection of a field defined
     *  at Gauss nodes.  Output is Mdiag += sum(integral N*N') and
     *  xfields += integral N*xfunc.
     */
    virtual void project_L2(Mesh* mesh, int eltid, int nfields,
                            double* Mdiag, double* xfields,
                            FieldEval& xfunc);

    /** Compute the stresses for an element at parent coordinates X.
     *  Returns an incremented stress pointer.  Use two calls to
     *  compute the stress associated with a complex displacement.
     */
    virtual double* stress(Mesh* mesh, int eltid, double* X, double* stress);

    /** Compute the time-averaged energy flux for an element at parent
     *  coordinates X.  Returns an incremented pointer.
     */
    virtual double* mean_power(Mesh* mesh, int eltid, double* X, double* EX);

    /** Read and write the id_slot map
     */
    int  num_id_slots()        { return id_slots.size();  }
    int  id_slot(int i) const  { return id_slots[i]; }
    int& id_slot(int i)        { return id_slots[i]; }
    void id_slot(int i, int x) { id_slots[i] = x;   }
    int  max_id_slot();

protected:
    std::vector<int> id_slots;    // id_slots[i] = slot for ith local var

    /** Default initializer -- just marks the mesh ID array. */
    void assign_ids_default(Mesh* mesh, int elt);

    void set_local_arrays(Mesh* mesh, int eltid, int* id, double* nx, int ndm);
    void set_local_id(Mesh* mesh, int eltid, int* id);
    static void set_local_x(Mesh* mesh, int eltid, double* nx, int ndm);

    void set_local_vec(Mesh* mesh, int v,         int eltid, double*   nu);
    void set_local_vec(Mesh* mesh, int v, int vi, int eltid, dcomplex* nu);
    void set_local_uvec(Mesh* mesh, int eltid, dcomplex c, dcomplex* nu);

    void set_local_u(Mesh* mesh, int eltid, dcomplex* nu);
    void set_local_v(Mesh* mesh, int eltid, dcomplex* nv);
    void set_local_a(Mesh* mesh, int eltid, dcomplex* na);
    void add_local_f(Mesh* mesh, int eltid, dcomplex* nodef1);

    void set_local_u(Mesh* mesh, int eltid, double* nu);
    void set_local_v(Mesh* mesh, int eltid, double* nv);
    void set_local_a(Mesh* mesh, int eltid, double* na);
    void add_local_f(Mesh* mesh, int eltid, double* nodef1);
};

/*@q*/
#endif /* ELEMENT_H */
