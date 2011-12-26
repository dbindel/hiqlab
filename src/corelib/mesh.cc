/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cassert>
#include <cstdio>
#include <algorithm>
#include <iostream>

#include "luasupport.h"
#include "tolua++.h"
#include "mesh.h"
#include "coordsorter.h"

using std::vector;
using std::fill;


/* ----------
 * Local utility functions and macros
 */


#define CHECK_LUA  if (!L) return


template<class T>
inline void clear(vector<T>& U)
{
    fill(U.begin(), U.end(), 0);
}


template<class T>
inline void clear(vector<T>& U, int M)
{
    U.resize(M);
    clear(U);
}


/*@T ----------
 * \subsection{Mesh construction, destruction, and memory management}
 *
 * Most of the data in the mesh object lives in C++ containers, which
 * are automatically managed.  An exception to this is the element types
 * and the Lua state.  We require that the user explicitly hand over
 * ownership of these types for two reasons.  First, it may not always
 * make sense for the mesh object to own these objects -- for example,
 * when using the Lua front-end, the Lua interpreter will typically outlive
 * the mesh!  Second, the mesh object isn't responsible for constructing
 * these objects, so unless there is some explicit registration mechanism,
 * the mesh will not necessarily know that the objects exist.
 *
 * It is a checked error to try to unassign or reassign ownership of
 * the Lua state after the first assignment.
 *
 *@c*/


Mesh::Mesh(int ndm) :
    ndm(ndm), maxnen(0), maxndf(0), numid(0), which_hist(0),
    lua_owned(NULL), L(NULL)
{
    NG[0] = 0;
    NG[1] = 0;
    NIX.push_back(0);
}


Mesh::~Mesh()
{
    for (unsigned i = 0; i < etypes_owned.size(); ++i)
        if (etypes_owned[i])
            delete etypes_owned[i];
    if (lua_owned)
        lua_close(lua_owned);
}


Element* Mesh::own(Element* e)
{
    etypes_owned.push_back(e);
    return e;
}


void Mesh::own(lua_State* L)
{
    assert(lua_owned == NULL);
    assert(L != NULL);
    lua_owned = L;
}


/*@T ----------
 * \subsection{Mesh scale management}
 *
 * We keep two distinct sets of scaling information associated with a
 * mesh.  
 *
 * First, there is the abstract scaling information associated with
 * different variable types (temperatures, displacements, etc).  The
 * user provides access to these scales by associating a
 * [[get_scale]] method with the Lua mesh object.  This method maps
 * a string argument to a numeric scaling value.  If the method isn't
 * there, or doesn't return anything, we return a default scale value
 * of 1.
 *
 * Second, there are scales associated with each primary and dual
 * variable in the system (including element and global variables).
 * These scales are stored in the [[d1]] and [[d2]] arrays,
 * respectively.  The default scales are all ones, but the defaults
 * would typically be overridden by the [[set_scaling_vector]]
 * method.  This method in turn calls [[set_nodal_u_scale]] and
 * [[set_nodal_f_scale]] to set all the scales associated with a
 * particular slot in the nodal variable array (since a given slot is
 * associated with a particular field -- see the element interface
 * description).
 *
 *@c*/


double Mesh::get_scale(const char* name)
{
    CHECK_LUA 1;
    lua_pushstring(L, name);
    double retval = 1;
    if (lua_method("get_scale", 1, 1) == 0) {
        retval = lua_tonumber(L,-1);
        lua_pop(L,1);
    }
    return retval;
}


void Mesh::clear_scaling_vector()
{
    int N = ID.size();
    for (int i = 0; i < N; ++i) {
        d1(i) = 1;
        d2(i) = 1;
    }
}


void Mesh::set_nodal_u_scale(int i, double s)
{
    for (int j = 0; j < numnp(); ++j)
        d1(i,j) = s;
}


void Mesh::set_nodal_f_scale(int i, double s)
{
    for (int j = 0; j < numnp(); ++j)
        d2(i,j) = s;
}


void Mesh::set_scaling_vector()
{
    CHECK_LUA;
    lua_method("set_scaling_vector", 0, 0);
}


void Mesh::get_scaling_vector(double* su, double* sf, int is_reduced)
{
    int N = (is_reduced ? ID.size() : numnp()*maxndf );
    dvec& D1 = vecs[VEC_D1];
    dvec& D2 = vecs[VEC_D2];
    for (unsigned i = 0; i < (unsigned) N; ++i) {
        int k = (is_reduced ? id(i) : i);
        if (k >= 0) {
            if (su) su[k] = D1[i];
            if (sf) sf[k] = D2[i];
        }
    }
}


/*@T ----------
 * \subsection{Basic mesh generation}
 *
 *@c*/


int Mesh::add_node(double* x, int n)
{
    int id = X.size() / ndm;
    for (int i = 0; i < n*ndm; ++i)
        X.push_back(x[i]);
    return id;
}


int Mesh::add_global(int n)
{
    int id = NG[1];
    NG[1] += n;
    return id;
}


int Mesh::add_element(int* e, Element* etype, int nen, int n)
{
    int id = etypes.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nen; ++j)
            IX.push_back(e[i*nen+j]);
        end_add_element(etype);
    }
    return id;
}


void Mesh::end_add_element(Element* etype)
{
    NIX.push_back(IX.size());
    etypes.push_back(etype);
    NB.push_back(0);
    NH.push_back(0);
    if (etype)
        etype->initialize(this, etypes.size()-1);
}


/*@T ----------
 * \subsection{Block mesh generation}
 *
 * The block mesh generators produce rectilinear coordinate-aligned
 * block meshes, which can subsequently be turned into more
 * interesting shape by mapping and tying operations.  For each
 * dimension, we specify the coordinate range ($[x_1, x_2]$) and the
 * number of {\em nodes} (not elements) in that direction.  It is a
 * checked error to try to specify a number of nodes that doesn't
 * yield an integer number of elements.
 *
 **@c*/


void Mesh::add_block(double* x1, double* x2, int* m,
                     Element* etype, int order)
{
    if (ndm == 1)
        add_block(x1[0], x1[0], m[0], etype, order);

    else if (ndm == 2)
        add_block(x1[0], x1[1],
                  x2[0], x2[1],
                  m[0],  m[1], etype, order);

    else if (ndm == 3)
        add_block(x1[0], x1[1], x1[2],
                  x2[0], x2[1], x2[2],
                  m[0],  m[1],  m[2], etype, order);
}


void Mesh::add_block(double x1, double x2, int nx,
                     Element* etype, int order)
{
    assert( nx > 0 && (nx-1)%order == 0 );

    int start_node = numnp();
    for (int ix = 0; ix < nx; ++ix) {
        double x[1];
        x[0] = ((nx-1-ix)*x1 + ix*x2) / (nx-1);
        add_node(x);
    }

    for (int ix = 0; ix < nx-order; ix += order) {
        for (int jx = 0; jx < order+1; ++jx)
            IX.push_back( (ix+jx) + start_node );
        end_add_element(etype);
    }
}


void Mesh::add_block(double x1, double y1,
                     double x2, double y2,
                     int nx, int ny,
                     Element* etype, int order)
{
    assert( nx > 0 && (nx-1)%order == 0 );
    assert( ny > 0 && (ny-1)%order == 0 );

    int start_node = numnp();
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            double x[2];
            x[0] = ((nx-1-ix)*x1 + ix*x2) / (nx-1);
            x[1] = ((ny-1-iy)*y1 + iy*y2) / (ny-1);
            add_node(x);
        }
    }

    for (int ix = 0; ix < nx-order; ix += order)
        for (int iy = 0; iy < ny-order; iy += order) {
            for (int jx = 0; jx < order+1; ++jx)
                for (int jy = 0; jy < order+1; ++jy)
                    IX.push_back( (iy+jy) + ny*(ix+jx) + start_node );
            end_add_element(etype);
        }
}


void Mesh::add_block(double x1, double y1, double z1,
                     double x2, double y2, double z2,
                     int nx, int ny, int nz,
                     Element* etype, int order)
{
    assert( nx > 0 && (nx-1)%order == 0 );
    assert( ny > 0 && (ny-1)%order == 0 );
    assert( nz > 0 && (nz-1)%order == 0 );

    int start_node = numnp();

    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                double x[3];
                x[0] = ((nx-1-ix)*x1 + ix*x2) / (nx-1);
                x[1] = ((ny-1-iy)*y1 + iy*y2) / (ny-1);
                x[2] = ((nz-1-iz)*z1 + iz*z2) / (nz-1);
                add_node(x);
            }
        }
    }

    for (int ix = 0; ix < nx-order; ix += order)
        for (int iy = 0; iy < ny-order; iy += order)
            for (int iz = 0; iz < nz-order; iz += order) {
                for (int jx = 0; jx < order+1; ++jx)
                    for (int jy = 0; jy < order+1; ++jy)
                        for (int jz = 0; jz < order+1; ++jz)
                            IX.push_back( (iz+jz) + nz*((iy+jy) + ny*(ix+jx)) +
                                          start_node );
                end_add_element(etype);
            }
}


/*@T ----------
 * \subsection{Tied meshes}
 *
 * We tie meshes together by approximately sorting the coordinates,
 * looking first at $x$, then $y$, then $z$.  Nodes that are
 * equivalent under this approximate sort (and therefore appear
 * together in the sorted node list) are tied.  For each equivalence
 * class, we pick one representative node, and map all references to
 * other nodes in the class into references to that representative.
 * 
 * The assumption in this algorithm is that all nodes that are
 * supposed to be tied together are within a given tolerance of each
 * other (and all nodes that are not supposed to be tied are farther
 * than that tolerance).  It's possible to build a pathological mesh
 * where this assumption fails; for example, if we have nodes at
 * $(0,0)$, $(0,\epsilon)$, $(0,-\epsilon)$ where $\epsilon$ is the
 * tolerance, then the relation we've defined is not transitive, and
 * we don't have an equivalence relation.  In this case, things that
 * should probably be tied together might not be.  I regard this as
 * an unchecked blunder on the part of the analyst.
 *
 **@c*/

void Mesh::tie(double tol, int start, int end)
{
    int np = numnp();
    if (start < 0) start = 0;
    if (end   < 0) end   = np;
    int nt = end-start;

    // Set scratch = start:end-1, then sort the indices according
    // to the order of the corresponding nodes.
    //
    CoordSorter sorter(ndm, *this, tol);
    vector<int> scratch(nt);
    for (int i = 0; i < nt; ++i) 
        scratch[i] = start+i;
    sort(scratch.begin(), scratch.end(), sorter);

    // map[i] := smallest index belonging to equivalence class of p[i]
    //
    vector<int> map(np);
    for (int i = 0; i < np; ++i) 
        map[i] = i;

    int inext;
    for (int i = 0; i < nt; i = inext) {
        int class_id = scratch[i];
        for (inext = i+1; inext < nt; ++inext) {
             if (sorter.compare(scratch[i], scratch[inext]) != 0)
                 break;
             if (scratch[inext] < scratch[i])
                class_id = scratch[inext];
        }
        for (int j = i; j < inext; ++j)
            map[scratch[j]] = class_id;
    }

    // Apply map to element connectivity array
    //
    for (int i = 0; i < (int) IX.size(); ++i)
        if (IX[i] >= 0)
            IX[i] = map[IX[i]];
}


/*@T ----------
 * \subsection{Initialization}
 *
 * Mesh initialization occurs after the geometry of the mesh has been
 * established and before any analysis takes place.  This is where we
 * decide how big all the major arrays will be, and set up the index
 * structures that are used for all of our assembly operations.
 *
 **@c*/


void Mesh::initialize()
{
    build_reverse_map();              // Build node-to-element map
    initialize_sizes();               // Figure out maxnen and maxndf
    int M  = initialize_indexing();   // Assign branch and global dofs
    int MH = initialize_history();    // Assign history storage
    initialize_arrays(M, MH);         // Allocate storage arrays

    assign_ids();
    apply_bc();
    set_scaling_vector();
}


void Mesh::initialize_sizes()
{
    for (int i = 0; i < etypes.size(); ++i) {
        int nen = get_nen(i);
        if (maxnen < nen)
            maxnen = nen;
    }
    for (int i = 0; i < etypes.size(); ++i) {
        if (etypes[i]) {
            int elt_maxndf = etypes[i]->max_id_slot() + 1;
            if (maxndf < elt_maxndf)
                maxndf = elt_maxndf;
        }
    }
}


void Mesh::build_reverse_map()
{
    int N = numnp();
    int M = etypes.size();
    clear(NIN,N+1);

    // -- NIN[i] := elements for node i
    for (unsigned i = 0; i < IX.size(); ++i)
        if (IX[i] >= 0)
            NIN[IX[i]]++;

    // -- NIN[i] := elements for nodes 0 .. i
    for (int i = 1; i < N+1; ++i)
        NIN[i] += NIN[i-1];

    // -- Build reverse map via bucket sort
    clear(IN, NIN[N]);
    for (int i = 0; i < M; ++i) {
        for (int j = NIX[i]; j < NIX[i+1]; ++j) {
            int node = IX[j];
            int slot = --NIN[node];
            IN[slot] = i;
        }
    }
}


int Mesh::initialize_indexing()
{
    // -- Assign index ranges for branch dofs
    int M = numnp()*maxndf;
    for (int i = 0; i < etypes.size(); ++i) {
        int start_M = M;
        M += NB[i];
        NB[i] = start_M;
    }
    NB.push_back(M);

    // -- Assign index ranges for global dofs
    NG[0] += M;
    NG[1] += M;
    M = NG[1];

    return M;
}


int Mesh::initialize_history()
{
    int MH = 0;
    for (int i = 0; i < etypes.size(); ++i) {
        int start_MH = MH;
        MH += NH[i];
        NH[i] = start_MH;
    }
    NH.push_back(MH);
    return MH;
}


void Mesh::initialize_arrays(int M, int MH)
{
    clear(vecs[VEC_U], M); clear(vecs[VEC_UI], M);
    clear(vecs[VEC_V], M); 
    clear(vecs[VEC_A], M); 
    clear(vecs[VEC_F], M); clear(vecs[VEC_FI], M);
    clear(vecs[VEC_BV],M); clear(vecs[VEC_BVI],M);
    clear(vecs[VEC_D1],M);
    clear(vecs[VEC_D2],M);
    clear(BC,M);
    clear(ID,M);
    clear(HIST, 2*MH);
}


/*@T ----------
 * \subsection{Assigning reduced indices}
 *
 * The [[assign_ids]] method sets up the map from the full index
 * set to the reduced index set.  Only the active variables (variables
 * that someone cares about that are not subject to a displacement BC)
 * are assigned valid reduced indices.  This method is called at
 * initialization, but it can also be called at subsequent phases in
 * the analysis if the user wishes to change the boundary conditions.
 *
 **@c*/


int Mesh::assign_ids()
{
    int N = etypes.size();   // Number of elements
    int M = ID.size();       // Number of potential dofs

    // -- Mark active dofs
    clear(ID);
    for (int i = 0; i < N; ++i) {
        if (etypes[i])
            etypes[i]->assign_ids(this, i);
    }
    for (int i = NG[0]; i < NG[1]; ++i)
        ID[i] = 1;

    // -- Assign IDs to nodal dofs
    numid = 0;
    for (int i = 0; i < M; ++i) {
        if (ID[i] != 0 && BC[i] != 'u')
            ID[i] = numid++;
        else
            ID[i] = -1;
    }

    return numid;
}


/*@T ----------
 * \subsection{Assembly loops}
 *
 * The [[assemble_dR]] routine assembles a linear combination of
 * mass, damping, and stiffness matrices.  The [[assemble_struct]]
 * routine is used to accumulate the nonzero structure of the matrix
 * for pre-allocation. 
 *
 * The [[assemble_R]] routine assembles the residual.  It's a
 * little quirky because it not only assembles the residual
 * contributions from the elements, but it also computes the
 * components of the extended residual that are associated with global
 * shape functions.  That is, we compute something like
 * \[
 *   R_g = G^T R_0
 * \]
 * where $R_g$ is the residual associated with the global shapes
 * described by the columns of $G$, and $R_0$ is the unreduced
 * residual without the global shape functions.
 *
 * The [[mean_power]] loop assembles a lumped $L^2$ projection
 * of the mean power as computed at each Gauss point.  It's a sort of
 * quirky routine -- why just the mean power, and not also other Gauss
 * point quantities like the stress?  We should probably rethink this
 * routine t some point.
 *
 **@c*/


void Mesh::assemble_struct(QStructAssembler* K, int reduced)
{
    QGlobalStructAssembler gassembler(this, *K);
    QReduceStructAssembler rassembler(this, gassembler);
    QStructAssembler* assembler = 
        (reduced ? 
         static_cast<QStructAssembler*>(&rassembler): 
         static_cast<QStructAssembler*>(&gassembler));

    assemble_struct_raw(assembler);
}


void Mesh::assemble_struct_raw(QStructAssembler* K)
{
    for (unsigned i = 0; i < etypes.size(); ++i) {
        if (etypes[i])
            etypes[i]->assemble_struct(this, i, K);
    }
}


void Mesh::assemble_dR(QAssembler* K, double cx, double cv, double ca,
                       int reduced)
{
    QGlobalAssembler gassembler(this, *K);
    QReduceAssembler rassembler(this, gassembler);
    QAssembler* assembler = 
        (reduced ? 
         static_cast<QAssembler*>(&rassembler) : 
         static_cast<QAssembler*>(&gassembler));

    assemble_dR_raw(assembler,cx,cv,ca);
}


void Mesh::assemble_dR_raw(QAssembler* K, double cx, double cv, double ca)
{
    for (unsigned i = 0; i < etypes.size(); ++i)
        if (etypes[i])
            etypes[i]->assemble_dR(this, i, K, cx, cv, ca);
}


void Mesh::assemble_R()
{
    set_f (NULL);
    set_fi(NULL);

    for (unsigned i = 0; i < etypes.size(); ++i)
        if (etypes[i])
            etypes[i]->assemble_R(this, i);

    dvec& F  = vecs[VEC_F];
    dvec& Fi = vecs[VEC_FI];
    if (numglobals() > 0) {
        for (int j = 0; j < numglobals(); ++j) {
            vector<double> scratch;
            scratch.resize(ID.size());
            shapeg(&scratch[0], 1, j, 0);
            double fgjr = 0;
            double fgji = 0;
            for (unsigned k = 0; k < ID.size(); ++k) {
                fgjr += F [k]*scratch[k];
                fgji += Fi[k]*scratch[k];
            }
            globalf (j) += fgjr;
            globalfi(j) += fgji;
        }
    }
}


void Mesh::mean_power(double* E)
{
    FieldEvalPower mean_power_func;
    vector<double> Mdiag(numnp());

    fill(Mdiag.begin(), Mdiag.end(),   0);
    fill(E,             E+numnp()*ndm, 0);

    for (unsigned i = 0; i < etypes.size(); ++i)
        if (etypes[i])
            etypes[i]->project_L2(this, i, ndm, &(Mdiag[0]), 
                                  E, mean_power_func);
    for (int i = 0; i < numnp(); ++i)
        if (Mdiag[i])
            for (int j = 0; j < ndm; ++j)
                E[i*ndm+j] /= Mdiag[i];
}


/*@T ----------
 * \subsection{Boundary condition manipulations}
 *
 * The boundary conditions are specified through the Lua fields 
 * [[bcfunc]], [[shapegbc]], and [[elementsbc]].  This is a little
 * bit of a relic, though, as all the actual work of setting up the
 * boundary conditions is now done through the [[form_bcs]]
 * function.  Since nothing in the C++ code cares directly about
 * anything but [[form_bcs]], why do we still have C++ setters for
 * these other functions?
 *
 **@c*/


void Mesh::set_bc(const char* funcname)
{
    lua_getglobal(L, funcname);
    lua_setfield("bcfunc");
}


void Mesh::set_globals_bc(const char* funcname)
{
    lua_getglobal(L, funcname);
    lua_setfield("shapegbc");
}


void Mesh::set_elements_bc(const char* funcname)
{
    lua_getglobal(L, funcname);
    lua_setfield("elementsbc");
}


void Mesh::apply_bc()
{
    apply_lua_bc();

    // Reassign IDs
    assign_ids();

    // Copy boundary data into place
    set_bc_u(vecs[VEC_U], vecs[VEC_BV] );
    set_bc_u(vecs[VEC_UI],vecs[VEC_BVI]);
    set_bc_f(vecs[VEC_F], vecs[VEC_BV] );
    set_bc_f(vecs[VEC_FI],vecs[VEC_BVI]);

    // Update node disps that are linked to globals
    dvec& U  = vecs[VEC_U];
    dvec& Ui = vecs[VEC_UI];
    for (int j = 0; j < numglobals(); ++j) {
        shapeg(&(U [0]), U [NG[0]+j], j, 0);
        shapeg(&(Ui[0]), Ui[NG[0]+j], j, 0);
    }
}


void Mesh::apply_lua_bc()
{
    CHECK_LUA;
    lua_method("form_bcs", 0, 0);
}


/*@T ----------
 * \subsection{Lua accessors}
 *
 * We use Lua functions to define the drive and sense vectors used in
 * transfer function evaluation, auxiliary fields that we might want
 * to plot, and the global shape functions.
 * 
 **@c*/


void Mesh::get_vector(const char* fname, double* v, int is_reduced)
{
    CHECK_LUA;
    QVecAssembler va(v, NULL, this, is_reduced);
    lua_getglobal(L, fname);
    tolua_pushusertype(L, this, "Mesh");
    tolua_pushusertype(L, &va, "QVecAssembler");
    lua_pcall2(L, 2, 0);
}


void Mesh::get_lua_fields(const char* funcname, int n, double* fout)
{
    lua_pushstring(L, funcname);
    lua_gettable(L, LUA_GLOBALSINDEX);
    if (!lua_isfunction(L,-1)) {
        printf("Error in get_lua_fields: [%s]\n", funcname);
        lua_pop(L,1);
        return;
    }

    int func = lua_gettop(L);
    int N    = numnp();
    fill(fout, fout+n*N, 0);

    for (int j = 0; j < N; ++j) {
        int t = lua_gettop(L)+1;
        lua_pushvalue(L, func);
        for (int i = 0; i < ndm; ++i)
            lua_pushnumber(L, x(i,j));
        if (lua_pcall2(L, ndm, n) == 0) {
            for (int i = 0; i < n && i+t < (int) lua_gettop(L)+1; ++i)
                fout[j*n+i] = lua_tonumber(L, i+t);
        }
        lua_pop(L, n);
    }
    lua_pop(L,1);
}


double Mesh::shapeg(int i, int j)
{
    CHECK_LUA 0;
    lua_pushnumber(L, i);
    lua_pushnumber(L, j);
    double retval = 0;
    if (lua_method("get_shapeg", 2, 1) == 0) {
        retval = lua_tonumber(L,-1);
        lua_pop(L,1);
    }
    return retval;
}


void Mesh::shapeg(double* v, double c, int j, int is_reduced, int vstride)
{
    if (!L || c == 0)
        return;

    QVecAssembler va(v, NULL, this, is_reduced, vstride);
    tolua_pushusertype(L, &va, "QVecAssembler");
    lua_pushnumber(L, c);
    lua_pushnumber(L, j);
    lua_method("get_shapeg_vec", 3, 0);
}


/*@T ----------
 * \subsection{Setting and getting [[U, V, A, F]]}
 *
 * There are a lot of variants of setting and getting the system state
 * vectors, simply because they're all potentially complex-valued.
 * These functions all copy to/from reduced vectors, with part of the
 * data merged in from the boundary value arrays.
 * 
 *@q*/


void Mesh::set_u(double* u, double* v, double* a)
{
    clear(vecs[VEC_U]);
    clear(vecs[VEC_V]);
    clear(vecs[VEC_A]);
    set_bc_u(vecs[VEC_U], vecs[VEC_BV]);
    set_data(u, vecs[VEC_U]);
    set_data(v, vecs[VEC_V]);
    set_data(a, vecs[VEC_A]);
}


void Mesh::set_ui(double* u)
{
    clear(vecs[VEC_UI]);
    set_bc_u(vecs[VEC_UI], vecs[VEC_BVI]);
    set_data(u, vecs[VEC_UI]);
}


void Mesh::set_f(double* f)
{
    clear(vecs[VEC_F]);
    set_bc_f(vecs[VEC_F], vecs[VEC_BV]);
    set_data(f, vecs[VEC_F]);
}


void Mesh::set_fi(double* f)
{
    clear(vecs[VEC_FI]);
    set_bc_f(vecs[VEC_FI], vecs[VEC_BVI]);
    set_data(f, vecs[VEC_FI]);
}


void Mesh::set_bc_u(vector<double>& U, vector<double>& BV)
{
    int N = BC.size();
    for (int i = 0; i < N; ++i)
        if (BC[i] == 'u')
            U[i] = BV[i];
}


void Mesh::set_bc_f(vector<double>& F, vector<double>& BV)
{
    int N = BC.size();
    for (int i = 0; i < N; ++i)
        if (BC[i] == 'f')
            F[i] = -BV[i];
}


void Mesh::set_data(double* u, vector<double>& U)
{
    if (!u)
        return;
    int M = ID.size();
    for (int i = 0; i < M; ++i) {
        if (ID[i] >= 0)
            U[i] = u[ID[i]];
    }
    for (int j = 0; j < numglobals(); ++j)
        shapeg(&(U[0]), U[NG[0]+j], j, 0);
}


void Mesh::get_reduced(double* u, int v)
{
    get_data(u, vecs[v]);
}


void Mesh::get_data(double* u, vector<double>& U)
{
    if (!u)
        return;
    int M = ID.size();
    for (int i = 0; i < M; ++i) {
        if (ID[i] >= 0)
            u[ID[i]] = U[i];
    }
    for (int j = 0; j < numglobals(); ++j)
        shapeg(u, -U[NG[0]+j], j, 1);
}


/*@T ----------
 * \subsection{Lua support methods}
 *
 * These routines support access to the Lua side of the mesh object.
 * The [[lua_getfield]], [[lua_setfield]], and [[lua_method]]
 * methods let us access fields and methods associated with the Lua
 * object that represents the mesh from tolua.
 *
 **@c*/


void Mesh::lua_getfield(const char* name)
{
    tolua_pushusertype(L, this, "Mesh");
    lua_pushstring(L, name);
    lua_gettable(L, -2);
    lua_remove(L, -2);
}


void Mesh::lua_setfield(const char* name)
{
    tolua_pushusertype(L, this, "Mesh");
    lua_insert(L,-2);
    lua_pushstring(L, name);
    lua_insert(L,-2);
    lua_settable(L,-3);
    lua_pop(L,1);
}


int Mesh::lua_method(const char* name, int nargin, int nargout)
{
    lua_getfield(name);
    if (!lua_isfunction(L,-1)) {
        lua_pop(L,nargin+1);
        return -1;
    }
    lua_insert(L,-1-nargin);
    tolua_pushusertype(L, this, "Mesh");
    lua_insert(L,-1-nargin);
    return lua_pcall2(L, nargin+1, nargout);
}
