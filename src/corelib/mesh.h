/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef MESH_H
#define MESH_H

extern "C" {
  #include <lua.h>
}

#include <vector>

#include "qcomplex.h"
#include "element.h"
#include "qassembly.h"

#define defq_mesh_accessor(type, name, array) \
type& name        (int i)        { return array[i];            } \
type& name        (int i, int j) { return array[inode(i,j)];   } \
type& branch##name(int i, int j) { return array[ibranch(i,j)]; } \
type& global##name(int i)        { return array[iglobal(i)];   }

#define defq_meshd_accessor(name, array) \
double& name        (int i)        { return vec(VEC_##array, i);        } \
double& name        (int i, int j) { return vec(VEC_##array,i,j);       } \
double& branch##name(int i, int j) { return branchvec(VEC_##array,i,j); } \
double& global##name(int i)        { return globalvec(VEC_##array,i);   }

#define defq_meshz_accessor(name, array) \
defq_meshd_accessor(name, array) \
defq_meshd_accessor(name##i, array##I)


/*@T -----------
 * \subsection{Mesh declarations}
 *@c*/

class Mesh {
public:

    enum VecList { 
        VEC_U = 0,  // Displacement
        VEC_V,      // Velocity
        VEC_A,      // Acceleration
        VEC_F,      // Force
        VEC_BV,     // Boundary values
        VEC_UI,     // Imaginary part of displacement
        VEC_FI,     // Imaginary part of force
        VEC_BVI,    // Imaginary part of boundary values
        VEC_D1,     // Primal variable scaling
        VEC_D2,     // Dual variable scaling
        VEC_LEN
    };

    Mesh(int ndm);
    virtual ~Mesh();

    /** Get Lua state used by (not necessarily owned by) system */
    lua_State* get_lua()                { return L; }
    void       set_lua(lua_State* argL) { L = argL; }

    /** Get out a named scale from the Lua state */
    double get_scale(const char* name);

    /** Set scaling vector D1 and D2 */
    void clear_scaling_vector();
    void set_nodal_u_scale(int i, double s);
    void set_nodal_f_scale(int i, double s);
    void set_scaling_vector();

    /** Get scaling vector D1 or D2 */
    void get_scaling_vector(double* su, double* sf, int is_reduced = 1);

    /** Add node(s), element(s), or global(s) to the mesh directly.
     *  In each case, return the identifier of the first item created.
     */
    int add_node(double* x, int n=1);
    int add_element(int* e, Element* etype, int nen, int n=1);
    int add_global(int n);

    /** Construct a Cartesian block. */
    void add_block(double* x1, double* x2, int* m,
                   Element* etype, int order=1);
    void add_block(double x1, double x2, int m,
                   Element* etype, int order=1);
    void add_block(double x1, double y1,
                   double x2, double y2,
                   int nx, int ny,
                   Element* etype, int order=1);
    void add_block(double x1, double y1, double z1,
                   double x2, double y2, double z2,
                   int nx, int ny, int nz,
                   Element* etype, int order=1);

    /** Tie mesh nodes in a range. */
    void tie(double tol, int start=-1, int end=-1);

    /** Allocate space for global arrays and assign initial index map. */
    void initialize();

    /** Reassign the index map. */
    int assign_ids();

    /** Assemble tangent matrix structure. */
    void assemble_struct(QStructAssembler* K, int reduced = 1);
    void assemble_struct_raw(QStructAssembler* K);

    /** Assemble global mass and stiffness matrices. */
    void assemble_dR(QAssembler* K, double cx=1, double cv=0, double ca=0,
                     int reduced = 1);
    void assemble_dR_raw(QAssembler* K, double cx=1, double cv=0, double ca=0);

    /** Assemble global reaction vector (goes into f) */
    void assemble_R();

    /** Assemble time-averaged energy flux at node points */
    void mean_power(double* E);

    /** Use Lua to set boundary conditions */
    void set_bc         (const char* funcname);
    void set_globals_bc (const char* funcname);
    void set_elements_bc(const char* funcname);

    /** Apply boundary conditions based on saved Lua functions */
    void apply_bc();

    /** Build alternate vectors (for drive or sense) using Lua */
    void get_vector(const char* fname, double* v, int is_reduced = 1);

    /** Evaluate a Lua function at every nodal point */
    void get_lua_fields(const char* funcname, int n, double* fout);

    /** Evaluate global shape functions via Lua funcs */
    double shapeg(int i, int j);
    void   shapeg(double* v, double c, int j,
                  int is_reduced, int vstride = 1);

    /** Copy i*omega*u and -omega^2*u into a and v */
    dcomplex harmonic_freq()  { return harmonic_freq_; }
    bool     is_harmonic()    { return harmonic_flag_; }
    void     clear_harmonic() { harmonic_flag_ = 0;    }
    
    void make_harmonic(dcomplex omega) {
        harmonic_flag_ = 1;
        harmonic_freq_ = omega;
    }

    void make_harmonic(double omega_r, double omega_i = 0) {
        make_harmonic(dcomplex(omega_r, omega_i));
    }

    int numnp()         { return X.size()  / ndm;    }
    int numelt()        { return etypes.size();      }
    int numglobals()    { return NG[1]-NG[0];        }
    int get_numid()     { return numid;              }
    int get_idsize()    { return ID.size();          }
    int get_ndm()       { return ndm;                }
    int get_ndf()       { return maxndf;             }
    int get_nen()       { return maxnen;             }
    int nbranch_id()    { return NB[numelt()]-NB[0]; }

    int       get_nen    (int i) { return NIX[i+1]-NIX[i]; }
    int       get_nne    (int i) { return NIN[i+1]-NIN[i]; }
    int       nbranch_id (int i) { return NB[i+1] -NB[i];  }
    int       nhist_id   (int i) { return NH[i+1] -NH[i];  }
    int&      nbranch    (int i) { return NB[i];           }
    int&      nhist      (int i) { return NH[i];           }
    Element*& etype      (int i) { return etypes[i];       }

    int&      ix(int i, int j) { return IX[NIX[j] + i]; }
    int&      in(int i, int j) { return IN[NIN[j] + i]; }
    double&   x (int i, int j) { return X[j*ndm + i];   }

    /** Get the global array indices for different variables */
    int inode  (int i, int j) { return j*maxndf+i; }
    int ibranch(int i, int j) { return NB[j]+i;    }
    int iglobal(int i)        { return NG[0]+i;    }

    /** Get and set the ID and BC arrays */
    defq_mesh_accessor( int,  id, ID )
    defq_mesh_accessor( char, bc, BC )

    /** Get and set the boundary values, u, v, a, f, and scaling */
    double& vec(int v, int i)              { return vecs[v][i];            }
    double& vec(int v, int i, int j)       { return vecs[v][inode(i,j)];   }
    double& branchvec(int v, int i, int j) { return vecs[v][ibranch(i,j)]; }
    double& globalvec(int v, int i)        { return vecs[v][iglobal(i)];   }

    defq_meshz_accessor( bv,  BV  )
    defq_meshz_accessor( u ,  U   )
    defq_meshd_accessor( v ,  V   )
    defq_meshd_accessor( a ,  A   )
    defq_meshz_accessor( f ,  F   )
    defq_meshd_accessor( d1,  D1  )
    defq_meshd_accessor( d2,  D2  )

    void set_u (double* u = NULL, double* v = NULL, double* a = NULL);
    void set_ui(double* u = NULL);

    void set_f (double* f = NULL);
    void set_fi(double* f = NULL);

    void get_reduced(double* u, int v);

    /** Access element history entries */
    double get_hist(int i, int j)           { return hist(i,j,0); }
    double put_hist(int i, int j, double x) { hist(i,j,1) = x;    }
    void   swap_hist()                      { which_hist = (which_hist+1)%2; }

    Element* own(Element* e);
    void     own(lua_State* L);

    void lua_getfield(const char* name);
    void lua_setfield(const char* name);
    int  lua_method(const char* name, int nargin, int nargout);

    /*@q ----------
     * Private data structures are documented in mesh.cc
     */
private:

    typedef std::vector<double>   dvec;
    typedef std::vector<int>      ivec;
    typedef std::vector<char>     cvec;
    typedef std::vector<Element*> evec;

    int ndm;           // Number of spatial dimensions
    int maxnen;        // Maximum number of nodes per element
    int maxndf;        // Maximum number of dofs per node
    int numid;         // Number of active dofs (not under BCs)
    int which_hist;    // Which set of history variables is active

    dvec X;            // Node positions
    ivec IX;           // Element connectivity
    ivec IN;           // Node-element connectivity
    ivec ID;           // Variable identifiers (nodal + branch)

    ivec NIX;          // Offsets of node lists in IX
    ivec NIN;          // Offsets of element lists in IN
    ivec NB;           // Offsets of branch variables in ID, etc
    int  NG[2];        // Offsets of global variables in ID, etc
    ivec NH;           // Number of history vars

    cvec BC;           // Boundary codes

    dvec HIST;         // History (both)

    dvec vecs[VEC_LEN];

    bool     harmonic_flag_;  // Is this a time-harmonic calculation?
    dcomplex harmonic_freq_;  // What frequency?

    evec etypes;       // Types associated with individual elements
    evec etypes_owned; // Types owned by the system

    lua_State* lua_owned;
    lua_State* L;          // Contains set_bc function

    double& hist(int i, int j, int k) {
        return HIST[NH[j]+i + numelt()*((which_hist+k)%2)];
    }

    void initialize_sizes();
    void build_reverse_map();
    int  initialize_indexing();
    int  initialize_history();
    void initialize_arrays(int M, int MH);

    void set_data(double* u, dvec& U);
    void get_data(double* u, dvec& U);

    void set_bc_u(dvec& U, dvec& BV);
    void set_bc_f(dvec& U, dvec& BV);

    void end_add_element(Element* etype);

    void apply_lua_bc();
    /*@c*/
};

#undef defq_mesh_accessor
#undef defq_meshz_accessor

#endif /* MESH_H */
