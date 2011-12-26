/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef PMLELEMENT_H
#define PMLELEMENT_H

extern "C" {
#include <lua.h>
}

#include "pmlfunc.h"
#include "element.h"


/** Base class for PML elements.
 */
class PMLElement : public Element {
 public:

    PMLElement(int nslots) :
        Element(nslots), func_(NULL), stretch(NULL),
        stretch_size(0) {}

    ~PMLElement();

    /** Let a Lua function compute the stretching function.  The Lua
     *  function should have the prototype
     *    - x,y   = stretch(x,y) or
     *    - x,y,z = stretch(x,y,z)
     *  The resulting PML will transform the x and y coordinates
     *  independently as xtilde = 1 - ix and ytilde = 1-iy.  Like any
     *  Lua function, the stretch function can access anything
     *  available in the environment in which it was defined.
     *
     * @param L  A Lua interpreter object
     * @param s  The location of the stretching function on the Lua stack
     */
    void set_stretch(lua_State* L, int s);
    void set_stretch(PMLFunc* func);

    /** Define a stretch function by defining it at the node points.
     *
     * @param stretch  Array of stretch values at nodes (ndm-by-nen).
     */
    void set_stretch(double* stretch, int ndm=0, int nen=0);

 protected:

    /** Determine whether a stretch function is defined
     */
    int has_stretch() { return (stretch != 0) || (func_ != 0); }

    /** Call the Lua function to evaluate the stretch vector.  Returns
     *  a zero stretch (identity transformation) if no Lua function is
     *  defined.
     *
     * @param xx   Position on input; stretch values on output.
     * @param ndm  Number of spatial dimensions (2 or 3)
     */
    void compute_stretch(double* xx, int ndm);

    /** Evaluate the stretch vector at a nodal point, either from a
     *  table lookup or using the Lua function.
     *
     * @param mesh    Mesh used for looking up node positions for Lua call
     * @param nodeid  Node identifier
     * @param xx      Output for stretch values
     * @param ndm     Spatial dimension
     */
    void compute_stretch(Mesh* mesh, int nodeid, double* xx, int ndm);

    /** Get the element stretch array.
     *
     * @param nodestretch  Output array of stretch values (ndm-by-nen)
     * @param nen          Number of element nodes
     * @param ndm          Spatial dimension
     */
    void set_local_stretch(Mesh* mesh, int eltid, dcomplex* nodestretch,
                           int nen, int ndm);

 private:
    PMLFunc* func_;
    double* stretch;
    int stretch_size;
};


#endif /* PMLELEMENT_H */
