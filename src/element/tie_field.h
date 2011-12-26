/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef TIE_FIELD_H
#define TIE_FIELD_H

extern "C" {
#include <lua.h>
}

#include "element.h"
#include "qmatrix.h"


/** General "tie" using Lagrange multipliers
 */
class TieFieldElement : public Element {
 public:

    /** Construct a new tie element
     */
    TieFieldElement(lua_State* L, int func);
    ~TieFieldElement();

    void initialize(Mesh* mesh, int eltid);
    void assign_ids(Mesh* mesh, int eltid);
    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    lua_State* L;
    int nbranch;
    void call_func(Mesh* mesh, int eltid, int mode,
                   double cx, QAssembler* K);
};


#endif /* TIE_FIELD_H */
