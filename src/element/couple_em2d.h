/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef COUPLE_EM2D_H
#define COUPLE_EM2D_H

#include "pmlelement.h"
#include "qmatrix.h"


/** Two-dimensional coupled electromechanical element
 */
class CoupleEM2d : public Element {
 public:

    /** Construct a coupled electromechanical element
     */
    CoupleEM2d(double kappa);
    ~CoupleEM2d();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double kappa;
    int nen;
};


#endif /* COUPLE_EM2D_H */
