/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 * $Id: coordmatrix_petsc.h,v 1.3 2006/06/18 04:07:11 tkoyama Exp $
 */
 
#ifndef COORDMATRIX_PETSC_H
#define COORDMATRIX_PETSC_H

#include "qcomplex.h"
#include "cscmatrix.h"
#include "coordmatrix.h"

#include <vector>

/** Matrix assembler for complex mass and stiffnesses.
 */
class CoordMatrix_Petsc : public CoordMatrix {
 public:

    /** Create a new assembler for mass and stiffness matrices
     *
     * @param N  dimension of the system
     */
    CoordMatrix_Petsc(int N);
    CoordMatrix_Petsc(int M, int N);
    ~CoordMatrix();

    /** Set functions to map row and column indices
     *
     * @param f_row  Row mapping function
     * @param f_col  Column mapping function
     *
     */

    void set_mapping_functions(int (* f_row)(int), int (* f_col)(int) );

    /** Add local mass and stiffness contributions
     *
     * @param eltid  Global identifiers for element dofs (n)
     * @param n      Number of element dofs
     * @param Ke     Element matrix (n-by-n)
     */
    void add(int* eltid, int n, dcomplex* Ke);
    void add(int* eltid, int n, double*   Ke);
    
 private:
    int (* f_row)(int);
    int (* f_col)(int);
};


#endif /* COORDMATRIX_PETSC_H */
