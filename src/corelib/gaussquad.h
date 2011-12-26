/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef GAUSSQUAD_H
#define GAUSSQUAD_H

/* Compute the ith point / weight in an npts-point Gauss quadrature 
 * rule on [-1,1]
 */
double gauss_point (int i, int npts);
double gauss_weight(int i, int npts);

#endif /* GAUSSQUAD_H */
