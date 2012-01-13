/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 * $Id: metisstubs.h,v 1.12 2006/06/19 16:56:11 dbindel Exp $
 */

#ifndef METISSTUBS_H
#define METISSTUBS_H
 
#include <mex.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "macros.h"
#include "defs.h"
#include "struct.h"

void METIS_PartGraphRecursive(int* n, idxtype* xadj, idxtype* adjncy, 
                                 idxtype* vwgt, idxtype* adjwgt, 
                         int* wgtflag, int* numflag, int* nparts, 
                         int* options, int* edgecut, idxtype* part); 

void METIS_PartGraphKway(int* n, idxtype* xadj, idxtype* adjncy, 
                                 idxtype* vwgt, idxtype* adjwgt, 
                         int* wgtflag, int* numflag, int* nparts, 
                         int* options, int* edgecut, idxtype* part); 
void METIS_EdgeND(int* n, idxtype* xadj, idxtype* adjncy, int* numflag, 
                  int* options, idxtype* perm, idxtype* iperm);
void METIS_NodeND(int* n, idxtype* xadj, idxtype* adjncy, int* numflag, 
                  int* options, idxtype* perm, idxtype* iperm);

#ifdef __cplusplus
}
#endif

// -- metis_PartGraph(Case for both weights)
void metis_PartGraph(int prout, int* xadj, int* adjncy,
                         int* part, int* edgecut,
                         int nparts, int wgtflag, int* options,
                         int nv, int ne, int ncon,
                         int* vwgt, int nvwgt, int* adjwgt, int nadjwgt);

// -- metis_PartGraph(Case for no weights)
void metis_PartGraph(int prout, int* xadj, int* adjncy,
                         int* part, int* edgecut,
                         int nparts, int wgtflag, int* options,
                         int nv, int ne);

// -- metis_PartGraph(Case for edge weights only)
void metis_PartGraph(int prout, int* xadj, int* adjncy,
                         int* part, int* edgecut,
                         int nparts, int wgtflag, int* options,
                         int nv, int ne,
                         int* adjwgt, int nadjwgt);

// -- metis_PartGraph(Case for node weights only)
void metis_PartGraph(int prout, int* xadj, int* adjncy,
                         int* part, int* edgecut,
                         int nparts, int wgtflag, int* options,
                         int nv, int ne, int ncon,
                         int* vwgt, int nvwgt);

void metis_ND(int prout, int* xadj, int* adjncy,
                         int* perm, int* iperm,
                         int* options,
                         int nv, int ne);

#endif /* METISSTUBS_H */
