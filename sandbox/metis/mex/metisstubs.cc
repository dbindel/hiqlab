#include "metisstubs.h"
#include <string.h>
#include <iostream>

void metis_PartGraph(int prout, int* xadj, int* adjncy,
                         int* part, int* edgecut,
                         int nparts, int wgtflag, int* options,
                         int nv, int ne, int ncon,
                         int* vwgt, int nvwgt, int* adjwgt, int nadjwgt)
{
    int numflag = 0;

    memset(part, 0, nv*sizeof(int));

    if (prout==0)
        METIS_PartGraphRecursive(&nv, xadj, adjncy, vwgt, adjwgt, 
                            &wgtflag, &numflag, &nparts, options, edgecut, part); 
    else if (prout==1)
        METIS_PartGraphKway(&nv, xadj, adjncy, vwgt, adjwgt, 
                            &wgtflag, &numflag, &nparts, options, edgecut, part); 
}

void metis_PartGraph(int prout, int* xadj, int* adjncy,
                         int* part, int* edgecut,
                         int nparts, int wgtflag, int* options,
                         int nv, int ne)
{
    metis_PartGraph(prout, xadj, adjncy,
                         part, edgecut,
                         nparts, wgtflag, options,
                         nv, ne, 0,
                         NULL, 0, NULL, 0);
}

void metis_PartGraph(int prout, int* xadj, int* adjncy,
                         int* part, int* edgecut,
                         int nparts, int wgtflag, int* options,
                         int nv, int ne,
                         int* adjwgt, int nadjwgt)
{
    metis_PartGraph(prout, xadj, adjncy,
                         part, edgecut,
                         nparts, wgtflag, options,
                         nv, ne, 0,
                         NULL, 0, adjwgt, nadjwgt);
}

void metis_PartGraph(int prout, int* xadj, int* adjncy,
                         int* part, int* edgecut,
                         int nparts, int wgtflag, int* options,
                         int nv, int ne, int ncon,
                         int* vwgt, int nvwgt)
{
    metis_PartGraph(prout, xadj, adjncy,
                         part, edgecut,
                         nparts, wgtflag, options,
                         nv, ne, ncon,
                         vwgt, nvwgt, NULL, 0);
}

void metis_ND(int prout, int* xadj, int* adjncy,
                         int* perm, int* iperm,
                         int* options,
                         int nv, int ne)
{
    int numflag = 0;

    memset(perm,  0, nv*sizeof(int));
    memset(iperm, 0, nv*sizeof(int));

    if (prout==0)
        METIS_NodeND(&nv,xadj,adjncy,&numflag,options,perm,iperm);
    else if (prout==1)
        METIS_EdgeND(&nv,xadj,adjncy,&numflag,options,perm,iperm);
}
