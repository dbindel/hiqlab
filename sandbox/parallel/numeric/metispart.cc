#include "metispart.h"
#include <cassert>

#define ME Metis

ME::ME(int n_, idxtype* xadj_, idxtype* adjncy_) 
  : n(n_), xadj(xadj_), adjncy(adjncy_), numflag(0), wgtflag(0),
    vwgt(NULL), adjwgt(NULL)
{
}

ME::~ME()
{
}

void ME::PartGraphRecursive(int nparts, idxtype* part, int* edgecut, int* options_)
{

    check_part_arguments();
    if (options_)
        options=options_;
    else {
        options=new int[5];
        options[0] = 0;
        options[1] = 3;
        options[2] = 1;
        options[3] = 1;
        options[4] = 0;
    }

    METIS_PartGraphRecursive(&n, xadj, adjncy, vwgt, adjwgt, 
                         &wgtflag, &numflag, &nparts, options, edgecut, part);

    if (options_==NULL) {
        delete[] options;
    }
}

void ME::PartGraphKway(int nparts, idxtype* part, int* edgecut, int* options_)
{

    check_part_arguments();
    if (options_)
        options=options_;
    else {
        options=new int[5];
        options[0] = 0;
        options[1] = 3;
        options[2] = 1;
        options[3] = 3;
        options[4] = 0;
    }

    METIS_PartGraphKway(&n, xadj, adjncy, vwgt, adjwgt, 
                         &wgtflag, &numflag, &nparts, 
                         options, edgecut, part); 

    if (options_==NULL) {
        delete[] options;
    }
}

void ME::EdgeND(idxtype* perm, idxtype* iperm, int* options_)
{
    if (options_)
        options=options_;
    else {
        options=new int[5];
        options[0] = 0;
        options[1] = 3;
        options[2] = 1;
        options[3] = 1;
        options[4] = 0;
    }

    METIS_EdgeND(&n, xadj, adjncy, &numflag, options, perm, iperm);

    if (options_==NULL) {
        delete[] options;
    }
}

void ME::NodeND(idxtype* perm, idxtype* iperm, int* options_)
{
    if (options_)
        options=options_;
    else {
        options=new int[5];
        options[0] = 0;
        options[1] = 3;
        options[2] = 1;
        options[3] = 1;
        options[4] = 0;
    }

    METIS_NodeND(&n, xadj, adjncy, &numflag, options, perm, iperm);

    if (options_==NULL) {
        delete[] options;
    }
}

void ME::check_part_arguments()
{
    if (wgtflag==1) {

        // -- Edge weights only 
        assert(adjwgt!=NULL);

    } else if (wgtflag==2) {

        // -- Node weights only 
        assert(vwgt!=NULL);

    } else if (wgtflag==3) {

        // -- Both weights
        assert(adjwgt!=NULL);
        assert(vwgt!=NULL);
    }
}
