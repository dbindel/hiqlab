/* HiQLab
 * Copyright (c): Regents of the University of California
 *
 */

#ifndef METISPART_H
#define METISPART_H

#include <vector>
 
extern "C" {
#include "macros.h"
#include "defs.h"
#include "struct.h"
#include "proto.h"
}

class Metis {

  public:
    Metis(int n_, idxtype* xadj_, idxtype* adjncy_);
    virtual ~Metis();

    void set_vwgt(  idxtype* vwgt_)   {vwgt=vwgt_;}
    void set_adjwgt(idxtype* adjwgt_) {adjwgt=adjwgt_;}
    void set_numflag(int numflag_)    {numflag=numflag_;}
    void set_wgtflag(int wgtflag_)    {wgtflag=wgtflag_;}
    void set_options(int* options_)   {options=options_;}

    int  get_n()       {return n;}
    int  get_ne()      {return xadj[n];}
    int  get_numflag() {return numflag;}
    int  get_wgtflag() {return wgtflag;}
    int* get_options() {return options;}

    void PartGraphRecursive(int nparts, idxtype* part, int* edgecut, int* options=NULL);
    void PartGraphKway(     int nparts, idxtype* part, int* edgecut, int* options=NULL);
    void EdgeND(idxtype* perm, idxtype* iperm, int* options=NULL);
    void NodeND(idxtype* perm, idxtype* iperm, int* options=NULL);

  private:

    int      n;
    int*     options;
    int      wgtflag;
    int      numflag;

    idxtype* xadj;
    idxtype* adjncy;
    idxtype* vwgt;
    idxtype* adjwgt;

    void check_part_arguments();
};

#endif /* METISPART_H */
