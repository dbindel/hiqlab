#ifndef _TRILINOS_INDEXMAP_H
#define _TRILINOS_INDEXMAP_H

class IndexMap {

 public:

    virtual int row(int id) = 0;
    virtual int col(int id) = 0;

};



#endif

