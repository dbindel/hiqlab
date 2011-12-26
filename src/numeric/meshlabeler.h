/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef MESHLABELER_H
#define MESHLABELER_H

#include "mesh.h"

class MeshLabeler {
public:
    MeshLabeler(Mesh* mesh) : mesh(mesh), numid(0) {}
    int label_nodal(int first_slot, int past_slot);
    void label_global();
    void label_branch();
private:
    Mesh* mesh;
    int numid;
};

#endif /* MESHLABELER_H */
