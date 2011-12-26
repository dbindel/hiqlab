/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include "mesh.h"
#include "meshlabeler.h"

#define ME MeshLabeler


int ME::label_nodal(int first_slot, int past_slot)
{
    int np = mesh->numnp();
    int numid_local = 0;
    for (int j = 0; j < np; ++j)
        for (int i = first_slot; i < past_slot; ++i)
            if (mesh->id(i,j) >= 0) {
                mesh->id(i,j) = numid++;
                ++numid_local;
            }
    return numid_local;
}


void ME::label_global()
{
    int ngb = mesh->numglobals();
    for (int j = 0; j < ngb; ++j)
        if (mesh->globalid(j) >= 0)
            mesh->globalid(j) = numid++;
}


void ME::label_branch()
{
    int numelt  = mesh->numelt();
    for (int j = 0; j < numelt; ++j) {
        int nb_id = mesh->nbranch_id(j);
        for (int i = 0; i < nb_id; ++i)
            if (mesh->branchid(i,j) >= 0)
                mesh->branchid(i,j) = numid++;
    }
}

