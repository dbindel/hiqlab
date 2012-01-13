#include "mesh_partitioner_luastubs.h"
#include "mesh_add_block.h"
#include "mesh_partitioner_metis.h"
#include "mesh_partitioner_conform.h"

#define ME PEPair_Data

ME::ME(Mesh_Partitioner* mp, Mesh_Add_Block* from, Mesh_Add_Block* to) 
       : numnp_t(0), nemap(0), npart_t(0)
{
    numnp_t = to->numnp();
    nemap   = new int[numnp_t];
    npart_t = new int[numnp_t];

    to->node_element_mapping(from,nemap);

    Mesh_Partitioner_METIS*   mpm;
    Mesh_Partitioner_Conform* mpc;
    if (mpm = dynamic_cast<Mesh_Partitioner_METIS*>(mp))
        mpm->get_node_partition(npart_t);
    else if (mpc = dynamic_cast<Mesh_Partitioner_Conform*>(mp))
        mpc->get_node_partition(npart_t);
}


ME::~ME()
{
    if (nemap)
        delete[] nemap;
    if (npart_t)
        delete[] npart_t;
}
