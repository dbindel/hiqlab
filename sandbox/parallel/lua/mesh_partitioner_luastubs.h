#ifndef MESH_PARTITIONER_LUASTUBS_H
#define MESH_PARTITIONER_LUASTUBS_H

class Mesh_Partitioner;
class Mesh_Add_Block;

class PEPair_Data {

  public:
    PEPair_Data(Mesh_Partitioner* mp, Mesh_Add_Block* from, Mesh_Add_Block* to);
    ~PEPair_Data();

    int  get_numnp_t() {return numnp_t;}
    int* get_nemap()   {return nemap;}
    int* get_npart_t() {return npart_t;}

  private:
    int  numnp_t;
    int* nemap;
    int* npart_t;
};

#endif /* MESH_PARTITIONER_LUASTUBS_H */
