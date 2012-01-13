#ifndef MESH_PARTITION_UTIL_H
#define MESH_PARTITION_UTIL_H

#include <cstdio>
#include <string>

#include "mesh_partitioner.h"
#include "mesh_partition.h"

class Mesh_Partition_Writer {

  public:
    Mesh_Partition_Writer() {}
    virtual ~Mesh_Partition_Writer();

    void write_mesh(Mesh_Partitioner* mp, const char* fname);
    void write_partition(Mesh_Partitioner* mp, int ni, const char* fname);
    void write_nemap(Mesh_Partitioner* mp, int* nemap, const char* fname);
    void write_partition_nemap(Mesh_Partitioner* mp, int ni, int* nemap, const char* fname);

  private:
    Mesh_Partitioner* meshp;
    FILE* fpa;
    FILE* fpb;
    int   bytecount;

    void open_file(int ni, const char* fname);
    void close_file();

    void write_basic_info(std::string fnameb);
    void write_node_info(int ni);
    void write_element_info(int ni);
    void write_array(const char* aname, int m, int n, int* data);

    void write_basicid_info(int ni);
    void write_nodeid_info(int ni);

    void write_nemap(int ni, int* nemap);
};

class Mesh_Partition_Reader {

  public:
    Mesh_Partition_Reader() {}
    virtual ~Mesh_Partition_Reader();

    void read_partition(Mesh_Partition* mesh, const char* fname);
    void read_partition_ids(Mesh_Partition* mesh, const char* fname);
    int* read_partition_nemap(Mesh_Partition* mesh, const char* fname);

  private:
    Mesh_Partition* meshnp;
    FILE* fpa;
    FILE* fpb;
    int   NDM, NDF, NEN;
    int   startid, end1id;

    void open_file(int ni, const char* fname);
    void close_file();

    void read_basic_info(std::string fnameb);
    void read_node_info();
    void read_element_info();
    int* read_array(const char* aname, int* m, int* n);

    void read_basicid_info();
    void read_nodeid_info();
};

#endif /* MESH_PARTITION_UTIL_H */
