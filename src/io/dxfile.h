/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef DXFILE_H
#define DXFILE_H

#include <cstdio>
#include <string>

#include "mesh.h"

class DXFile {
public:
    DXFile(const char* name);
    ~DXFile();

    void start_field(const char* name = NULL);
    void field_component(const char* name, int id);
    void end_field();
    void start_array(int m, int n, float* A);
    void start_array(int m, int n, double* A);
    void start_array(int m, int n, int* A);
    void start_array(int m, int n, float* A, int numid, int* id);
    void start_array(int m, int n, double* A, int numid, int* id);
    void write(float x);
    void write(double x);
    void write(int x);
    void array_attribute(const char* key, const char* value);
    void end_array();
    void writemesh(Mesh* mesh);

private:
    std::string basename;
    FILE* fid1;
    FILE* fid2;
    int count;
    int objcount;

    int writer_state;  // 0 for neutral, 1 in field, 2 in array
    char item_type;    // i for int, f for float
    int  item_left;    // Number of items remaining in array

    void change_state(int state);
    void writemesh2d(Mesh* mesh);
    void writemesh3d(Mesh* mesh);
};

#endif /* DXFILE_H */
