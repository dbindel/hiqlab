/* HiQLab
 * Copyright (c): Regents of the University of California
 */


#include "../../config.h"
#if HAVE_INTTYPES_H
  #include <inttypes.h>
#else
  typedef HIQINT32 int32_t;
#endif

#include "dxfile.h"
#include "shapes.h"

#include <cassert>
#include <iostream>
#include <vector>

using namespace std;


inline string stringify(int x)
{
    ostringstream o;
    o << x;
    return o.str();
}


static const char* endian_string()
{
    float one = 1.0;
    const char* s = (const char*) &one;
    if (s[0] == 0)
        return "lsb";
    else
        return "msb";
}


DXFile::DXFile(const char* name) :
    basename(name)
{
    string dxname  = basename;
    string binname = basename;
    dxname  += ".dx";
    binname += ".bin";
    fid1 = fopen(dxname.c_str(),  "w" );
    fid2 = fopen(binname.c_str(), "wb");
    count    = 0;
    objcount = 0;
    writer_state = 0;
}


DXFile::~DXFile()
{
    change_state(0);
    fprintf(fid1, "end\n");
    fclose(fid2);
    fclose(fid1);
}


void DXFile::start_field(const char* name)
{
    change_state(1);
    ++objcount;
    if (name)
        fprintf(fid1, "object \"%s\" class field\n", name);
    else
        fprintf(fid1, "object %d class field\n", objcount);
}


void DXFile::field_component(const char* name, int id)
{
    assert(writer_state == 1);
    fprintf(fid1, "component \"%s\" value %d\n", name, id);
}


void DXFile::end_field()
{
    assert(writer_state == 1);
    fprintf(fid1, "#\n");
}


void DXFile::start_array(int m, int n, float* A)
{
    change_state(2);
    item_type = 'f';
    item_left = m*n;

    ++objcount;
    fprintf(fid1, "object %d class array type float\n", objcount);
    if (m == 1 || n == 1)
        fprintf(fid1, "rank 0 items %d\n", m*n);
    else
        fprintf(fid1, "rank 1 shape %d items %d\n", m, n);
    fprintf(fid1, "%s binary data file %s.bin,%d\n",
            endian_string(), basename.c_str(), count*4);
    if (A) {
        fwrite(A, m*n, sizeof(float), fid2);
        item_left = 0;
    }
    count += m*n;
}


void DXFile::start_array(int m, int n, double* A)
{
    if (A) {
        vector<float> A2(m*n);
        for (int i = 0; i < m*n; ++i){
            A2[i] = (float) A[i];
        }
        start_array(m, n, &(A2[0]));
    } else {
        start_array(m, n, (float*) NULL);
    }
}


void DXFile::start_array(int m, int n, int* A)
{
    change_state(2);
    item_type = 'i';
    item_left = m*n;

    ++objcount;
    fprintf(fid1, "object %d class array type int\n", objcount);
    if (m == 1 || n == 1)
        fprintf(fid1, "rank 0 items %d\n", m*n);
    else
        fprintf(fid1, "rank 1 shape %d items %d\n", m, n);
    fprintf(fid1, "%s binary data file %s.bin,%d\n",
            endian_string(), basename.c_str(), count*4);
    if (A) {
        for (int i = 0; i < m*n; ++i) {
            int32_t Atmp = (int32_t) A[i];
            fwrite(&Atmp, 1, sizeof(int32_t), fid2);
        }
        item_left = 0;
    }
    count += m*n;
}


void DXFile::start_array(int m, int n, float* A, int numid, int* id)
{
    change_state(2);
    item_type = 'f';
    item_left = numid*n;

    ++objcount;
    fprintf(fid1, "object %d class array type float\n", objcount);
    if (numid == 1 || n == 1)
        fprintf(fid1, "rank 0 items %d\n", numid*n);
    else
        fprintf(fid1, "rank 1 shape %d items %d\n", numid, n);
    fprintf(fid1, "%s binary data file %s.bin,%d\n",
            endian_string(), basename.c_str(), count*4);
    if (A) {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < numid; ++j)
                fwrite(&(A[m*i+id[j]]), 1, sizeof(float), fid2);
        item_left = 0;
    }
    count += n*numid;
}

void DXFile::start_array(int m, int n, double* A, int numid, int* id)
{
    if (A) {
        vector<float> A2(m*n);
        for (int i = 0; i < m*n; ++i){
            A2[i] = (float) A[i];
        }
        start_array(m, n, &(A2[0]), numid, id);
    } else {
        start_array(m, n, (float*) NULL, numid, id);
    }
}


void DXFile::write(double x)
{
    write((float) x);
}


void DXFile::write(float x)
{
    assert(writer_state == 2);
    assert(item_left > 0);
    assert(item_type == 'i' || item_type == 'f');

    if (item_type == 'f') {
        fwrite(&x, 1, sizeof(float), fid2);
        --item_left;
    } else {
        write((int) x);
    }
}


void DXFile::write(int x)
{
    assert(writer_state == 2);
    assert(item_left > 0);
    assert(item_type == 'i' || item_type == 'f');

    if (item_type == 'i') {
        int32_t xtmp = (int32_t) x;
        fwrite(&xtmp, 1, sizeof(int32_t), fid2);
        --item_left;
    } else {
        write((float) x);
    }
}


void DXFile::array_attribute(const char* key, const char* value)
{
    assert(writer_state == 2);
    fprintf(fid1, "attribute \"%s\" string \"%s\"\n", key, value);
}


void DXFile::end_array()
{
    assert(writer_state == 2);
    assert(item_left == 0);
    fprintf(fid1, "#\n");
}


void DXFile::change_state(int state)
{
    if (writer_state == 1)
        end_field();
    if (writer_state == 2)
        end_array();

    writer_state = state;
}


void DXFile::writemesh(Mesh* mesh)
{
    assert(writer_state == 0);
    assert(objcount == 0);

    if (mesh->get_ndm() == 2)
        writemesh2d(mesh);
    else
        writemesh3d(mesh);

    start_field("mesh");
    field_component("positions", 1);
    field_component("connections", 2);
    field_component("displacements", 3);
    field_component("displacementsi", 4);
    for (int i = 0; i < mesh->get_ndf(); ++i) {
        string dataname = "field";
        dataname += stringify(i);
        field_component(dataname.c_str(), 5+i);
    }
    for (int i = 0; i < mesh->get_ndf(); ++i) {
        string dataname = "fieldi";
        dataname += stringify(i);
        field_component(dataname.c_str(), 5+mesh->get_ndf()+i);
    }
    end_field();
}


void DXFile::writemesh2d(Mesh* mesh)
{
    start_array( mesh->get_ndm(), mesh->numnp(), &(mesh->x(0,0)) );
    array_attribute( "dep", "positions" );
    end_array();

    int order = order2d(mesh->get_nen());
    int np = order+1;

    start_array( 4, mesh->numelt() * order*order, (int*) NULL );
    array_attribute( "element type", "quads"     );
    array_attribute( "ref",          "positions" );

    // For each element
    for (int i = 0; i < mesh->numelt(); ++i)

        // For each sub-panel of the element
        for (int ii = 0; ii < order; ++ii)
            for (int jj = 0; jj < order; ++jj)

                // Write the sub-panel connectivity
                for (int ll = 0; ll < 2; ++ll)
                    for (int mm = 0; mm < 2; ++mm) {
                        int j = (jj+ll)*np + (ii+mm);
                        write(mesh->ix(j, i));
                    }

    end_array();

    // Displacement vector
    int ind[2] = {0, 1};
    start_array( mesh->get_ndf(), mesh->numnp(), &(mesh->u(0,0)), 2, ind);
    array_attribute( "dep", "positions" );
    end_array();

    // Displacement vector(Imag)
    start_array( mesh->get_ndf(), mesh->numnp(), &(mesh->ui(0,0)), 2, ind);
    array_attribute( "dep", "positions" );
    end_array();

    // Fields
    for (int i = 0; i < mesh->get_ndf(); ++i) {
        int ind[1] = {i};
        start_array( mesh->get_ndf(), mesh->numnp(), &(mesh->u(0,0)), 1, ind);
        array_attribute( "dep", "positions" );
        end_array();
    }

    // Fields(Imaginary part)
    for (int i = 0; i < mesh->get_ndf(); ++i) {
        int ind[1] = {i};
        start_array( mesh->get_ndf(), mesh->numnp(), &(mesh->ui(0,0)), 1, ind);
        array_attribute( "dep", "positions" );
        end_array();
    }
}


void DXFile::writemesh3d(Mesh* mesh)
{
    start_array( mesh->get_ndm(), mesh->numnp(), &(mesh->x(0,0)) );
    array_attribute( "dep", "positions" );
    end_array();

    int order = order3d(mesh->get_nen());
    int np = order+1;

    start_array( 8, mesh->numelt() * order*order*order, (int*) NULL );
    array_attribute( "element type", "cubes"     );
    array_attribute( "ref",          "positions" );

    // For each element
    for (int i = 0; i < mesh->numelt(); ++i)

        // For each sub-panel of the element
        for (int ii = 0; ii < order; ++ii)
            for (int jj = 0; jj < order; ++jj)
                for (int kk = 0; kk < order; ++kk)

                    // Write the sub-panel connectivity
                    for (int ll = 0; ll < 2; ++ll)
                        for (int mm = 0; mm < 2; ++mm)
                            for (int nn = 0; nn < 2; ++nn) {
                                int j = ((kk+ll)*np+(jj+mm))*np+(ii+nn);
                                write(mesh->ix(j,i));
                            }

    end_array();

    // Displacement vector
    int ind[3] = {0, 1, 2};
    start_array( mesh->get_ndf(), mesh->numnp(), &(mesh->u(0,0)), 3, ind);
    array_attribute( "dep", "positions" );
    end_array();

    // Displacement vector(Imaginary part)
    start_array( mesh->get_ndf(), mesh->numnp(), &(mesh->ui(0,0)), 3, ind);
    array_attribute( "dep", "positions" );
    end_array();

    // Fields
    for (int i = 0; i < mesh->get_ndf(); ++i) {
        int ind[1] = {i};
        start_array( mesh->get_ndf(), mesh->numnp(), &(mesh->u(0,0)), 1, ind);
        array_attribute( "dep", "positions" );
        end_array();
    }

    // Fields(Imaginary part)
    for (int i = 0; i < mesh->get_ndf(); ++i) {
        int ind[1] = {i};
        start_array( mesh->get_ndf(), mesh->numnp(), &(mesh->ui(0,0)), 1, ind);
        array_attribute( "dep", "positions" );
        end_array();
    }
}
