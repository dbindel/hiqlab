#include "mesh_partition_util.h"
#include <string>
#include <cstdio>
#include <sstream>
#include <cassert>

using std::string;

/** local utility functions
 */
inline static string stringify(int x)
{
    std::ostringstream o;
    o << x;
    return o.str();
}


/** Mesh_Partition_Writer methods
 */
#define ME Mesh_Partition_Writer

ME::~ME()
{
}


void ME::open_file(int ni, const char* fname)
{
    string fnamea = string(fname) + stringify(ni) + string(".hq");
    string fnameb = string(fname) + stringify(ni) + string(".bin");

    if ( (fpa = fopen(fnamea.c_str(),"w")) ==NULL ) {
        printf("Failed in opening file:%s\n",fnamea.c_str());
        exit(1);
    } 

    if ( (fpb = fopen(fnameb.c_str(),"wb")) ==NULL ) {
        printf("Failed in opening file:%s\n",fnameb.c_str());
        exit(1);
    } 

    write_basic_info(fnameb);
}


void ME::close_file()
{
    fclose(fpa);
    fclose(fpb);
}


void ME::write_mesh(Mesh_Partitioner* mp, const char* fname)
{
    int nparts = mp->get_num_partitions();
    for (int i = 0; i < nparts; ++i)
        write_partition(mp,i,fname);
}


void ME::write_partition(Mesh_Partitioner* mp, int ni, const char* fname)
{
    meshp = mp;    
    open_file(ni,fname);

    bytecount = 0;
    write_node_info(ni);
    write_element_info(ni);
    write_basicid_info(ni);
    write_nodeid_info(ni);

    close_file();
}


void ME::write_nemap(Mesh_Partitioner* mp, int* nemap, const char* fname)
{
    int nparts = mp->get_num_partitions();
    for (int i = 0; i < nparts; ++i)
        write_partition_nemap(mp,i,nemap,fname);
}


void ME::write_partition_nemap(Mesh_Partitioner* mp, int ni, int* nemap, const char* fname)
{
    meshp = mp;    
    open_file(ni,fname);

    bytecount = 0;
    write_nemap(ni,nemap);

    close_file();
}


void ME::write_basic_info(string fnameb)
{
    fprintf(fpa,"ndm: %d ndf: %d nen: %d filename: %s\n",
            meshp->get_ndm(),meshp->get_ndf(),meshp->get_nen(),fnameb.c_str());
}


void ME::write_basicid_info(int ni)
{
    fprintf(fpa,"startid: %d end1id: %d\n",
            meshp->get_startid(ni),meshp->get_end1id(ni));
}


void ME::write_node_info(int ni)
{
    int *nodes, *xmap, *parts;
    int numnp = meshp->get_numnp(ni);

    // -- read and write nodes
    nodes  = new int[numnp];
    meshp->get_nodes  (ni,numnp,nodes);
    write_array("nodes",1,numnp,nodes);

    // -- read, write, and delete xmap
    xmap   = new int[numnp];
    meshp->get_tie_map(ni,numnp,nodes,xmap);
    write_array("xmap" ,1,numnp,xmap);
    delete[] xmap;

    // -- read, write, and delete parts
    parts  = new int[numnp];
    meshp->get_node_partition(ni,numnp,nodes,parts);
    write_array("npart",1,numnp,parts);
    delete[] parts;

    // -- delete nodes and clean up
    delete[] nodes;
}


void ME::write_nodeid_info(int ni)
{
    int *nodes, *nodeids;
    int ndf   = meshp->get_ndf();
    int numnp = meshp->get_numnp(ni);
    nodes  = new int[numnp];
    nodeids= new int[numnp*ndf];

    meshp->get_nodes   (ni,numnp,nodes);
    meshp->get_node_ids(ni,numnp,nodes,nodeids);

    write_array("nodeids",ndf,numnp,nodeids);

    // -- clean up
    delete[] nodes;
    delete[] nodeids;
}


void ME::write_element_info(int ni)
{
    int *elts, *parts;
    int numelt = meshp->get_numelt(ni);
    elts   = new int[numelt];
    parts  = new int[numelt];

    meshp->get_elements         (ni,numelt,elts);
    meshp->get_element_partition(ni,numelt,elts,parts);

    write_array("elts", 1,numelt,elts);
    write_array("epart",1,numelt,parts);

    // -- clean up
    delete[] elts;
    delete[] parts;
}


void ME::write_array(const char* aname, int m, int n, int* data)
{
    fprintf(fpa,"%s %d %d %d\n",aname,m,n,bytecount);
    if( (fwrite(data, sizeof(int), m*n, fpb)!=m*n) ) {
        printf("Could not successfully write data [%s] of length [%d]\n",aname,m*n);
        fclose(fpa);
        fclose(fpb);
        exit(1);
    }
    bytecount += m*n*sizeof(int);
}


void ME::write_nemap(int ni, int* nemap)
{
    int *nodes, *pnemap;
    int numnp = meshp->get_numnp(ni);
    nodes  = new int[numnp];
    pnemap = new int[numnp];

    meshp->get_nodes(ni,numnp,nodes);
    for (int i = 0; i < numnp; ++i)
        pnemap[i] = nemap[nodes[i]];

    write_array("nemap", 1, numnp, pnemap);

    // -- clean up
    delete[] nodes;
    delete[] pnemap;
}
#undef ME


/** Mesh_Partition_Reader methods
 */
#define ME Mesh_Partition_Reader

ME::~ME()
{
}


void ME::open_file(int ni, const char* fname)
{
    string fnamea = string(fname) + stringify(ni) + string(".hq");
    string fnameb = string(fname) + stringify(ni) + string(".bin");

    if ( (fpa = fopen(fnamea.c_str(),"r")) ==NULL ) {
        printf("Failed in opening file:%s\n",fnamea.c_str());
        exit(1);
    } 

    if ( (fpb = fopen(fnameb.c_str(),"rb")) ==NULL ) {
        printf("Failed in opening file:%s\n",fnameb.c_str());
        exit(1);
    } 

    read_basic_info(fnameb);
}


void ME::close_file()
{
    fclose(fpa);
    fclose(fpb);
}


void ME::read_partition(Mesh_Partition* mesh, const char* fname)
{
    meshnp = mesh;
    open_file(mesh->get_pid(),fname);

    read_node_info();
    read_element_info();

    close_file();
}


int* ME::read_partition_nemap(Mesh_Partition* mesh, const char* fname)
{
    meshnp = mesh;
    open_file(mesh->get_pid(),fname);

    int numnp, m;
    int* nemap = read_array("nemap", &m, &numnp);

    close_file();

    return nemap;
}


void ME::read_basic_info(string fnameb)
{
    char ndm_s[5], ndf_s[5], nen_s[5], fname_s[10];
    char fnameb_s[100];

    if( fscanf(fpa, "%s   %d %s   %d %s   %d %s        %s", 
        ndm_s, &NDM, ndf_s, &NDF, nen_s, &NEN, fname_s, fnameb_s)!=8) {
        printf("Failed to read basic info from file[%s]\n",fnameb.c_str());
        fclose(fpa);
        fclose(fpb);
        exit(1);
    }

    assert(strcmp(fnameb.c_str(),fnameb_s)==0);
}


void ME::read_node_info()
{
    int *nodes, *xmap, *parts, *nodeids;
    int numnp, m;

    nodes  = read_array("nodes",     &m, &numnp);
    xmap   = read_array("xmap" ,     &m, &numnp);
    parts  = read_array("npart",     &m, &numnp);

    meshnp->set_nodes(numnp,nodes);
    meshnp->set_ntie_map(numnp,xmap);
    meshnp->set_npart_map(numnp,parts);

    // -- clean up
    delete[] nodes;
    delete[] xmap;
    delete[] parts;
}


void ME::read_element_info()
{
    int *elts, *parts;
    int numelt, m;

    elts   = read_array("elts" ,     &m, &numelt);
    parts  = read_array("epart",     &m, &numelt);

    meshnp->set_elements(numelt,elts);
    meshnp->set_epart_map(numelt,parts);

    // -- clean up
    delete[] elts;
    delete[] parts;
}


int* ME::read_array(const char* aname, int* m, int* n)
{
    int sbyte;
    char aname_s[100];
    int* data;

    fscanf(fpa,"%s %d %d %d",aname_s,m,n,&sbyte);
    assert(strcmp(aname,aname_s)==0);

    data = new int[(*m)*(*n)];
    fseek(fpb, sbyte, SEEK_SET);
    if( (fread(data, sizeof(int), (*m)*(*n), fpb)!=(*m)*(*n)) ) {
        printf("Could not successfully read data [%s] of length [%d]\n",aname,(*m)*(*n));
        fclose(fpa);
        fclose(fpb);
        exit(1);
    }

    return data;
}


void ME::read_partition_ids(Mesh_Partition* mesh, const char* fname)
{
    meshnp = mesh;
    open_file(mesh->get_pid(),fname);


    read_basicid_info();
    read_nodeid_info();

    close_file();
}


void ME::read_basicid_info()
{
    // -- Skip 5 times till id info
    int sbyte, m, n;
    char aname_s[100];
    for (int i = 0; i < 5; ++i)
    fscanf(fpa,"%s %d %d %d",aname_s,&m,&n,&sbyte);

    // -- obtain basic id info
    char startid_s[9], end1id_s[8];
    fscanf(fpa,"%s %d %s %d",
            startid_s, &startid, end1id_s, &end1id);
}


void ME::read_nodeid_info()
{
    int *nodeids;
    int numnp, m;

    nodeids = read_array("nodeids",   &m, &numnp);

    meshnp->set_ids_range(startid,end1id);
    meshnp->set_node_ids(numnp, NDF, nodeids);

    // -- clean up
    delete[] nodeids;
}
