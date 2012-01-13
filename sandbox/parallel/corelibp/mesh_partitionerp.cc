#include <cassert>
#include <iostream>

#include "mpi.h"

#include "mesh_partitionerp.h"

#define ME Parallel_Mesh_Partitioner

ME::ME(Mesh* mesh, int nparts) : mesh(mesh), nparts(nparts)
{

    MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocess);
std::cout << "Process:" << process_id << "/" << numprocess << "\n";

    // -- Distribute partitions onto processes
    assert(nparts>=numprocess);
    PP.resize(numprocess+1);
    PP[0] = 0;
    int npp = nparts/numprocess;
    int res = nparts - npp*numprocess;
    for (int i = 0; i < numprocess; ++i)
        PP[i+1] = npp;
    for (int i = 0; i < res; ++i)
        PP[i+1] += 1;
    for (int i = 0; i < numprocess; ++i)
        PP[i+1]+=PP[i];

    // -- Create necessary Mesh_Partition on each processor
    for (int i = 0; i < PP[process_id+1]-PP[process_id]; ++i) {

        // -- Construct Mesh_Partition
        Mesh_Partition* mnp = new Mesh_Partition(mesh->get_ndm());
        // -- Set partition id and lua state
        mnp->set_pid(i);
        mnp->set_lua(mesh->get_lua());

        meshes.push_back(mnp);
    }

    // -- Create partitioner on root process
    if (process_id==0)
        meshp = new Mesh_Partitioner(mesh,nparts);
}


ME::~ME()
{
    if (process_id==0)
        delete meshp;
    for (int i = 0; i < PP[process_id+1]-PP[process_id]; ++i)
        delete meshes[i];
}


void ME::tie(double tol, int start, int end)
{
    if (process_id==0)
        meshp->tie(tol,start,end);
}


void ME::partition_nodes(int ptype)
{
    if (process_id==0)
        meshp->partition_nodes(ptype);
}


void ME::partition_elements()
{
    if (process_id==0)
        meshp->partition_elements();
}

void ME::send_nodes_info(int sproc, int numpart)
{
    int tag = 0;
    MPI_Request requests[3];
    MPI_Status  statuss[3];

    // -- Construct arrays to send to each process
    int snumnp = meshp->get_numnp(numpart);
    int snodes[snumnp];
    int sxmap [snumnp];
    int sparts[snumnp];
    meshp->get_nodes(numpart,snodes,sxmap,sparts);

    // -- Non-blocking send but wait before exiting
    tag = 0;
    MPI_Isend(snodes, snumnp, MPI_INT, sproc, tag, MPI_COMM_WORLD,&requests[0]);
    tag = 1;
    MPI_Isend(sxmap , snumnp, MPI_INT, sproc, tag, MPI_COMM_WORLD,&requests[1]);
    tag = 2;
    MPI_Isend(sparts, snumnp, MPI_INT, sproc, tag, MPI_COMM_WORLD,&requests[2]);
    MPI_Waitall(3,requests,statuss);
}


void ME::receive_nodes_info(int rproc, int numpart)
{
    // -- Check if receiving info is for this process
    assert(numpart>=PP[process_id] && numpart <PP[process_id+1]);

    int tag = 0;
    MPI_Request requests[3];
    MPI_Status  statuss[3];
    MPI_Status  status;

    // -- Blocking Probe to find out the length of the receiving array
    int rnumnp;
    tag = 0;
    MPI_Probe(rproc,tag,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status,MPI_INT,&rnumnp);

    // -- Construct arrays to receive from each process
    int rnodes[rnumnp];
    int rxmap [rnumnp];
    int rparts[rnumnp];

    // -- Non-blocking recieve but wait before exiting
    tag = 0;
    MPI_Irecv(rnodes,rnumnp,MPI_INT,rproc,tag,MPI_COMM_WORLD,&requests[0]);
    tag = 1;
    MPI_Irecv(rxmap ,rnumnp,MPI_INT,rproc,tag,MPI_COMM_WORLD,&requests[1]);
    tag = 2;
    MPI_Irecv(rparts,rnumnp,MPI_INT,rproc,tag,MPI_COMM_WORLD,&requests[2]);
    MPI_Waitall(3,requests,statuss);

    // -- Set data on a partition
    meshes[numpart-PP[process_id]]->set_nodes(rnumnp,rnodes);
    meshes[numpart-PP[process_id]]->set_ntie_map(rnumnp,rxmap);
    meshes[numpart-PP[process_id]]->set_npart_map(rnumnp,rparts);
}


void ME::sendreceive_nodes_info_root()
{
    for (int i = 0; i < PP[1]; ++i) {

        int numpart = i;

        // -- Get nodes
        int snumnp = meshp->get_numnp(numpart);
        int snodes[snumnp];
        int sxmap [snumnp];
        int sparts[snumnp];
        meshp->get_nodes(numpart,snodes,sxmap,sparts);

        // -- Set nodes
        meshes[numpart-PP[process_id]]->set_nodes(snumnp,snodes);
        meshes[numpart-PP[process_id]]->set_ntie_map(snumnp,sxmap);
        meshes[numpart-PP[process_id]]->set_npart_map(snumnp,sparts);
    }
}


void ME::sendreceive_nodes_info()
{
    if (process_id==0) {

        // -- send/receive info for root
        sendreceive_nodes_info_root();

        // -- send info for other processes
        for (int i = 1; i < numprocess; ++i) {
            for (int j = PP[i]; j < PP[i+1]; ++j)
                send_nodes_info(i,j);
        }

    } else
        for (int j = PP[process_id]; j < PP[process_id+1]; ++j)
            receive_nodes_info(0,j);
}

void ME::send_elements_info(int sproc, int numpart)
{
    int tag = 0;
    MPI_Request requests[2];
    MPI_Status  statuss[2];

    // -- Construct arrays to send to each process
    int snumelt = meshp->get_numelt(numpart);
    int selts[snumelt];
    int sparts[snumelt];
    meshp->get_elements(numpart,selts,sparts);

    // -- Non-blocking send but wait before exiting
    tag = 0;
    MPI_Isend(selts,  snumelt, MPI_INT, sproc, tag, MPI_COMM_WORLD,&requests[0]);
    tag = 1;
    MPI_Isend(sparts, snumelt, MPI_INT, sproc, tag, MPI_COMM_WORLD,&requests[1]);
    MPI_Waitall(2,requests,statuss);
}


void ME::receive_elements_info(int rproc, int numpart)
{
    // -- Check if receiving info is for this process
    assert(numpart>=PP[process_id] && numpart <PP[process_id+1]);

    int tag = 0;
    MPI_Request requests[2];
    MPI_Status  statuss[2];
    MPI_Status  status;

    // -- Blocking Probe to find out the length of the receiving array
    int rnumelt;
    tag = 0;
    MPI_Probe(rproc,tag,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status,MPI_INT,&rnumelt);

    // -- Construct arrays to receive from each process
    int relts [rnumelt];
    int rparts[rnumelt];

    // -- Non-blocking recieve but wait before exiting
    tag = 0;
    MPI_Irecv(relts ,rnumelt,MPI_INT,rproc,tag,MPI_COMM_WORLD,&requests[0]);
    tag = 1;
    MPI_Irecv(rparts,rnumelt,MPI_INT,rproc,tag,MPI_COMM_WORLD,&requests[1]);
    MPI_Waitall(2,requests,statuss);

    // -- Set data on a partition
    meshes[numpart-PP[process_id]]->set_elements(rnumelt,relts);
    meshes[numpart-PP[process_id]]->set_epart_map(rnumelt,rparts);
}


void ME::sendreceive_elements_info_root()
{
    for (int i = 0; i < PP[1]; ++i) {

        int numpart = i;

        // -- Get elements
        int snumelt = meshp->get_numelt(numpart);
        int selts [snumelt];
        int sparts[snumelt];
        meshp->get_elements(numpart,selts,sparts);

        // -- Set elements
        meshes[numpart-PP[process_id]]->set_elements(snumelt,selts);
        meshes[numpart-PP[process_id]]->set_epart_map(snumelt,sparts);
    }
}


void ME::sendreceive_elements_info()
{
    if (process_id==0) {

        // -- send/receive info for root
        sendreceive_elements_info_root();

        // -- send info for other processes
        for (int i = 1; i < numprocess; ++i) {
            for (int j = PP[i]; j < PP[i+1]; ++j)
                send_elements_info(i,j);
        }

    } else
        for (int j = PP[process_id]; j < PP[process_id+1]; ++j)
            receive_elements_info(0,j);
}

void ME::allgather_idg_info()
{
    IDP.resize(numprocess+1);
    int ptnumid = 0;
    int nmeshes = meshes.size();

    // -- Compute the number of ids on this process
    for (int i = 0; i < nmeshes; ++i)
        ptnumid += meshes[i]->potential_numid();
    IDP[process_id] = ptnumid;

    // -- allgather info
    MPI_Allgather(&(IDP[process_id]),1,MPI_INT,&(IDP[0]),1,MPI_INT,MPI_COMM_WORLD);
}

void ME::assign_onp_idg()
{
    int nmeshes = meshes.size();
    int ptnumid = IDP[process_id];
    for (int i = 0; i < nmeshes; ++i) {
        meshes[i]->number_own_id(ptnumid, ptnumid+meshes[i]->potential_numid());
        ptnumid+=meshes[i]->potential_numid();
    }
}

void ME::assign_offp_idg()
{
}
