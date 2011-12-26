-- Loading the Lua input file
require 'meshes1.lua'

-- Set up to read partition
mnpr  = Mesh_Partition_Reader:new()
nmesh = MPI_Comm_rank(PETSC_COMM_WORLD);

-- Construct meshes(numbered coarse to fine)
meshes={}
for i = 1,nlevels do
    meshes[i] = Mesh_Partition:new(nmesh,ndm)
    dense_x = dense_xa[i]
    dense_y = dense_ya[i]
    meshes[i]:set_partition(mnpr,fnamea[i],mesh_func)
end

-- Assemble prolongators
mat_type = 2
is_reduced=1
P = {}
for i = 1,nlevels-1 do
    P[i] = Mesh_Partition_assemble_P_petsc2(meshes[i],meshes[i+1],fpnamea[i],mat_type,is_reduced)
end

-- print matrices
--pv = PETSC_VIEWER_STDOUT_WORLD
pv = PetscViewerASCIIOpen('petscmat2.m')
PetscViewerSetFormat(pv,'matlab')
for i = 1,nlevels-1 do
    MatView(P[i],pv)
end
PetscViewerDestroy(pv)


-- clean up
for i = 1,nlevels-1 do
    MatDestroy(P[i])
end

for i = 1,nlevels do
    meshes[i]:delete()
end
mnpr:delete()
