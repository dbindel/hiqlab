-- Loading the Lua input file
require 'test_part_mesh.lua'

rank = MPI_Comm_rank()
size = MPI_Comm_size()

-- Specify parameters
is_removed_unused_nodes = 0
--order = 1;
--dense = 0.25e-6;
--dense_x=  0.25e0
dense_y=  1e0

-- Construct mesh and partition into number of processes
nparts= size
mesh  = Mesh:new(2);

-- Run the Lua input file for serial partitioner
meshp = Parallel_Mesh_Partitioner:new(mesh,nparts)
if (rank==0) then
    mesh = test_part_mesh(mesh)
end

-- Partition 
meshp:tie(1e-8)
meshp:partition_nodes(0)
meshp:partition_elements()
meshp:sendreceive_nodes_info()
meshp:sendreceive_elements_info()

print('Rank[',rank,']:',meshp:num_mymeshes(),'\n')

-- MUST SEE THE INFO
mnp = {}
for i = 1,meshp:num_mymeshes() do
    print('\n\n')
    print('--- Rank', rank, '  Mesh',i,' info ---')
    mnp[i] = meshp:get_mesh(i-1)
    mnp[i] = test_part_mesh(mnp[i])

    mnp[i]:initialize()
    mesh_view(mnp[i])
end

meshp:delete()
mesh:delete()
