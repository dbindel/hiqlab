-- Loading the Lua input file
require 'meshes1.lua'

-- Construct fine(to) mesh
dense_x=  dense_xa[2]
dense_y=  dense_ya[2]
mesh_t = Mesh_Add_Block:new(2)
mesh_t = mesh_func(mesh_t)
--mesh_t:initialize()
mesh_t:initialize_minimal()

-- Serial partitioner for mesh_t
-- For the serial case there is no requirement to partition
meshp_t = Mesh_Partitioner_METIS:new(nparts,mesh_t)
meshp_t:partition()

-- Construct coarse(from) mesh
dense_x=  dense_xa[1]
dense_y=  dense_ya[1]
mesh_f = Mesh_Add_Block:new(2)
mesh_f = mesh_func(mesh_f)
--mesh_f:initialize()
mesh_f:initialize_minimal()

-- Construct mapping between meshes
npart_t= meshp_t:get_node_partition()
nemap  = mesh_t:node_element_mapping(mesh_f)

-- Serial partitioner mesh_f conforming to partition mesh_t
-- For the serial case there is no requirement to partition
meshp_f = Mesh_Partitioner_Conform:new(nparts,mesh_f)
meshp_f:build_pepairs(mesh_t:numnp(),npart_t,nemap)
meshp_f:partition()

-- Construct prolongator serially
P = assemble_P(mesh_f,mesh_t,mesh_t:numnp(),nemap,1)
P:dump('P1.txt')
P:delete()

-- cleanup
meshp_t:delete()
meshp_f:delete()
mesh_f:delete()
mesh_t:delete()
