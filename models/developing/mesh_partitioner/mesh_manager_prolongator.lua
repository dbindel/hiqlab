-- Loading the Lua input file
require 'meshes1.lua'

-- Specify parameters
order = ordera[1];

-- Construct fine(to) mesh
dense_x=  dense_xa[2]
dense_y=  dense_ya[2]
mesh_t = Mesh_Add_Block:new(ndm)
mesh_t = mesh_func(mesh_t)
mesh_t:initialize()
--mesh_t:initialize_minimal()

-- Serial partitioner mesh_t
meshp_t = Mesh_Partitioner_METIS:new(nparts,mesh_t)
meshp_t:partition()
npart_t= meshp_t:get_node_partition()

-- Construct coarse(from) mesh
dense_x=  dense_xa[1]
dense_y=  dense_ya[1]
mesh_f = Mesh_Add_Block:new(ndm)
mesh_f = mesh_func(mesh_f)
mesh_f:initialize()
--mesh_f:initialize_minimal()

-- Construct mapping between meshes
nemap = mesh_t:node_element_mapping(mesh_f)

-- Serial partitioner mesh_f conforming to partition mesh_t
meshp_f = Mesh_Partitioner_Conform:new(nparts,mesh_f)
meshp_f:build_pepairs(mesh_t:numnp(),npart_t,nemap)
meshp_f:partition()

-- Construct mesh_manager and meshes on each process and add to Mesh_Manager
dense_x=  dense_xa[2]
dense_y=  dense_ya[2]
mm_t = Mesh_Manager:new()
mm_t:create_meshes(meshp_t,mesh_func)

-- Construct mesh_manager and meshes on each process and add to Mesh_Manager
dense_x=  dense_xa[1]
dense_y=  dense_ya[1]
mm_f = Mesh_Manager:new()
mm_f:create_meshes(meshp_f,mesh_func)

-- Construct prolongator serially
P1 = assemble_P(mesh_f,mesh_t,mesh_t:numnp(),nemap,1)

-- Construct prolongator parallely
pmm = Prolongator_Mesh_Manager:new(mm_f,mm_t)
pmm:set_nemap(mesh_t:numnp(),nemap)
P   = pmm:assemble_P(1)

-- print matrix
P1:dump('P1.txt')
P:dump('P.txt')

-- clean up
P1:delete()
P:delete()

pmm:delete()
mm_f:delete_meshes()
mm_f:delete()
mm_t:delete_meshes()
mm_t:delete()

meshp_t:delete()
meshp_f:delete()
mesh_f:delete()
mesh_t:delete()
