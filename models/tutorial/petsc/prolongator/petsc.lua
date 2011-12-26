-- Loading the Lua input file
require 'meshes1.lua'

-- Construct fine(to) mesh
dense_x=  dense_xa[2]
dense_y=  dense_ya[2]
mesh_t = Mesh_Add_Block:new(2)
mesh_t = mesh_func(mesh_t)
mesh_t:initialize()
--mesh_t:initialize_minimal()

-- Serial partitioner mesh_t
meshp_t = Mesh_Partitioner_METIS:new(nparts,mesh_t)
meshp_t:partition()

-- Construct coarse(from) mesh
dense_x=  dense_xa[1]
dense_y=  dense_ya[1]
mesh_f = Mesh_Add_Block:new(2)
mesh_f = mesh_func(mesh_f)
mesh_f:initialize()
--mesh_f:initialize_minimal()

-- Construct mapping between meshes
npart_t= meshp_t:get_node_partition()
nemap  = mesh_t:node_element_mapping(mesh_f)

-- Serial partitioner mesh_f conforming to partition mesh_t
meshp_f = Mesh_Partitioner_Conform:new(nparts,mesh_f)
meshp_f:build_pepairs(mesh_t:numnp(),npart_t,nemap)
meshp_f:partition()

-- Construct prolongator with PETSC
mat_type = 2
is_reduced=1
Pz = Mesh_assemble_P_petsc1(mesh_f,mesh_t,mesh_t:numnp(),nemap,mat_type,is_reduced)

-- print results
--pv = PETSC_VIEWER_STDOUT_WORLD
pv = PetscViewerASCIIOpen('petscmat1.m')
PetscViewerSetFormat(pv,'matlab')
MatView(Pz,pv)
PetscViewerDestroy(pv)

-- clean up
MatDestroy(Pz)
meshp_t:delete()
meshp_f:delete()
mesh_f:delete()
mesh_t:delete()
