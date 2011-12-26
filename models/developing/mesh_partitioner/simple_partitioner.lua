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
npart_f = meshp_f:get_node_partition()


---[[
-- assemble stiffness and rhs of mesh_t
Kt= mesh_t:assemble_dR(1,0,0)
mesh_t:assemble_R()
Ft= QArray:new(mesh_t:get_numid(),1,2)
mesh_t:get_f(Ft)

-- assemble stiffness and rhs of mesh_f
Kf= mesh_f:assemble_dR(1,0,0)
mesh_f:assemble_R()
Ff= QArray:new(mesh_f:get_numid(),1,2)
mesh_f:get_f(Ff)

-- Construct prolongator
Pft = assemble_P(mesh_f,mesh_t,mesh_t:numnp(),nemap,1)

-- print data
Kt:dump('Kt.txt')
Kf:dump('Kf.txt')
Ft:dump('Ft.txt','%16.16f')
Ff:dump('Ff.txt','%16.16f')
Pft:dump('Pft.txt')

--[[
-- write mesh in dxfile
dxf = DXFile:new('mesh_f')
dxf:writemesh(mesh_f)
dxf:delete()
dxf = DXFile:new('mesh_t')
dxf:writemesh(mesh_t)
dxf:delete()
--]]

--[[
-- Write out to partition info
mnpw = Mesh_Partition_Writer:new()
mnpw:write_mesh(meshp_t,'mesh_t')
mnpw:write_mesh(meshp_f,'mesh_f')
mnpw:write_nemap1(meshp_t,mesh_t:numnp(),nemap,'Pft')
mnpw:delete()
--]]

-- clean up
Kt:delete()
Kf:delete()
Ft:delete()
Ff:delete()
Pft:delete()

meshp_t:delete()
meshp_f:delete()
mesh_f:delete()
mesh_t:delete()
