-- Loading the Lua input file
require 'meshes1.lua'

-- Specify parameters
order   = ordera[3]
dense_x = dense_xa[3]
dense_y = dense_ya[3]

-- Run the Lua input file for serial partitioner
mesh = Mesh:new(ndm);
mesh = mesh_func(mesh)
mesh:initialize()
--mesh:initialize_minimal()

-- Run mesh partitioner
meshp = Mesh_Partitioner_METIS:new(nparts,mesh)
meshp:partition()

--[[
print('\n\n')
print('--- Mesh entire info ---')
mesh_view(mesh)
--]]

-- Write info about mesh
mnpw = Mesh_Partition_Writer:new()
mnpw:write_mesh(meshp,fnamea[3])
mnpw:delete()

-- Clean up writing environment
meshp:delete()
mesh:delete()
