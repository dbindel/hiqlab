-- Loading the Lua input file
require 'meshes1.lua'

-- Run the Lua input file for serial partitioner
mesh = Mesh:new(ndm);
mesh = mesh_func(mesh)
mesh:initialize()
--mesh:initialize_minimal()

-- Run mesh partitioner
meshp = Mesh_Partitioner_METIS:new(nparts,mesh)
meshp:partition()

-- Write info about mesh
mnpw = Mesh_Partition_Writer:new()
mnpw:write_mesh(meshp,'mesh')
mnpw:delete()

-- Assemble stiffness and dump
Kz = mesh:assemble_dR(1,0,0)
Kz:dump('K1.txt')
Kz:delete()

-- Assemble rhs and dump
mesh:assemble_R()
Fz= QArray:new(mesh:get_numid(),1,2)
mesh:get_f(Fz)
Fz:dump('F1.txt','%16.16f')
Fz:delete()

--[[
print('\n\n')
print('--- Mesh entire info ---')
mesh_view(mesh)
--]]

-- Clean up writing environment
meshp:delete()
mesh:delete()





