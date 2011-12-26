-- Loading the Lua input file
require 'mesh1.lua'

-- Specify parameters
order  =  1
dense_x=  2.0e0
dense_y=  1.0e-0

-- Run the Lua input file for single Mesh
mesh = Mesh:new(2);
mesh = mesh_func(mesh)
mesh:initialize()

-- assemble stiffness in CSCMatrix
K1= mesh:assemble_dR(1,0,0)

-- assemble rhs in QArray
mesh:assemble_R()
F1= QArray:new(mesh:get_numid(),1,2)
mesh:get_f(F1)

-- solve system
K1umf = UMFMatrix:new(K1)
U1    = K1umf:solve(F1)

-- Print results
K1:dump('K1.txt')
F1:dump('F1.txt','%16.16f')
U1:dump('U1.txt','%16.16f')

-- clean up
U1:delete()
F1:delete()
K1umf:delete()
K1:delete()

-- save mesh in DXFile
--[[
dxf = DXFile:new('mesh1')
dxf:writemesh(mesh)
dxf:delete()
--]]

-- clean up
mesh:delete()
