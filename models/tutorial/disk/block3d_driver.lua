-- Loading the Lua input file
runfile = loadfile('circular_disk3d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
numid = mesh:get_numid()
print('Numid:',numid)

--[[
-- Assemble stiffness
K = mesh:assemble_dR(1,0,0)

-- Assemble the loading vector
QArray_type = 0
F = QArray:new(numid, 1, QArray_type)
mesh:assemble_R()
mesh:get_f(F)
F:mul(-1)

-- Assemble LHS
U = QArray:new(numid, 1, QArray_type)

-- Set up the linear system and solve
K:solve(U, F)
mesh:set_u(U)

-- Check residual
R = QArray:new(numid, 1, QArray_type)
K:apply(U, R)
R:sub(F)
resid = R:normf()
print('Norm of residual:',resid)
--]]

-- Write file for plotting
---[[
dxf = DXFile:new('disk3d')
dxf:writemesh(mesh)
dxf:delete()
--]]

-- Print results
--[[
U:print()
F:print()
--]]
--[[
-- Delete objects
R:delete()
U:delete()
F:delete()
K:delete()
--]]
