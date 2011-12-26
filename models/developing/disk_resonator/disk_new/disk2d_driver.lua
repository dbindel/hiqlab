-- Loading the Lua input file
runfile = loadfile('diskmesh.lua')

-- Specify parameters
f0    =20;
order = 1;
dense = 0.25e-6;
rbd   = 6e-6;
rpml  = 8e-6;

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print('Numid:',numid)
print('Numelt:',mesh:numelt())

-- Compute eignvalues
w0      = 7.157e+08*2*pi;
--w0      = 2.889953e+08*2*pi;
nev= 1
ncv= 10
--dr = {}
--di = {}
status, dr, di, vr, vi = compute_eigs_vecs(mesh, w0, nev, ncv)
--status, dr, di, Vz = compute_eigs_slepc(mesh, w0, nev)

-- Assemble distributed stiffness
mattype    = 0;
is_reduced = 1;
--Kz = mesh:assemble_dR_petsc(1, 0, 0, mattype, is_reduced)

-- Print eigs
for k = 1, nev do
    print(k, ':', dr[k]/2e6/pi, '  :', di[k]/2e6/pi, 'MHz')
end

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
mesh:set_u(vr)
mesh:set_ui(vi)
dxf = DXFile:new('disk2d')
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

mesh:delete()
