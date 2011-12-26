-- Loading the Lua input file
runfile = loadfile('disk3d.lua')

-- Specify parameters
f0    = 0;
order = 1;
dense = 1e-6; -- 1e-6
--xo    = 0.2e-6
--rbd   = 6e-6;
--rpml  = 8e-6;

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print('Numid:',numid)

---[[
-- Compute eignvalues
--w0      = 7.406e+08*2*pi;
--w0 = 755.7050479892e6*2*pi;
--w0 = 765.7050479892e6*2*pi;
--w0 = 765.7050479892e6*2*pi;
--w0 = 766.97516471285e6*2*pi;
w0      = 2.889953e+08*2*pi;
nev= 1 
ncv= 5
--dr = {}
--di = {}
status, dr, di, vr, vi = compute_eigs_vecs(mesh, w0, nev, ncv)
--status, dr, di, Vz = compute_eigs_anasazi_complex(mesh, w0, nev);
--status, dr, di, Vz = compute_eigs_arpack(mesh, w0, nev);
--mesh:set_u(vr)
--mesh:set_ui(vi)

-- Print eigs
for k = 1, nev do
    local Q = sqrt(dr[k]^2+di[k]^2)/2/di[k]
    print(k, ':', dr[k]/2e6/pi, '  :', di[k]/2e6/pi, 'MHz',  '  Q:',Q)
end

--]]

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
--[[
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

mesh:delete()
