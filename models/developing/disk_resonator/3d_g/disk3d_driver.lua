-- Loading the Lua input file
runfile = loadfile('disk3d.lua')

-- Specify parameters
f0    =20;
order = 1;
dense = 1e-6;

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print('Numid:',numid)
print('Numelt:',mesh:numelt())

---[[
-- Compute eignvalues
w0      = 2.889953e+08*2*pi;
nev= 1
ncv= 10
--status, dr, di, vr, vi = compute_eigs_vecs(mesh, w0, nev, ncv)
--status, dr, di, Vz = compute_eigs_slepc(mesh, w0, nev)
status, dr, di, Vz = compute_eigs_anasazi_complex(mesh, w0, nev);

-- Print eigs
for k = 1, nev do
    print(k, ':', dr[k]/2e6/pi, '  :', di[k]/2e6/pi, 'MHz')
end

--]]

-- Write file for plotting
--[[
mesh:set_u(vr)
mesh:set_ui(vi)
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
