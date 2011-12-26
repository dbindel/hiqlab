-- Loading the Lua input file
runfile = loadfile('disk3d.lua')

-- Specify parameters
f0    = 0;
order = 1;
dense = 1e-6;

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
numelt = mesh:numelt()
print('Numid:',numid)
print('Numelt:',numelt)

-- Assemble distributed stiffness
--mattype    = 2;
--is_reduced = 1;
--Kz = mesh:assemble_dR_petsc(1, 0, 0, mattype, is_reduced)

---[[
nev = 1
w0  = 2.889953e+08*2*pi;
status, dr, di, Vz = compute_eigs_slepc(mesh, w0, nev)

-- Print eigenvalues
for i = 1,nev do
   print(dr[i],' ',di[i])
end

-- Print results
--[[
ToMatrixMarketFile('Vr.mm',Vr)
ToMatrixMarketFile('Vi.mm',Vi)
--]]

--]]

-- Delete objects
MatDestroy(Kz)
mesh:delete()
