-- Loading the Lua input file
f0 = 20
runfile = loadfile('pml1d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()

-- Assemble distributed stiffness
nev = 2
w0  = 0
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

-- Delete objects
VecDestroy(Vz)
mesh:delete()
