-- Loading the Lua input file
f0 = 00
--order=1
--dense=5
--runfile = loadfile('/home/tkoyama/programs/hiqlab/models/developing/disk_resonator/3d/plate3d.lua')
--runfile = loadfile('block3d.lua')
runfile = loadfile('pml1d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print('Numid:',numid)

-- Assemble distributed stiffness
nev = 2
w0  = 0
status, dr, di, Vz = compute_eigs_arpack(mesh, w0, nev);

-- Print eigenvalues
for i = 1,nev do
   print(dr[i],' ',di[i])
end

-- Print results
--[[
Kz = Mesh_assemble_dRz_trilinos(mesh,1,0,0,1,0)
Mz = Mesh_assemble_dRz_trilinos(mesh,0,0,1,1,0)
ToMatrixMarketFile('Mr.mm',Mz)
ToMatrixMarketFile('Mi.mm',Mz:get_Az())
ToMatrixMarketFile('Kr.mm',Kz)
ToMatrixMarketFile('Ki.mm',Kz:get_Az())
ToMatrixMarketFile('Vr.mm',Vz)
ToMatrixMarketFile('Vi.mm',Vz:get_Vz())
--]]

-- Delete objects
Vz:delete()
mesh:delete()
