-- Start clock
ts   = os.clock()

-- Loading the Lua input file
f0 =  20
runfile = loadfile('pml1d.lua')
--runfile = loadfile('block3d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print('Numid:',numid)

-- Assemble distributed stiffness
nev     = 2
w0      = 0
Kshiftz = Mesh_assemble_dRz_trilinos(mesh,1,0,-w0*w0,1,0)
Mz      = Mesh_assemble_dRz_trilinos(mesh,0,0,     1,1,0)

-- Construct Operator
Op      = Amesos_Operator_Complex:new(Kshiftz, Mz, 6)
pl      = AnasaziDefaultPL()
dr      = {}
di      = {}
Vz      = Epetra_MultiVector_Complex:create(numid,nev)
status  = compute_eigs_anasazi1(Op, nev, pl, dr, di, Vz);
undo_spectral_trans(nev,w0,0,1,dr,di)
print('Status:',status)

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
pl:delete()
Op:delete()
Vz:delete()
Mz:delete()
Kshiftz:delete()
mesh:delete()

-- Stop clock
te   = os.clock()
print('Elapsed time:',te-ts)
