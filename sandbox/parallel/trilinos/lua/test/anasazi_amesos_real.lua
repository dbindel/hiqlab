-- Start clock
ts   = os.clock()

-- Loading the Lua input file
f0 =  00
runfile = loadfile('pml1d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print('Numid:',numid)

-- Assemble distributed stiffness
nev    = 2
w0     = 0
Kshift = Mesh_assemble_dR_trilinos(mesh,1,0,-w0*w0,1)
M      = Mesh_assemble_dR_trilinos(mesh,0,0,     1,1)

-- Construct Operator
Op      = Amesos_Operator:new(qCrs2RowMatrix(Kshift), 
                              qCrs2RowMatrix(M), 1) -- KLU
pl      = AnasaziDefaultPL()
dr      = {}
di      = {}
V       = Epetra_MultiVector:new(numid,nev)
status  = compute_eigs_anasazi1(Op, nev, pl, dr, di, V);
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
V:delete()
M:delete()
Kshift:delete()
mesh:delete()

-- Stop clock
te   = os.clock()
print('Elapsed time:',te-ts)
