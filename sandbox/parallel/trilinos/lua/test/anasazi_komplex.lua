-- Start clock
ts   = os.clock()

-- Loading the Lua input file
f0 =20;
runfile = loadfile('pml1d.lua')

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
Kk      = Complex2SuperCrsMatrix(Kshiftz, 0)

-- Set AztecOO for Anasazi
ELP= Epetra_LinearProblem:new()
ELP:SetOperator(qCrs2RowMatrix(Kk:GetCrsMatrix()))
AZS = AztecOO:new(ELP)
AZS:SetAztecOption(AZ_solver, AZ_gmres)
AZS:SetAztecOption(AZ_output, 32)
AZS:SetAztecOption(AZ_kspace, 160)

-- Construct preconditioner for AztecOO
AZS:SetAztecOption(AZ_precond, AZ_dom_decomp)
AZS:SetAztecOption(AZ_subdomain_solve, AZ_ilu)

-- Construct Operator to pass into Anasazi
AZOp = AztecOO_Operator_Komplex:new(AZS,Kshiftz:OperatorDomainMap(),Mz)
AZOp:SetMaxIter(1000)
AZOp:SetTol(1e-9)

-- Compute with Anasazi
pl      = AnasaziDefaultPL()
pl:set_int("Verbosity", Anasazi_Warning)
dr      = {}
di      = {}
Vz      = Epetra_MultiVector_Complex:create(numid,nev)
status  = compute_eigs_anasazi1(AZOp, nev, pl, dr, di, Vz);
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
AZOp:delete()
AZS:delete()
ELP:delete()
Kk:delete()
Kshiftz:delete()
Mz:delete()
Vz:delete()
mesh:delete()

-- Stop clock
te   = os.clock()
print('Elapsed time:',te-ts)
