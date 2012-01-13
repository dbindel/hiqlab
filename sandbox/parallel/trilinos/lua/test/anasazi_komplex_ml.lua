-- Start clock
ts   = os.clock()

-- Loading the Lua input file
f0 =20;
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
Kk      = Complex2SuperCrsMatrix(Kshiftz, 0)
--ToMatrixMarketFile('K.mm', Kk:GetCrsMatrix())

-- Set AztecOO for Anasazi
ELP= Epetra_LinearProblem:new()
ELP:SetOperator(qCrs2RowMatrix(Kk:GetCrsMatrix()))
AZS = AztecOO:new(ELP)
AZS:SetAztecOption(AZ_solver, AZ_gmres)
AZS:SetAztecOption(AZ_output, 32)
AZS:SetAztecOption(AZ_kspace, 160)

-- Construc ML preconditioner for AztecOO
mlpl = ParameterList:new()
ML_Epetra_SetDefaults(mlpl, "SA")
mlpl:set_int(   "output", 2)
mlpl:set_string("prec type", "full-MGV")
mlpl:set_int(   "print unused",  1)
mlpl:set_int(   "max levels", 10)
mlpl:set_string("increasing or decreasing", "increasing")
mlpl:set_string("aggregation: type", "ParMETIS")
mlpl:set_string("smoother: type", "Aztec")
mlpl:set_double("smoother: damping factor", 0.67)
mlpl:set_int(   "smoother: sweeps", 2)
mlpl:set_string("smoother: pre or post", "both")
mlpl:set_string("coarse: type", "Amesos-KLU")

-- Construct Aztec Options and Params for ML smoother
Aop = AztecOptionsParams:new()
Aop:SetOption(AZ_precond, AZ_dom_decomp)
Aop:SetOption(AZ_subdomain_solve, AZ_ilu)
Aop:SetPL(mlpl)
MLPrec = qCreate_ML_Epetra_MultiLevelPreconditioner(Kk:GetCrsMatrix(),mlpl)
AZS:SetPrecOperator(MLPrec)

-- Construct Operator to pass into Anasazi
AZOp = AztecOO_Operator_Komplex:new(AZS,Mz:OperatorDomainMap(),Mz)
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
Vz:delete()
Aop:delete()
AZOp:delete()
AZS:delete()
mlpl:delete()
MLPrec:delete()
ELP:delete()
Kk:delete()
Kshiftz:delete()
Mz:delete()
mesh:delete()

-- Stop clock
te   = os.clock()
print('Elapsed time:',te-ts)
