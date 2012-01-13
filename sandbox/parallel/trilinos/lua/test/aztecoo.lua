-- Start clock
ts   = os.clock()

-- Loading the Lua input file
f0 = 0
runfile = loadfile('block3d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print('Numid:',numid)

-- Assemble distributed stiffness
K = mesh:assemble_dR_trilinos(1,0,0,1)

-- Assemble RHS
F = mesh:assemble_R_trilinos()
F:Scale(-1)

-- Assemble LHS
U = Epetra_Vector:new(numid)

-- Set up the linear system 
ELP= Epetra_LinearProblem:new(qCrs2RowMatrix(K),U,F)

-- Solve linear system in Anasazi with AztecOO
AZS = AztecOO:new(ELP)
AZS:SetAztecOption(AZ_solver, AZ_gmres)
AZS:SetAztecOption(AZ_output, 32)
AZS:SetAztecOption(AZ_kspace, 160)

-- Construct preconditioner for AztecOO
AZS:SetAztecOption(AZ_precond, AZ_dom_decomp)
AZS:SetAztecOption(AZ_subdomain_solve, AZ_ilu)

-- Solve problem 
AZS:Iterate(1000, 1e-9)
AZS:NumIters()
AZS:TrueResidual()

-- Check residual
R = Epetra_Vector:new(numid)
K:Multiply(false, U, R)
R:Update( -1, F, 1);

resid = {}
R:Norm2(resid)
print('Norm of residual:',resid[1])

-- Print results
--[[
ToMatrixMarketFile('K.mm', K)
ToMatrixMarketFile('U.mm', U)
ToMatrixMarketFile('F.mm', F)
--]]

-- Delete objects
R:delete()
AZS:delete()
ELP:delete()
U:delete()
F:delete()
K:delete()
mesh:delete()

-- Stop clock
te   = os.clock()
print('Elapsed time:',te-ts)
