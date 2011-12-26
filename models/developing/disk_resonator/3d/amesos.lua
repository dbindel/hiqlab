ts   = os.clock()
-- Loading the Lua input file
f0 = 0
order = 1;
dense = 1e-6; -- 1e-6
runfile = loadfile('disk3d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print('Numid:',numid)

-- Assemble distributed stiffness
K = mesh:assemble_dR_trilinos(1,0,0)

-- Assemble the loading vector
F = mesh:assemble_R_trilinos()
F:Scale(-1)

-- Assemble LHS
U = Epetra_Vector:new(numid)

-- Set up the linear system 
ELP= Epetra_LinearProblem:new(qCrs2RowMatrix(K),U,F)

-- Solve problem with Amesos
ABS = Amesos_BaseSolver:new(ELP, 'mumps')
ABS:SymbolicFactorization()
ABS:NumericFactorization()
ABS:Solve()

-- Check residual
R = Epetra_Vector:new(numid)
K:Multiply(false, U, R)
R:Update( -1, F, 1);

resid = {}
R:Norm2(resid)
print('Norm of residual:',resid[1])

-- Print results
---[[
ToMatrixMarketFile('K.mm', K)
ToMatrixMarketFile('U.mm', U)
ToMatrixMarketFile('F.mm', F)
--]]

-- Delete objects
ABS:delete()
ELP:delete()
R:delete()
U:delete()
F:delete()
K:delete()

te   = os.clock()
print(te-ts)
