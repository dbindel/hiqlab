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

-- Assemble distributed stiffness
Kz = Mesh_assemble_dRz_trilinos(mesh,1,0,0,1,0)
Kk = Complex2SuperCrsMatrix(Kz, 0)

-- Assemble the loading vector
Fz = Mesh_assemble_Rz_trilinos(mesh,0)
Fz:Scale(-1,0)
Fzm= qSingle2MultiComplexVector(Fz)
Fk = Complex2SuperMultiVector(Fz, 0)

-- Assemble LHS
Uz = Epetra_MultiVector_Complex:create(numid)
Uk = Epetra_Vector:new(numid*2)

-- Set up the linear system 
ELP= Epetra_LinearProblem:new(qCrs2RowMatrix(Kk:GetCrsMatrix()),Uk,Fk:GetMultiVector())

-- Solve problem with Amesos
ABS = Amesos_BaseSolver:new(ELP, 'mumps')
ABS:SymbolicFactorization()
ABS:NumericFactorization()
ABS:Solve()
Multi2ComplexMultiVector(Uk, Uz, 0)

-- Check error
Rz = Epetra_MultiVector_Complex:create(numid)
Kz:Multiply(false, Uz, Rz)
Rz:Update( -1, 0, Fzm, 1, 0);

Rk = Epetra_Vector:new(numid*2)
Kk:GetCrsMatrix():Multiply(false, Uk, Rk)
Rk:Update( -1, Fk:GetMultiVector(), 1);

resid = {}
Rz:Norm2(resid)
print('Norm of residual:',resid[1])

-- Print results
--[[
ToMatrixMarketFile('Kk.mm', Kk:GetCrsMatrix())
ToMatrixMarketFile('Uk.mm', Uk)
ToMatrixMarketFile('Fk.mm', Fk:GetMultiVector())
ToMatrixMarketFile('Kr.mm', Kz)
ToMatrixMarketFile('Ki.mm', Kz:get_Az())
ToMatrixMarketFile('Ur.mm', Uz)
ToMatrixMarketFile('Ui.mm', Uz:get_Vz())
ToMatrixMarketFile('Fr.mm', Fz)
ToMatrixMarketFile('Fi.mm', Fz:get_Vz())
--]]

-- Delete objects
ABS:delete()
ELP:delete()
Rk:delete()
Uk:delete()
Fk:delete()
Kk:delete()
Rz:delete()
Uz:delete()
Fz:delete()
Fzm:delete()
Kz:delete()
mesh:delete()

-- Stop clock
te   = os.clock()
print('Elapsed time:',te-ts)
