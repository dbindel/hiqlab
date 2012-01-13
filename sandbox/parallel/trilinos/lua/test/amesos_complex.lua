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

-- Assemble the loading vector
Fz = Mesh_assemble_Rz_trilinos(mesh,0)
Fz:Scale(-1,0)

-- Assemble LHS
Uz = Epetra_Vector_Complex:create(numid)

-- Set up the linear system 
ELP= Epetra_LinearProblem_Complex:new(Kz,Uz,Fz)

-- Solve problem with Amesos
ABS = Amesos_BaseSolver:new(ELP, 'mumps')
ABS:SymbolicFactorization()
ABS:NumericFactorization()
ABS:Solve()

-- Check error
Rz = Epetra_Vector_Complex:create(numid)
Kz:Multiply(false, Uz, Rz)
Rz:Update( -1, 0, Fz, 1, 0);

resid = {}
Rz:Norm2(resid)
print('Norm of residual:',resid[1])

-- Print results
--[[
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
Rz:delete()
Uz:delete()
Fz:delete()
Kz:delete()
mesh:delete()

-- Stop clock
te   = os.clock()
print('Elapsed time:',te-ts)
