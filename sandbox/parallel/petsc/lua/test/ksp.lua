require '../petsc-util.lua'

-- Loading the Lua input file
f0 = 0;
runfile = loadfile('pml1d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()

-- Assemble distributed stiffness
mattype    = 2;
is_reduced = 1;
Kz = mesh:assemble_dR_petsc(1, 0, 0, mattype, is_reduced)

-- Assemble the loading vector
Fz = mesh:assemble_R_petsc(is_reduced)
VecScale(Fz,-1,0)

-- Assemble LHS
Uz = VecCreate()
VecSetSizes(Uz, PETSC_DECIDE, numid)
VecSetFromOptions(Uz)

-- Create KSP
ksp = KSPCreate()
KSPSetOperators(ksp, Kz, Kz, "DIFFERENT_NONZERO_PATTERN")
KSPSetTolerances(ksp, 1.0e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT)
KSPSetFromOptions(ksp)
--pc  = KSPGetPC(ksp)
--PCSetType(pc, 0)

-- Solve system
KSPSolve(ksp,Fz,Uz)

--petsc_vector_dump(Uz,'tt.txt')

rnorm = KSPGetResidualNorm(ksp)
its   = KSPGetIterationNumber(ksp)

-- Delete objects
KSPDestroy(ksp)
VecDestroy(Uz)
VecDestroy(Fz)
MatDestroy(Kz)
mesh:delete()

print('Residual:',rnorm, ' NumIterations:', its)
--[[
--]]
