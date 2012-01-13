ts   = os.clock()
require '../petsc-util.lua'

-- Loading the Lua input file
runfile = loadfile('../../lua_trilinos/test/block3d.lua')
--runfile = loadfile('pml1d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print("Numid:",numid)

-- Assemble distributed stiffness
mattype    = 3;
is_reduced = 2;
Kz = mesh:assemble_dR_petsc(1, 0, 0, mattype, is_reduced)

-- Assemble the RHS
Fz = mesh:assemble_R_petsc(is_reduced)
VecScale(Fz,-1)

-- Assemble LHS
Uz = MatGetVecX(Kz)
VecSetFromOptions(Uz)

-- Create KSP
ksp = KSPCreate()
KSPSetOperators(ksp, Kz, Kz, "DIFFERENT_NONZERO_PATTERN")
KSPSetTolerances(ksp, 1.0e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT)
PetscOptionsSetValue("-ksp_type","gmres")
pc  = KSPGetPC(ksp)

-- HYPRE PC
---[[
PCSetType(pc, "hypre")
PCHYPRESetType(pc, "boomeramg")
PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter","2")
PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold","0.25")
PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_all","symmetric-SOR/Jacobi")
PCSetFromOptions(pc)
--]]

-- PROMETHEUS PC
--[[
PCSetType(pc, "prometheus")
PCSetCoordinatesFromMesh(pc,mesh)
PetscOptionsSetValue("-ksp_monitor","")
--PetscOptionsSetValue("-log_summary","")
--PetscOptionsSetValue("-out_verbose","2") 
PetscOptionsSetValue("-pc_mg_type","multiplicative")
PetscOptionsSetValue("-pc_mg_cycles","1")  
PetscOptionsSetValue("-prometheus_preduce_base","500") 
PetscOptionsSetValue("-prometheus_top_grid_limit","2500")
PetscOptionsSetValue("-aggmg_avoid_resmooth","")
PetscOptionsSetValue("-aggmg_smooths","1")
PetscOptionsSetValue("-prometheus_repartition","")
--PetscOptionsSetValue("-options_left","")
PCSetFromOptions(pc)
--]]
PetscOptionsPrint()

-- Solve system
KSPSetFromOptions(ksp)
--KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD)
KSPSolve(ksp,Fz,Uz)

--petsc_vector_dump(Uz,'tt2.txt')

rnorm = KSPGetResidualNorm(ksp)
its   = KSPGetIterationNumber(ksp)

-- Delete objects
KSPDestroy(ksp)
VecDestroy(Uz)
VecDestroy(Fz)
MatDestroy(Kz)
mesh:delete()

print('Residual:',rnorm, ' NumIterations:', its)

te   = os.clock()
print(te-ts)
